"""Conserved AceK-IDH model and a general Bernstein observation enclosure.

Primary network: Dexter & Gunawardena 2013, Table 1. Coefficients: Table 2.
No experimental curves are manufactured here. Source-rate and observation
scenarios are distinguished in the result metadata.
"""
from pathlib import Path
import csv
import json
from math import comb
import numpy as np
from numpy.polynomial import Polynomial as Poly
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, linprog

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
NAMES=['Ioo','Iop','Ipp','EIoo','EIop','EIpp','EIooE','EIopE','EIppE']
ENZ=np.array([0,0,0,1,1,1,2,2,2])
COUNT=np.array([0,1,2,0,1,2,0,1,2])
# Each entry is (source,target,multiplied_by_free_enzyme), k1 through k20.
EDGES=[(0,3,1),(3,0,0),(3,1,0),(1,4,1),(4,1,0),(4,0,0),(4,2,0),
       (2,5,1),(5,2,0),(5,1,0),(4,7,1),(7,4,0),(7,5,0),(7,3,0),
       (5,8,1),(8,5,0),(3,6,1),(6,3,0),(8,4,0),(6,4,0)]


def source_rates():
    k=np.zeros(21)
    k[[1,4,11,17,8,15]]=100.
    k[[2,5,12,18]]=23.
    k[[9,16]]=40.
    k[[3,7,13,20]]=.098
    k[[6,10,14,19]]=.065
    return k


def generator(e,k,scales=None):
    q=np.zeros((9,9))
    for j,(a,b,binding) in enumerate(EDGES,1):
        v=k[j]*(e if binding else 1.)
        q[b,a]+=v;q[a,a]-=v
    if scales is not None:q=q*np.array(scales)[None,:]
    return q


def distribution(e,k,scales=None):
    q=generator(e,k,scales)
    q[-1,:]=1.
    rhs=np.zeros(9);rhs[-1]=1.
    p=np.linalg.solve(q,rhs)
    return p


def coefficients(k):
    a=k[17]*k[20]/(k[18]+k[20])
    b=k[11]*k[14]/(k[12]+k[13]+k[14])
    c=k[11]*k[13]/(k[12]+k[13]+k[14])
    d=k[15]*k[19]/(k[16]+k[19])
    n0,d0=Poly([k[6],b]),Poly([k[3],a])
    n2,d2=Poly([k[7],c]),Poly([k[10],d])
    weights=[(k[2]*n0+k[6]*d0)*d2/k[1],
             (k[3]*n0*d2+k[5]*d0*d2+k[10]*n2*d0)/k[4],
             (k[7]*d2+k[9]*n2)*d0/k[8]]
    return np.array([np.pad(v.coef,(0,3-len(v.coef))) for v in weights])


def bernstein_vertices(coef,lower,upper):
    """Positive-polynomial enclosure; accepts arbitrary state/degree counts."""
    assert 0<=lower<=upper and np.all(coef>=0)
    n=coef.shape[1]-1
    shifted=np.zeros_like(coef)
    for j in range(n+1):
        for l in range(j,n+1):
            shifted[:,j]+=coef[:,l]*comb(l,j)*lower**(l-j)*(upper-lower)**j
    bernstein=np.zeros_like(coef)
    for k in range(n+1):
        for j in range(k+1):
            bernstein[:,k]+=shifted[:,j]*comb(k,j)/comb(n,j)
    nonzero=bernstein.sum(axis=0)>0
    return bernstein[:,nonzero]/bernstein[:,nonzero].sum(axis=0)


def observation_interval(vertices,epsilon,reporter,observations=()):
    """LP with free-vertex masses and bound modification-state masses."""
    n,m=vertices.shape
    mixture=np.c_[vertices,np.eye(n)]
    c=np.asarray(reporter)@mixture
    au=[np.r_[np.zeros(m),np.ones(n)]];bu=[epsilon]
    for weights,low,high in observations:
        row=np.asarray(weights)@mixture
        au.extend([row,-row]);bu.extend([high,-low])
    kwargs=dict(A_ub=au,b_ub=bu,A_eq=np.ones((1,n+m)),b_eq=[1.],bounds=(0,None),method='highs')
    fits=[linprog(sign*c,**kwargs) for sign in [1.,-1.]]
    if not all(f.success for f in fits):
        return None
    return [float(fits[0].fun),float(-fits[1].fun)]


def closed_stationary(total_target,total_enzyme,k,scales=None):
    def balance(e):
        return e+total_target*(ENZ@distribution(e,k,scales))-total_enzyme
    e=brentq(balance,max(1e-15,total_enzyme*1e-12),total_enzyme,xtol=1e-14)
    return e,total_target*distribution(e,k,scales)


def main():
    OUT.mkdir(exist_ok=True,parents=True)
    k=source_rates();coef=coefficients(k)
    rng=np.random.default_rng(20260916)
    errors=[];column_errors=[];records=[]
    for e in np.geomspace(1e-7,10,101):
        p=distribution(e,k);free=p[:3]/p[:3].sum()
        polynomial=coef@np.array([1,e,e*e]);polynomial/=polynomial.sum()
        errors.append(float(max(abs(free-polynomial))))
        scales=np.r_[np.ones(3),np.exp(rng.uniform(-2,2,6))]
        scaled=distribution(e,k,scales);scaledfree=scaled[:3]/scaled[:3].sum()
        column_errors.append(float(max(abs(free-scaledfree))))
    assert max(errors)<1e-8 and max(column_errors)<1e-8
    # Source-rate target concentration and enzyme totals below are prospective.
    # Historical activity units/mg are not converted into this micromolar grid.
    for it in [.1,1.,10.]:
        for et in [it*.001,it*.01,it*.1]:
            vertices=bernstein_vertices(coef,0,et)
            generic=observation_interval(vertices,min(1.,et/it),[1,0,0])
            assert generic is not None
            values=[];rhsmax=0.;consmax=0.
            for rep in range(100):
                scales=np.r_[np.ones(3),np.exp(rng.uniform(-2,2,6))]
                ef,y=closed_stationary(it,et,k,scales)
                total0=y[[0,3,6]].sum()/it
                values.append(float(total0))
                rhsmax=max(rhsmax,float(np.max(abs(generator(ef,k,scales)@y))))
                consmax=max(consmax,abs(y.sum()-it),abs(ef+ENZ@y-et))
                assert generic[0]-1e-9<=total0<=generic[1]+1e-9
            records.append(dict(IDH_dimer_uM=it,AceK_uM=et,bound_fraction_cap=et/it,
                            interval_low=generic[0],interval_high=generic[1],width=generic[1]-generic[0],
                            sampled_low=min(values),sampled_high=max(values),
                            max_rhs=rhsmax,max_conservation_error=consmax))
    # Prospective matched phosphate-count gate. These limits are declared
    # design thresholds, not measured data or fitted confidence intervals.
    gates=[]
    for et in [.02,.05,.1]:
        it=1.;vertices=bernstein_vertices(coef,0,et)
        interval=observation_interval(vertices,et,[1,0,0],[(np.arange(3),1.435,1.475)])
        generic=observation_interval(vertices,et,[1,0,0])
        assert interval is not None
        accepted=0;max_violation=0.
        for rep in range(1000):
            scales=np.r_[np.ones(3),np.exp(rng.uniform(-2,2,6))]
            ef,y=closed_stationary(it,et,k,scales)
            m=COUNT@y
            if 1.435<=m<=1.475:
                accepted+=1
                val=y[[0,3,6]].sum()
                max_violation=max(max_violation,interval[0]-val,val-interval[1])
        assert accepted>0 and max_violation<1e-9
        gates.append(dict(IDH_dimer_uM=it,AceK_uM=et,phosphate_count_lower=1.435,
                          phosphate_count_upper=1.475,generic_interval=generic,
                          conditioned_interval=interval,generic_width=generic[1]-generic[0],
                          conditioned_width=interval[1]-interval[0],compatible_samples=accepted,
                          total_samples=1000,max_violation=max_violation))
    # Independent conserved mass-action integration from unmodified IDH.
    ode=[]
    for it,et in [(1.,.001),(1.,.1)]:
        ef,ref=closed_stationary(it,et,k)
        y0=np.zeros(10);y0[0]=et;y0[1]=it
        def rhs(t,y):
            dy=generator(y[0],k)@y[1:]
            return np.r_[-ENZ@dy,dy]
        sol=solve_ivp(rhs,(0,2e6),y0,method='BDF',rtol=1e-9,atol=1e-12)
        assert sol.success
        error=float(max(abs(sol.y[:,-1]-np.r_[ef,ref])))
        assert error<1e-7
        ode.append(dict(IDH_dimer_uM=it,AceK_uM=et,max_species_error=error,
                        enzyme_conservation_error=float(max(abs(sol.y[0]+ENZ@sol.y[1:]-et)))))
    # Published rounded data: summarize the overexpression ratio without
    # fitting its in vivo activity to unrelated in vitro concentration units.
    data=list(csv.DictReader((ROOT/'data/laporte1985_table2.csv').open()))
    pair=[row for row in data if row['plasmid'] in ['pBR322','pCK505']]
    overexpression=float(pair[1]['total_activity_units_mg'])/float(pair[0]['total_activity_units_mg'])
    active_ratio=float(pair[1]['activity_during_acetate_growth_units_mg'])/float(pair[0]['activity_during_acetate_growth_units_mg'])
    result=dict(passed=True,free_coefficient_matrix=coef.tolist(),
                normalized_columns=(coef/coef.sum(axis=0)).tolist(),
                coefficient_rank=int(np.linalg.matrix_rank(coef)),
                polynomial_check_max=max(errors),column_scaling_max=max(column_errors),
                checked_free_enzyme_points=101,closed_stationary_pairs=900,ode_checks=ode,
                observation_gate_checks=gates,
                published_overexpression_fold=overexpression,published_active_activity_ratio=active_ratio,
                source='Dexter & Gunawardena 2013 Tables 1-2, with published kinetic values; LaPorte 1985 Table II',
                scope='Mathematical and numerical validation of a source-parameterized model. Column scalings are a constructed dwell-time family, not a measured kinetic confidence set. In vivo table retained as an independent observation constraint, not yet a quantitative matched-assay validation.')
    (OUT/'verification.json').write_text(json.dumps(result,indent=2))
    with (OUT/'concentration_enclosures.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
