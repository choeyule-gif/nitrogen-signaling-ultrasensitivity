"""Ordered and random substrate-binding realizations of AceK kinetic summaries.

The six source rows are measured parameter summaries, not raw observations.
The ordered model reproduces the published sequential initial-rate law exactly.
The random model approaches that law in a quantified rapid-binding limit; the
finite random network is solved directly, not replaced with its limit.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'


def stationary(q):
    a=q.copy();a[-1]=1.
    b=np.zeros(len(q));b[-1]=1.
    return np.linalg.solve(a,b)


def law(a,b,row):
    ka,kma,kmb=[float(row[key]) for key in ['K_ATP_uM','Km_ATP_uM','Km_IDH_form_uM']]
    return a*b/(ka*kmb+kma*b+kmb*a+a*b)


def network(a,b,row,kind,gamma=1000.):
    ka,kma,kmb=[float(row[key]) for key in ['K_ATP_uM','Km_ATP_uM','Km_IDH_form_uM']]
    # Time is scaled by each assay's kcat. Printed Vmax converts back to rates.
    n=3 if kind=='ordered' else 4
    q=np.zeros((n,n))
    def edge(i,j,v):q[j,i]+=v;q[i,i]-=v
    if kind=='ordered':
        # States E, EA, EAB; B cannot bind the free enzyme.
        edge(0,1,a/kma);edge(1,0,ka/kma)
        edge(1,2,2*b/kmb);edge(2,1,1.)
        edge(2,0,1.)
    else:
        # States E, EA, EB, EAB; thermodynamic consistency in binding cycle.
        kb=ka*kmb/kma
        edge(0,1,gamma*a/ka);edge(1,0,gamma)
        edge(0,2,gamma*b/kb);edge(2,0,gamma)
        edge(1,3,gamma*b/kmb);edge(3,1,gamma)
        edge(2,3,gamma*a/kma);edge(3,2,gamma)
        edge(3,0,1.)
    return q


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((ROOT/'data/miller1996_table1.csv').open()))
    records=[];design=[];odes=[]
    for row in rows:
        ka,kma,kmb=[float(row[key]) for key in ['K_ATP_uM','Km_ATP_uM','Km_IDH_form_uM']]
        kb=ka*kmb/kma
        # Dimensionless coverage relative to the reported kinetic constants.
        aa=np.geomspace(.01,100,25)*kma
        bb=np.geomspace(.01,100,25)*kmb
        for gamma in [1.,10.,100.,1000.]:
            oe=[];re=[];absolute=[];residual=[]
            for a in aa:
                for b in bb:
                    ref=law(a,b,row)
                    qo=network(a,b,row,'ordered');qr=network(a,b,row,'random',gamma)
                    po=stationary(qo);pr=stationary(qr)
                    oe.append(abs(po[-1]-ref))
                    re.append(abs(pr[-1]/ref-1))
                    absolute.append(abs(pr[-1]-ref))
                    residual.append(float(np.max(abs(qr@pr))))
            assert max(oe)<1e-10
            if gamma==1000:assert max(re)<.01
            records.append(dict(protein=row['protein'],activity=row['activity'],
                binding_to_catalysis_ratio=gamma,grid_points=len(aa)*len(bb),
                max_ordered_absolute_error=max(oe),max_random_relative_error=max(re),
                max_random_absolute_error_Vmax=max(absolute),max_stationary_residual=max(residual)))
        for aratio in [0.,.01,.1,1.,10.]:
            a=aratio*kma
            for bratio in [.1,1.,10.]:
                b=bratio*kb
                po=stationary(network(a,b,row,'ordered'))
                pr=stationary(network(a,b,row,'random',1000.))
                ob=po[-1];rb=pr[2:].sum()
                design.append(dict(protein=row['protein'],activity=row['activity'],ATP_uM=a,IDH_form_uM=b,
                    ATP_over_Km=aratio,protein_over_apparent_Kd=bratio,
                    ordered_protein_bound=ob,random_protein_bound=rb,occupancy_contrast=rb-ob,
                    ordered_rate_nmol_mg_min=float(row['Vmax_nmol_mg_min'])*po[-1],
                    random_rate_nmol_mg_min=float(row['Vmax_nmol_mg_min'])*pr[-1]))
        # Independent solution of the linear finite-network initial-value
        # problem. A matrix exponential avoids unnecessary late-time BDF steps.
        for kind in ['ordered','random']:
            q=network(kma,kb,row,kind,1000.)
            y=np.zeros(len(q));y[0]=1.
            ev=np.linalg.eigvals(q)
            gap=min(-v.real for v in ev if v.real<-1e-8)
            endtime=40/gap
            evolved=expm(q*endtime)@y
            error=float(np.max(abs(evolved-stationary(q))))
            assert error<1e-8
            assert abs(evolved.sum()-1)<1e-8 and evolved.min()>-1e-10
            odes.append(dict(protein=row['protein'],activity=row['activity'],model=kind,max_species_error=error,
                             solver='matrix exponential of the finite reaction generator',time_in_inverse_kcat=endtime))
    for filename,rr in [('kinetic_limit.csv',records),('binding_design.csv',design)]:
        with (OUT/filename).open('w') as f:
            writer=csv.DictWriter(f,fieldnames=rr[0]);writer.writeheader();writer.writerows(rr)
    report=dict(passed=True,source_kinetic_rows=6,finite_network_ODE_checks=odes,
                ordered='Exact steady-state realization of Miller Eq. 1 in rescaled catalytic time.',
                random='Finite mass-action random binding, converging to the same law as binding/catalysis ratio grows.',
                evidence='Published kinetic-parameter summaries across WT and two mutants; no raw replicate fit or new experimental validation.',
                biological_limit='These are compatible initial-rate mechanisms, not validated full AceK phosphorylation cycles. Protein binding under matched ATP/ADP and product conditions is a prospective discriminator. Dimer-level binding and double-AceK complexes are audited separately.')
    (OUT/'kinetic_ambiguity_verification.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
