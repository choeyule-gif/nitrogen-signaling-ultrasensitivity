"""Independent stoichiometric balances for a finite enzyme-complex realization.

Construct steady states at buffered free G and free-target totals, rather than
claiming these solutions share closed-system totals across G.
"""
import json
import numpy as np
from analyze import ROOT, OUT, mean, dist, enzyme_rates

def main():
    report=json.loads((OUT/'summary.json').read_text());p=np.array(report['independent_binding']['parameters']);r=p[1]
    k1=report['finite_reaction_realization']['binding_K1_mM'];k2=report['finite_reaction_realization']['binding_K2_mM']
    cases=[]
    for g in [.02,.08,.53,1.,10.]:
        for refr in [False,True]:
            a,b,ar,ah,beta,arate,d=enzyme_rates(g,p)
            alpha=arate if refr else ah
            e=np.array([1.,g/k1,g*g/(k1*k2)]);e*=.001/e.sum()
            theta=mean(g,p);st=dist(theta,r if refr else 0).copy()
            if refr: st[3]-=r
            de=np.zeros(3);ds=np.zeros(4);dc=[];cb=[]
            # kon=2*specificity, koff=kcat=1 gives the desired specificity.
            # Forward complexes carry i UMP; reverse complexes carry i+1.
            for j in range(3):
                for i in range(3):
                    for src,dst,spec in [(i,i+1,(3-i)*alpha[j]),(i+1,i,(i+1)*beta[j])]:
                        bound=spec*e[j]*st[src]
                        association=2*spec*e[j]*st[src];dissociation=bound;catalysis=bound
                        de[j]+=-association+dissociation+catalysis
                        ds[src]+=-association+dissociation;ds[dst]+=catalysis
                        dc.append(association-dissociation-catalysis);cb.append(bound)
            for j,k in [(0,k1),(1,k2)]:
                flux=e[j]*g-e[j+1]*k;de[j]-=flux;de[j+1]+=flux
            residual=float(np.max(np.abs(np.r_[de,ds,dc])))
            assert residual<1e-12 and min(e.min(),st.min(),min(cb))>=-1e-13
            cases.append({'G_mM':g,'refractory':refr,'maximum_species_balance_residual':residual,'enzyme_total':float(e.sum()+sum(cb)),'target_total':float(st.sum()+sum(cb)+(r if refr else 0))})
    rng=np.random.default_rng(20260913); random_error=0.
    for _ in range(1000):
        aa=np.exp(rng.uniform(-3,3,3));bb=np.exp(rng.uniform(-3,3,3));ww=np.exp(rng.uniform(-3,3,3))
        rr=float(np.min(aa/(aa+bb))*rng.uniform(.01,.99));ap=(1-rr)*aa-rr*bb
        assert ap.min()>0
        th=float(aa@ww/((aa+bb)@ww));qq=float(ap@ww/((ap+bb)@ww))
        random_error=max(random_error,abs(th-(rr+(1-rr)*qq)))
    assert random_error<1e-12
    out={'scope':'Explicit ligand-binding and enzyme-target complexes; conditional free-input steady states, not a closed-total titration fit. Zero catalytic branches are omitted reactions. Buffered nucleotide cofactors absorbed into effective catalytic rates.','cases':cases,'maximum_residual':max(c['maximum_species_balance_residual'] for c in cases),'random_positive_parameter_cases':1000,'random_equivalence_error':random_error}
    (OUT/'network_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
