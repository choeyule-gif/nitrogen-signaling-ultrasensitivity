"""Decision benchmark for A11, retaining source-model assumptions explicitly.
No interval is a confidence region and no numerical example is native validation.
"""
from pathlib import Path
import sys,csv,json
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'closed_cascade/code'))
from model import Cascade
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)

def bounds(G,eps,J):
    Pt,St,Et=5.,10.,.01 # Fixed GS sequestration allowance epsilon_S=0.001.
    Dt=eps*Pt-2*Et
    if Dt<=0:raise ValueError('negative GlnD total')
    m=Cascade(G=G,J=J,Pt=Pt,St=St,Dt=Dt,Et=Et)
    low=m.obs(m.construct(Pt*(1-eps),1.))['GSfree']
    high=m.obs(m.construct(Pt,1.))['GSfree']
    nu=1 if J==0 else 12
    broad=min(1.,np.tanh(nu*(-np.log1p(-eps))/4)+Et/St)
    conditional=min(1.,abs(high-low)+Et/St)
    return m,broad,conditional

def main():
    doses=np.unique(np.r_[np.geomspace(20,10000,41),np.linspace(450,650,21)])
    epsvals=np.geomspace(.00401,.8,25);rows=[]
    for J in [0.,.4]:
      for eps in epsvals:
       for G in doses:
        m,univ,cond=bounds(G,eps,J)
        y=m.steady();o=m.obs(y)
        n=Cascade(G=G,J=J,Pt=5,St=10,Dt=m.Dt,Et=.01,lam=4,sigma=4)
        z=n.steady();v=n.obs(z);gap=abs(o['GStotal']-v['GStotal'])
        assert gap<=cond+1e-10 and cond<=univ+1e-10
        rows.append(dict(G_uM=float(G),epsilon_P=float(eps),epsilon_S=.001,J=J,actual_gap=gap,universal_bound=float(univ),condition_bound=float(cond),condition_certifies=cond<=.01,universal_certifies=univ<=.01))
    with (OUT/'bound_decision.csv').open('w') as f:
      w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    report={'scope':'Specified source topology, fixed Et/St=.001; condition bound uses disclosed regulatory coefficients and free PII composition. Not a measured uncertainty region.','rows':len(rows),'decisions':[]}
    for J in [0.,.4]:
      rr=[r for r in rows if r['J']==J]
      chosen=[r for r in rr if r['condition_certifies']and not r['universal_certifies']]
      grouped=[]
      for e in epsvals:
       sub=[r for r in rr if r['epsilon_P']==e]
       grouped.append(dict(epsilon_P=float(e),max_universal=max(r['universal_bound']for r in sub),max_condition=max(r['condition_bound']for r in sub),max_actual=max(r['actual_gap']for r in sub)))
      report['decisions'].append(dict(J=J,newly_certified_conditions=len(chosen),all_dose_condition_limit=max([r['epsilon_P']for r in grouped if r['max_condition']<=.01],default=None),all_dose_universal_limit=max([r['epsilon_P']for r in grouped if r['max_universal']<=.01],default=None),sweep=grouped))
    (OUT/'bound_decision_summary.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({**report,'decisions':[{k:v for k,v in d.items()if k!='sweep'}for d in report['decisions']]},indent=2))
if __name__=='__main__':main()
