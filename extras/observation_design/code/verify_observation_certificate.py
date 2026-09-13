from pathlib import Path
import json,sys
import numpy as np
from model_source import Cascade
from observation_certificate import certificate
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'
rng=np.random.default_rng(8201);rows=[];maxerr=0.
for G in [1300.,2000.,10000.]:
 for kp in [.18,1.]:
  for gain in [.001,.01,.1,1,10]:
   for et in [.1,.04,.02,.01]:
    cfg=dict(G=G,Pt=5.,St=.5,Dt=.5,Et=et,KP=kp,beta=175/7.5,beta_prime=1315/7.5,gain=gain)
    m1=Cascade(**cfg);m2=Cascade(**cfg,lam=4,sigma=1)
    y,z=m1.steady(),m2.steady();a,b=m1.obs(y),m2.obs(z)
    mean_upper=max(float(a['ptotal']@np.arange(4)),float(b['ptotal']@np.arange(4)))
    proof=certificate(Pt=5,St=.5,Dt=.5,Et=et,m_upper=mean_upper,G=G/1000,KG=15.6,KP=kp,alpha=.17,beta=175/7.5,beta_prime=1315/7.5)
    gap=abs(a['GStotal']-b['GStotal']);assert gap<=proof['bound']+1e-9
    maxerr=max(maxerr,float(np.max(abs(m1.rhs(0,y)))),float(np.max(abs(m2.rhs(0,z)))))
    rows.append(dict(G_uM=G,KP_uM=kp,gain=gain,Et_uM=et,PII_count_upper=mean_upper,actual_gap=gap,observation_bound=proof['bound'],generic_bound=proof['generic_bound']))
assert maxerr<1e-8
(OUT/'certificate_full_reaction_checks.json').write_text(json.dumps(dict(cases=len(rows),max_RHS_residual=maxerr,rows=rows,scope='Finite conserved source-topology realizations; the observation bound does not identify a native GlnD mechanism. GlnE sigma fixed at one to retain the source reference regulatory coefficients.'),indent=2))
print(len(rows),'paired finite realizations, maximum RHS residual',maxerr)
