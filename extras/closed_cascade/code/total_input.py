from model import *
from scipy.optimize import brentq
from pathlib import Path
import json,csv
OUT=Path(__file__).resolve().parents[1]/'results'
rows=[]
for scale in [1,10]:
 for totalG in np.geomspace(10,2000,31):
  args=dict(Pt=5*scale,St=10*scale,Dt=.5*scale,Et=.25*scale)
  vals=[]
  for lam,sigma in [(1,1),(4,4)]:
   def eq(G,ret=False):
    m=Cascade(G=G,lam=lam,sigma=sigma,**args);y=m.steady();o=m.obs(y)
    return (m,y,o) if ret else G+o['boundG']-totalG
   lower=max(1e-9,totalG-2*args['Dt']-args['Et'])
   freeG=brentq(eq,lower,totalG,xtol=1e-10)
   m,y,o=eq(freeG,True);vals.append((freeG,o))
  g1,o1=vals[0];g2,o2=vals[1]
  rows.append(dict(pool_scale=scale,total_G_uM=totalG,free_G_reference=g1,free_G_transformed=g2,free_PII_TV=.5*abs(o1['pfree']-o2['pfree']).sum(),total_GS_gap=abs(o1['GStotal']-o2['GStotal'])))
with (OUT/'total_glutamine.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
s=dict(conditions=len(rows),max_free_PII_TV=max(rows,key=lambda r:r['free_PII_TV']),max_free_G_difference=max(abs(r['free_G_reference']-r['free_G_transformed']) for r in rows),scope='Conserved regulatory glutamine, no metabolic consumption; nucleotide pools buffered; steady roots only, no 121-species ODE validation.')
(OUT/'total_input_summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
