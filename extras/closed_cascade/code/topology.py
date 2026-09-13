from model import *
from pathlib import Path
import csv,json
OUT=Path(__file__).resolve().parents[1]/'results'
rows=[]
for eta in np.r_[0,np.geomspace(.0001,1,20)]:
 for G in np.geomspace(20,10000,41):
  a=Cascade(G=G,alpha2=2.17,productive=eta)
  b=Cascade(G=G,alpha2=8.85,productive=eta)
  # Matched FREE protein pools isolate loss of silent-state compensation.
  oa=a.obs(a.construct(5,10));ob=b.obs(b.construct(5,10))
  # Protein-conserving counterpart, same totals.
  ya=a.steady();yb=b.steady();ca=a.obs(ya);cb=b.obs(yb)
  rows.append(dict(eta=eta,G_uM=G,fixed_free_GS_gap=abs(oa['GSfree']-ob['GSfree']),closed_GS_gap=abs(ca['GStotal']-cb['GStotal']),RHS=max(np.max(np.abs(a.rhs(0,ya))),np.max(np.abs(b.rhs(0,yb))))))
with (OUT/'productive_topology.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
# Recognition-weight uncertainty within the specified class.
maxviolation=0.;maxrhs=0.
for chi in [.5,1,2]:
 for G in np.geomspace(20,10000,15):
  for J in [0,.4]:
   a=Cascade(G=G,J=J,recognition=chi,Dt=.25,Et=.125)
   b=Cascade(G=G,J=J,recognition=chi,Dt=.25,Et=.125,lam=4,sigma=4)
   ya=a.steady();yb=b.steady();oa=a.obs(ya);ob=b.obs(yb)
   bound=np.tanh((1 if J==0 else 12)*(-np.log(.9))/4)+.0125
   maxviolation=max(maxviolation,abs(oa['GStotal']-ob['GStotal'])-bound)
   maxrhs=max(maxrhs,np.max(np.abs(a.rhs(0,ya))),np.max(np.abs(b.rhs(0,yb))))
summary=dict(conditions=len(rows),max_fixed_free_GS_gap=max(rows,key=lambda x:x['fixed_free_GS_gap']),zero_productivity_max_free_gap=max(r['fixed_free_GS_gap'] for r in rows if r['eta']==0),max_RHS=max(r['RHS'] for r in rows),recognition_conditions=90,recognition_bound_violation=maxviolation,recognition_RHS=maxrhs)
(OUT/'topology_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
# Productive EPU does not abolish every closure bound: for fixed alpha2 and
# paired lambda/sigma, r=a+b/T+cT still has |d log r/d log T| <= 1.
viol=0.;maxgap=0.;rhs=0.;n=0;extra=[]
for eta in [.1,1]:
 for direct in [0,1]:
  for J in [0,.4]:
   for eps in [.01,.1,.8]:
    for G in np.geomspace(100,2000,15):
     kw=dict(G=G,direct=direct,J=J,Dt=2.5*eps,Et=1.25*eps,productive=eta)
     a=Cascade(**kw);b=Cascade(**kw,lam=4,sigma=4)
     ya=a.steady();yb=b.steady();oa=a.obs(ya);ob=b.obs(yb)
     bound=min(1,np.tanh((1 if J==0 else 12)*(-np.log1p(-eps))/4)+.125*eps)
     gap=abs(oa['GStotal']-ob['GStotal']);viol=max(viol,gap-bound);maxgap=max(maxgap,gap)
     rhs=max(rhs,np.max(np.abs(a.rhs(0,ya))),np.max(np.abs(b.rhs(0,yb))));n+=1
     extra.append(dict(eta=eta,direct=direct,J=J,epsilon=eps,G_uM=G,total_GS_gap=gap,bound=bound))
summary['productive_compensation_bound']=dict(conditions=n,max_violation=viol,max_total_GS_gap=maxgap,max_RHS=rhs)
(OUT/'topology_summary.json').write_text(json.dumps(summary,indent=2));print('productive closure',summary['productive_compensation_bound'])

with (OUT/'productive_compensation.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=extra[0]);w.writeheader();w.writerows(extra)
assert viol<1e-10 and rhs<1e-8
