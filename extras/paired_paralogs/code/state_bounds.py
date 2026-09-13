"""Partial identification from independent peptide-level observations.
No independent-site assumption; each trimer has three sites. SE intervals are
working t intervals for published means, not reconstructed replicate likelihoods.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.stats import t
from scipy.optimize import linprog
R=Path(__file__).resolve().parents[1]
with (R/'results/observations.csv').open() as f:obs=list(csv.DictReader(f))
rows=[]
for strain,row in sorted({(x['strain'],x['source_row']) for x in obs}):
 d={x['observable']:x for x in obs if x['strain']==strain and x['source_row']==row}
 for prot,mod in [('GlnB','UMP'),('GlnK','UMP'),('GS','AMP')]:
  names=[prot+'_total',prot+'_'+mod,prot+'_unmodified']
  if not all(k in d for k in names):continue
  T,U,V=[float(d[k]['mean_copy_cell']) for k in names]
  if T<=0 or U+V<=0:continue
  nsite=12 if prot=='GS' else 3
  try:
   uncertainty=[(float(d[k]["SE_copy_cell"]),int(float(d[k]["n"]))) for k in names]
  except (ValueError,TypeError):continue
  ranges=[]
  # Simultaneous marginal intervals by Bonferroni; dependence across assays is not assumed absent.
  for k in names:
   x=d[k];m=float(x['mean_copy_cell']);se=float(x['SE_copy_cell']);n=int(float(x['n']));h=t.ppf(1-.05/6,n-1)*se if n>1 else np.inf
   ranges.append((max(0,m-h),max(0,m+h)))
  (tl,th),(ul,uh),(vl,vh)=ranges
  for estimator,q,lo,hi in [('modified_over_total',U/T,ul/th if th else 0,uh/tl if tl else 1),('paired_peptide_fraction',U/(U+V),ul/(ul+vh) if ul+vh else 0,uh/(uh+vl) if uh+vl else 1)]:
   feasible=0<=q<=1
   lo=max(0,lo);hi=min(1,hi)
   if lo>hi:continue
   # Tight moment bounds; independently solved linear programs test formulas.
   idx=np.arange(nsite+1)/nsite
   aub=np.array([idx,-idx]);bub=np.array([hi,-lo]);aeq=np.ones((1,nsite+1))
   res=[]
   for state in [0,nsite]:
    c=np.eye(nsite+1)[state]
    low=linprog(c,A_ub=aub,b_ub=bub,A_eq=aeq,b_eq=[1],bounds=(0,None),method='highs')
    high=linprog(-c,A_ub=aub,b_ub=bub,A_eq=aeq,b_eq=[1],bounds=(0,None),method='highs')
    assert low.success and high.success
    res.extend([low.fun,-high.fun])
   expected=[max(0,1-nsite*hi),1-lo,max(0,nsite*lo-(nsite-1)),hi]
   assert np.max(np.abs(np.array(res)-expected))<1e-9
   rows.append(dict(strain=strain,time_min=d[names[0]]['time_min'],protein=prot,estimator=estimator,theta=q,point_physically_admissible=feasible,theta_lo=lo,theta_hi=hi,p0_lo=res[0],p0_hi=res[1],pn_lo=res[2],pn_hi=res[3],point_pn_min=max(0,nsite*q-(nsite-1)) if feasible else '',point_p0_max=1-q if feasible else '',source_row=row))
with (R/'results/state_bounds.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
summary=dict(records=len(rows),lp_optimizations=4*len(rows),infeasible_point_ratios=sum(not x['point_physically_admissible'] for x in rows),majority_fully_modified_certified=[{k:x[k] for k in ['strain','time_min','protein','estimator','pn_lo']} for x in rows if x['pn_lo']>.5],interpretation='Bounds on latent oligomer fractions require the stated peptide observation/calibration model. They are not measurements of individual oligomer states. Simultaneous t intervals use source SE/n with no covariance estimate; low n can make intervals uninformative.')
(R/'results/state_bounds_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
