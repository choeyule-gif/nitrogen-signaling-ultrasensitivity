from model import *
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix
import json,csv
OUT=Path(__file__).resolve().parents[1]/'results';OUT.mkdir(exist_ok=True)
rng=np.random.default_rng(20260914)
rows=[];maxrhs=0.;maxmass=0.;maxfree=0.;viol=0.
for direct in [0.,1.]:
 for J in [0.,.4]:
  for G in np.unique(np.r_[np.geomspace(20,10000,41),np.linspace(450,650,41)]):
   for eps in np.geomspace(.0001,.8,25):
    args=dict(G=G,direct=direct,J=J,Pt=5.,St=10.,Dt=eps*5/2,Et=eps*5/4)
    m1=Cascade(**args);m2=Cascade(**args,lam=4,sigma=4)
    y1=m1.steady();y2=m2.steady();o1=m1.obs(y1);o2=m2.obs(y2)
    gapP=.5*np.abs(o1['ptotal']-o2['ptotal']).sum();gapS=abs(o1['GStotal']-o2['GStotal']);gapSF=abs(o1['GSfree']-o2['GSfree'])
    L=-np.log1p(-eps);nu=1 if J==0 else 12
    bound=min(1,(np.tanh(nu*L/4) if direct else 0)+args['Et']/args['St'])
    # Condition-specific envelope from pool extremes, independent of root locations.
    lo=m1.obs(m1.construct(5*(1-eps),1))['GSfree'];hi=m1.obs(m1.construct(5,1))['GSfree']
    tight=min(1,abs(lo-hi)+args['Et']/args['St'])
    viol=max(viol,gapP-eps,gapS-bound,gapS-tight)
    maxfree=max(maxfree,float(np.max(np.abs(o1['pfree']-o2['pfree']))))
    maxrhs=max(maxrhs,float(np.max(np.abs(m1.rhs(0,y1)))),float(np.max(np.abs(m2.rhs(0,y2)))))
    maxmass=max(maxmass,float(np.max(np.abs(m1.cons@y1/np.array([5,10,args['Dt'],args['Et']])-1))))
    rows.append(dict(G_uM=G,epsilon=eps,direct=direct,J=J,PII_TV=gapP,GS_total_gap=gapS,GS_free_gap=gapSF,GS_bound=bound,GS_condition_bound=tight))
with (OUT/'cascade_sweep.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('sweep',len(rows),'rhs',maxrhs,'bound_violation',viol,flush=True)
# Random conservation/residual and two widely separated starts, varying pool ratios.
randomerr=0.;rooterr=0.
for k in range(60):
 m=Cascade(G=10**rng.uniform(1,4),lam=rng.uniform(1,6),sigma=rng.uniform(1,5),J=rng.uniform(0,.4),Pt=10**rng.uniform(-.3,1),St=10**rng.uniform(0,2),Dt=10**rng.uniform(-3,-.3),Et=10**rng.uniform(-3,-.3))
 y=m.steady();z=m.steady((.02,.02));rooterr=max(rooterr,float(np.max(np.abs(y-z))))
 randomerr=max(randomerr,float(np.max(np.abs(m.rhs(0,y)))))
# Independent integration of full reactions: no reduced steady equations in RHS.
odes=[]
for lam,sigma,J in [(1,1,0),(4,4,0),(4,4,.1)]:
 m=Cascade(G=530,lam=lam,sigma=sigma,J=J,Dt=.5,Et=.25)
 eq=m.steady();y0=np.zeros(120);y0[m.P[0]]=5;y0[m.S[0]]=10;y0[m.D[0]]=.5;y0[m.E[0]]=.25
 dep=np.zeros((len(m.ks),120))
 for k,react in enumerate(m.rs):dep[k,react]=1
 sparsity=csr_matrix((np.abs(m.N)@dep)>0)
 sol=solve_ivp(m.rhs,[0,2e6],y0,method='BDF',rtol=2e-8,atol=2e-11,jac_sparsity=sparsity)
 err=float(np.max(np.abs(sol.y[:,-1]-eq)))
 assert sol.success and err<2e-5
 odes.append(dict(lam=lam,sigma=sigma,J=J,max_concentration_error=err,rhs_residual=float(np.max(np.abs(m.rhs(0,sol.y[:,-1])))),minimum_concentration=float(sol.y.min())))
 print('ode',odes[-1],flush=True)
summary=dict(species=120,reactions=Cascade().N.shape[1],paired_sweep_conditions=len(rows),random_networks=60,max_steady_RHS=maxrhs,max_relative_mass_error=maxmass,max_free_PII_distribution_gap=maxfree,bound_violation=viol,random_RHS=randomerr,multistart_max_difference=rooterr,independent_ODE=odes)
for direct in [0,1]:
 for J in [0,.4]:
  subset=[r for r in rows if r['direct']==direct and r['J']==J]
  best=max(subset,key=lambda r:r['GS_total_gap'])
  summary[f'direct{direct}_J{J}']=best
  summary[f'max_free_direct{direct}_J{J}']=max(subset,key=lambda r:r['GS_free_gap'])
assert maxfree<1e-12 and maxrhs<1e-8 and viol<1e-10 and rooterr<1e-6
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
