"""Shared ladder parameter fitted jointly across the three cascade summaries.
This probes beyond the fixed illustration values; all upstream means stay matched.
"""
import json
import numpy as np
from scipy.optimize import differential_evolution, least_squares
from analyze import OUT,ROOT,SEED,prepare_input,cascade_score
r=json.loads((OUT/'summary.json').read_text());p=np.array(r['reference']['parameters'])
t=np.loadtxt(ROOT/'data/cascade_comparison.csv',delimiter=',',skiprows=1);mids,obs=t[:,1],t[:,2]
rows=[]
for kind,bounds in [('baseline',[(-1.5,1.5)]),('power',[(-1.5,1.5),(np.log(.3),np.log(5))]),('direct',[(-1.5,1.5),(0,np.log(100)),(np.log(.001),np.log(100))]),('gs_ladder',[(-1.5,1.5),(-.8,.8)])]:
    def obj(v):
        try:
            inp=prepare_input(p,v[0],grid=1501)
            pred=cascade_score(inp,mids,kind,v[1:])
            return float(np.sum(np.log(pred/obs)**2))
        except ValueError:return 1e6
    runs=[differential_evolution(obj,bounds,seed=SEED+k,popsize=8,maxiter=65,tol=2e-7,polish=True) for k in [0,1]]
    best=min(runs,key=lambda z:z.fun)
    def residual(v):return np.log(cascade_score(prepare_input(p,v[0],grid=6001),mids,kind,v[1:])/obs)
    polished=least_squares(residual,best.x,bounds=np.array(bounds).T,xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=300,diff_step=1e-4)
    v=polished.x
    inp=prepare_input(p,v[0],grid=12001)
    pred,kap,curves=cascade_score(inp,mids,kind,v[1:],True)
    row={'kind':kind,'shared_parameters_J_then_route':v.tolist(),'predicted':pred.tolist(),'max_abs_relative_error':float(np.max(abs(pred/obs-1))),'objective':float(np.sum(np.log(pred/obs)**2)),'optimizer_objectives':[float(z.fun) for z in runs],'optimizer_success':[bool(z.success) for z in runs],'kappa':kap,'bounds':bounds,'dense_polish_success':bool(polished.success)}
    rows.append(row);print(json.dumps(row),flush=True)
r['joint_shared_ladder']=rows
(OUT/'summary.json').write_text(json.dumps(r,indent=2)+'\n')
