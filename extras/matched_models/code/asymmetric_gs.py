"""Shared two-parameter asymmetric effective GS ladder; same calibration budget."""
import json
import numpy as np
from scipy.optimize import differential_evolution,least_squares
from analyze import OUT,ROOT,SEED,prepare_input,cascade_score,savecsv
r=json.loads((OUT/'summary.json').read_text());p=np.array(r['reference']['parameters'])
t=np.loadtxt(ROOT/'data/cascade_comparison.csv',delimiter=',',skiprows=1);mids,obs=t[:,1],t[:,2]
inp=prepare_input(p,grid=3001)
bounds=[(-.5,.5),(-.5,.5)]
def obj(v):
    try:return float(np.sum(np.log(cascade_score(inp,mids,'gs_asym',v)/obs)**2))
    except ValueError:return 1e6
fits=[differential_evolution(obj,bounds,seed=SEED+i,popsize=10,maxiter=100,tol=1e-8,polish=True) for i in [0,1]]
best=min(fits,key=lambda z:z.fun)
inp=prepare_input(p,grid=12001)
fit=least_squares(lambda v:np.log(cascade_score(inp,mids,'gs_asym',v)/obs),best.x,bounds=np.array(bounds).T,diff_step=1e-4,xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=400)
pred,k,curves=cascade_score(inp,mids,'gs_asym',fit.x,True)
prefr,kr,curver=cascade_score(prepare_input(p,refractory=True,grid=12001),mids,'gs_asym',fit.x,True)
entry={'kind':'gs_asym','shared_parameters_quadratic_cubic':fit.x.tolist(),'predicted':pred.tolist(),'reported':obs.tolist(),'max_abs_relative_error':float(np.max(abs(pred/obs-1))),'objective':float(np.sum(np.log(pred/obs)**2)),'optimizer_success':[bool(z.success) for z in fits],'optimizer_objectives':[float(z.fun) for z in fits],'polish_success':bool(fit.success),'bounds':bounds,'kappa_independent':k,'kappa_refractory':kr,'max_curve_equivalence_error':float(np.max(abs(curves-curver))),'max_coefficient_equivalence_error':float(np.max(abs(pred-prefr)))}
r['asymmetric_gs']=entry
savecsv('cascade_independent_gs_asym.csv','G_mM,Y_0p5,Y_5,Y_36',np.column_stack([np.exp(inp[0]),curves.T]))
savecsv('cascade_refractory_gs_asym.csv','G_mM,Y_0p5,Y_5,Y_36',np.column_stack([np.exp(inp[0]),curver.T]))
(OUT/'summary.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(entry,indent=2))
