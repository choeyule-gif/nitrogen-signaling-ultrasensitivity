from native_trimers import *
report=[];pred=[]
for lam in [1.,2.,4.]:
 def residual(v):
  p,_=simulate('WT',lam,v);return (p[2:,4]-p[2:,2])/p[2:,3]
 starts=[np.log([1,1,1.5,1,6]),np.log([5,5,270,.01,1])]
 fits=[least_squares(residual,v,bounds=([-4,-4,-4,-5,-3],[4,4,6,5,3]),diff_step=.005,ftol=3e-5,xtol=3e-5,gtol=3e-5,max_nfev=120) for v in starts]
 fit=min(fits,key=lambda x:x.fun@x.fun)
 for strain in ['WT','delta_glnK','delta_glnB']:
  z,qa=simulate(strain,lam,fit.x);sc=z[2 if strain=='WT' else 1:];rec=dict(lam=lam,strain=strain,parameters=np.exp(fit.x).tolist(),weighted_SSE=float(np.sum(((sc[:,4]-sc[:,2])/sc[:,3])**2)),RMSE_fraction=float(np.sqrt(np.mean((sc[:,4]-sc[:,2])**2))),optimizer_success=bool(fit.success),**qa);report.append(rec);print(rec,flush=True)
  pred.extend(dict(lam=lam,strain=strain,time_min=t,protein='GlnB' if j==0 else 'GlnK',observed_fraction=q,SE_upper=se,predicted_fraction=p,bound_fraction=s if strain!='delta_glnK' else '') for t,j,q,se,p,s in z)
(R/'results/paralog_profile_summary.json').write_text(json.dumps(report,indent=2))
with (R/'results/paralog_profile_predictions.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=pred[0]);w.writeheader();w.writerows(pred)
