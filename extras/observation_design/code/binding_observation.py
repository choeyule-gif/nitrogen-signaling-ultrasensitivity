"""End-to-end binding assay benchmark; noise levels are sensitivity scenarios.
Raw titration observations -> finite-depletion fit -> estimated-variance test.
The source model's alpha=.17 is a reference hypothesis, not a measured certainty.
"""
from pathlib import Path
import argparse,csv,json
import numpy as np
from scipy.optimize import least_squares,brentq
from scipy.stats import t
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
P=np.array([0,.005,.01,.02,.05,.1,.2,.5,1,2,5.])
G=np.array([0.,15.6,31.2]) # mM; protein and binding scales are micromolar.
ET=.04;ALPHA0=.17;KG_RANGE=(7.8,15.6)

def halfpoints(alpha,kg,kp=.18):
 g=G/kg;return kp*(1+g)/(1+g/alpha)

def occupancy(pt,k,e):
 c=pt+k+e;return 2*pt/(c+np.sqrt(np.maximum(0,c*c-4*e*pt)))

def mixed_occupancy(pt,k,e,weight):
 # Two independently binding enzyme conformations; ligand mass is conserved.
 if pt==0:return 0.
 def b(p):return weight*p/(k+p)+(1-weight)*p/(5*k+p)
 pf=brentq(lambda p:p+e*b(p)-pt,0,pt)
 return b(pf)

def fit_curves(obs,e_measured,pii_doses=None):
 if pii_doses is None:pii_doses=P[None,:]
 def residual(z):
  pred=z[3]+z[4]*occupancy(pii_doses,np.exp(z[:3,None]),e_measured)
  return (pred[:,:,None]-obs).ravel()
 fit=least_squares(residual,np.r_[np.log([.2,.1,.1]),0.,1.],bounds=(np.r_[np.log([1e-4]*3),-.2,.5],np.r_[np.log([10.]*3),.2,1.5]),ftol=1e-8,gtol=1e-8,xtol=1e-8,max_nfev=150)
 good=fit.success and np.all(fit.active_mask==0)
 h=np.exp(fit.x[:3]);score=np.mean(np.log(h[1:]/h[0]))
 return score,h,good,float(np.sqrt(np.mean(fit.fun**2)))

def experiment(rng,alpha,kg,scenario,return_raw=False):
 # One shared stock calibration and one enzyme-total estimate per experiment.
 stock=np.exp(rng.normal(0,.10)) if scenario!='ideal' else 1.
 em=ET*np.exp(rng.normal(0,.20)) if scenario!='ideal' else ET
 # Systematic nonzero-G dose bias is NOT averaged away by repeated experiments.
 bias=-.10 if scenario in ['dose_bias','dose_bias_corrected'] else 0.
 truep=np.tile(P,(3,1))*stock;truep[1:]*=np.exp(bias)
 h=halfpoints(alpha,kg)
 if scenario=='hidden_conformation':
  vals=np.array([[mixed_occupancy(p,h[j],ET,.9 if j==0 else .6)for p in truep[j]]for j in range(3)])
 else:vals=occupancy(truep,h[:,None],ET)
 background=rng.normal(0,.01);gain=np.exp(rng.normal(0,.03))
 mu=background+gain*vals
 # Two technical readings at each dose, absolute and signal-dependent noise.
 obs=mu[:,:,None]+rng.normal(size=(3,len(P),2))*np.sqrt(.025**2+(.02*vals[:,:,None])**2)
 fitp=np.tile(P,(3,1))
 if scenario=='dose_bias_corrected':
  # Independent concentration-standard readings for each nonzero-G solution.
  fitp[1:]*=np.exp(bias+rng.normal(0,.02,(2,1)))
 score,est,good,rmse=fit_curves(obs,em,fitp)
 if return_raw:return dict(alpha=alpha,kg=kg,scenario=scenario,P_nominal_uM=P.tolist(),G_mM=G.tolist(),E_measured_uM=em,observations=obs.tolist(),H_estimated_uM=est.tolist(),score=score,valid_fit=bool(good),fit_rmse=rmse)
 return score,good,rmse,np.log(est[1:]/est[0])

def wilson(k,n):
 z=1.95996398454;p=k/n;den=1+z*z/n;c=(p+z*z/(2*n))/den;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return [float(c-h),float(c+h)]

def run(trials,n):
 rng=np.random.default_rng(9142026)
 # Null envelope independently profiles KG within a disclosed source-anchored range.
 null=max(np.mean(np.log(halfpoints(ALPHA0,k)[1:]/.18))for k in KG_RANGE)
 rows=[];trialrows=[]
 for scenario in ['ideal','calibration_noise','dose_bias','dose_bias_corrected','hidden_conformation']:
  for kg in [7.8,11.7,15.6]:
   for alpha in [.17,.25,.34,.50,.68]:
    reject=naive_reject=invalid=0;residual=[];scores_all=[];correlations=[]
    for trial in range(trials):
     batch=[experiment(rng,alpha,kg,scenario)for i in range(n)]
     valid=all(x[1]for x in batch);invalid+=not valid
     scores=np.array([x[0]for x in batch]);scores_all.append(scores.mean());residual.extend(x[2]for x in batch)
     se=scores.std(ddof=1)/np.sqrt(n)
     # Shared-baseline covariance is retained in each experiment's two-ratio score.
     statistic=(scores.mean()-null)/se if se>0 else -np.inf
     decision=bool(valid and statistic>t.ppf(.95,n-1));reject+=decision
     ratio=np.array([x[3]for x in batch]);cov=np.cov(ratio.T,ddof=1);corr=cov[0,1]/np.sqrt(cov[0,0]*cov[1,1]);correlations.append(corr)
     naive_se=np.sqrt(np.trace(cov)/(4*n));naive_reject+=bool(valid and (scores.mean()-null)/naive_se>t.ppf(.95,n-1))
     trialrows.append(dict(scenario=scenario,KG_mM=kg,alpha=alpha,trial=trial,mean_score=float(scores.mean()),score_se=float(se),naive_score_se=float(naive_se),baseline_correlation=float(corr),rejected=int(decision),valid_fit=int(valid)))
    row=dict(scenario=scenario,KG_mM=kg,alpha=alpha,experiments_per_trial=n,trials=trials,rejections=reject,rejection_rate=reject/trials,naive_independence_rejection_rate=naive_reject/trials,mean_ratio_correlation=float(np.mean(correlations)),interval_low=wilson(reject,trials)[0],interval_high=wilson(reject,trials)[1],invalid_trials=invalid,mean_score=float(np.mean(scores_all)),null_score=float(null),mean_fit_rmse=float(np.mean(residual)))
    rows.append(row);print(json.dumps(row),flush=True)
 with (OUT/'binding_observation_trials.csv').open('w')as f:
  w=csv.DictWriter(f,fieldnames=trialrows[0]);w.writeheader();w.writerows(trialrows)
 with (OUT/'binding_observation.csv').open('w')as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 (OUT/'binding_observation_example.json').write_text(json.dumps(experiment(np.random.default_rng(41),.34,15.6,'calibration_noise',True),indent=2))
 (OUT/'binding_observation_design.json').write_text(json.dumps(dict(seed=9142026,P_uM=P.tolist(),G_mM=G.tolist(),E_uM=ET,alpha_reference=ALPHA0,KG_nuisance_mM=KG_RANGE,technical_replicates=2,independent_experiments=n,monte_carlo_trials=trials,test='one-sided Student t on each independent experiment mean log-ratio score; shared baseline retained',scope='Known source-model reference alpha, independently profiled KG range; measurement/noise settings are scenarios, not estimated experimental uncertainty. Misspecified scenarios test limits.'),indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=100);p.add_argument('--experiments',type=int,default=8);a=p.parse_args();run(a.trials,a.experiments)
