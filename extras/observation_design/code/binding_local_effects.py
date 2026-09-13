"""Resolve the small-effect transition under the corrected observation design."""
import csv,json
import numpy as np
from scipy.stats import t
from binding_observation import experiment,halfpoints,ALPHA0,KG_RANGE,OUT,wilson
rng=np.random.default_rng(9142027);rows=[];trials=[];n=8;N=100
null=max(np.mean(np.log(halfpoints(ALPHA0,k)[1:]/.18))for k in KG_RANGE)
for kg in [7.8,11.7,15.6]:
 for alpha in [.18,.19,.20,.22]:
  decisions=[];rm=[];invalid=0
  for trial in range(N):
   batch=[experiment(rng,alpha,kg,'dose_bias_corrected')for _ in range(n)]
   s=np.array([x[0]for x in batch]);valid=all(x[1]for x in batch);se=s.std(ddof=1)/np.sqrt(n)
   reject=bool(valid and (s.mean()-null)/se>t.ppf(.95,n-1));decisions.append(reject);invalid+=not valid;rm.extend(x[2]for x in batch)
   trials.append(dict(KG_mM=kg,alpha=alpha,trial=trial,mean_score=s.mean(),score_se=se,rejected=int(reject),valid_fit=int(valid)))
  ci=wilson(sum(decisions),N);row=dict(KG_mM=kg,alpha=alpha,trials=N,rejection_rate=np.mean(decisions),interval_low=ci[0],interval_high=ci[1],invalid_trials=invalid,mean_fit_rmse=np.mean(rm));rows.append(row);print(json.dumps(row),flush=True)
for name,data in [('binding_local_effects.csv',rows),('binding_local_effects_trials.csv',trials)]:
 with (OUT/name).open('w')as f:
  w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
