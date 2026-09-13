"""Independent equation, estimation and saved-trial checks for revision 4."""
from pathlib import Path
import csv,json,importlib.util
import numpy as np
from scipy.optimize import brentq
from scipy.stats import t
from binding_observation import occupancy,halfpoints,fit_curves,P,ET,KG_RANGE,ALPHA0
from observation_certificate import certificate
from model_source import Cascade
R=Path(__file__).resolve().parents[1];out=R/'results';rng=np.random.default_rng(23841);checks={}
# Check stable quadratic against independently solved free-ligand mass balance.
worst=0.
for _ in range(500):
 p,k,e=10**rng.uniform(-3,2,3);pf=brentq(lambda x:x+e*x/(k+x)-p,0,p,xtol=1e-14);worst=max(worst,abs(occupancy(p,k,e)-pf/(k+pf)));assert abs(occupancy(p,k,e)-pf/(k+pf))<1e-9
 assert abs(occupancy(k+e/2,k,e)-.5)<1e-10
checks['quadratic_vs_mass_balance_max_error']=worst
# Refit exact observations from a non-starting parameter choice with finite depletion.
h=halfpoints(.23,10.);mu=.035+1.07*occupancy(P[None,:],h[:,None],ET);obs=np.repeat(mu[:,:,None],2,axis=2);score,hh,valid,rmse=fit_curves(obs,ET);assert valid and np.max(abs(hh-h))<1e-7;checks['exact_curve_recovery_max_error']=float(max(abs(hh-h)))
# The null is monotone in both alpha and KG over its declared region.
null=np.mean(np.log(halfpoints(.17,15.6)[1:]/.18))
for alpha in np.geomspace(.001,.17,60):
 for kg in np.linspace(*KG_RANGE,50):assert np.mean(np.log(halfpoints(alpha,kg)[1:]/.18))<=null+1e-12
checks['null_envelope_points']=3000
# Verify saved aggregates from all trial summaries and reconstruct covariance error formula.
a=list(csv.DictReader((out/'binding_observation.csv').open()));tr=list(csv.DictReader((out/'binding_observation_trials.csv').open()));assert len(a)==75 and len(tr)==7500
for row in a:
 d=[x for x in tr if all(x[k]==row[k]for k in ['scenario','KG_mM','alpha'])];assert len(d)==100;reject=sum(int(x['rejected'])for x in d);assert reject==int(row['rejections']);assert reject/100==float(row['rejection_rate'])
 for x in d:
  rej=int(x['valid_fit'])and (float(x['mean_score'])-null)/float(x['score_se'])>t.ppf(.95,7);assert bool(rej)==bool(int(x['rejected']))
for name in ['binding_local_effects']:
 ag=list(csv.DictReader((out/f'{name}.csv').open()));tt=list(csv.DictReader((out/f'{name}_trials.csv').open()));assert len(tt)==1200
 for row in ag:
  d=[x for x in tt if all(x[k]==row[k]for k in ['KG_mM','alpha'])];assert len(d)==100 and sum(int(x['rejected'])for x in d)/100==float(row['rejection_rate'])
checks['independently_recounted_trials']=len(tr)+len(tt)
for _ in range(500):
 ratios=rng.normal(size=(8,2));ratios[:,1]+=.8*ratios[:,0];s=ratios.mean(axis=1).var(ddof=1);c=np.cov(ratios.T);assert abs(s-c.sum()/4)<1e-12
checks['shared_baseline_covariance_checks']=500
# Recompute selected certificates and gate sensitivity without relying on stored bound values.
q=json.loads((out/'observation_certificate.json').read_text());cfg=q['parameters'];r=certificate(**cfg);assert r['bound']<.05<r['generic_bound'];assert certificate(**(cfg|{'KP':180}))['bound']>.05
L=cfg['Pt']*(1-r['epsilon_P']);U=cfg['Pt'];ell=q['ell_limit_for_five_percent'];direct=cfg['Et']/cfg['St']+np.tanh(3*np.log((1+ell/L)/(1+ell/U)));assert abs(direct-.05)<1e-12
checks['coefficient_threshold_recomputed']=direct
# Source-parameter extension must preserve baseline model defaults exactly.
old=R.parent/'closed_cascade/code/model.py'
if not old.exists():old=R/'baseline/S1_Data/extras/closed_cascade/code/model.py'
spec=importlib.util.spec_from_file_location('previous_cascade',old);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
for params in [{},{'lam':4,'sigma':2,'J':.2,'G':2000}]:
 a=Cascade(**params);b=mod.Cascade(**params);ya=a.steady();yb=b.steady();assert np.allclose(ya,yb,rtol=1e-12,atol=1e-12);assert np.max(abs(a.rhs(0,ya)-b.rhs(0,ya)))<1e-12
checks['baseline_default_model_cases']=2
# A fresh source-constrained pair independently checks conservation and a count-based certificate.
a=Cascade(G=2000,Pt=5,St=.5,Dt=.5,Et=.02,beta=175/7.5,beta_prime=1315/7.5,gain=.3);b=Cascade(G=2000,Pt=5,St=.5,Dt=.5,Et=.02,beta=175/7.5,beta_prime=1315/7.5,gain=.3,lam=4);ya=a.steady();yb=b.steady();oa=a.obs(ya);ob=b.obs(yb);M=max(oa['ptotal']@np.arange(4),ob['ptotal']@np.arange(4));v=certificate(**(cfg|{'m_upper':M}));assert abs(oa['GStotal']-ob['GStotal'])<=v['bound'];assert max(abs(a.rhs(0,ya)))<1e-9 and max(abs(b.rhs(0,yb)))<1e-9
checks['fresh_full_reaction_pair_passed']=True
report=dict(passed=True,checks=checks);(out/'revision_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
