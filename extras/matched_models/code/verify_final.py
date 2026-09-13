"""Scientific cross-checks on the complete generated analysis, independent where possible."""
import json
import numpy as np
from scipy.special import logsumexp
from analyze import OUT,SEED,Ladder,distribution,hill,cascade_score,prepare_input
s=json.loads((OUT/'summary.json').read_text());p=np.array(s['reference']['parameters'])
rng=np.random.default_rng(SEED+9);errors=[]
for J in [-.7,0,.7]:
    L=Ladder(3,J)
    for g in np.geomspace(.02,10,41):
        q=float(hill(g,p)/3);u=np.exp(L.exact_inverse(q));c=np.exp(L.logc)
        Q=np.zeros((4,4))
        for i in range(3):
            forward=u*c[i+1]/c[i];reverse=1.
            Q[i+1,i]=forward;Q[i,i]-=forward;Q[i,i+1]=reverse;Q[i+1,i+1]-=reverse
        A=Q.copy();A[-1,:]=1;b=np.array([0.,0,0,1.])
        pk=np.linalg.solve(A,b)
        errors.append(float(np.max(abs(pk-L.probs(np.log(u))))))
assert max(errors)<1e-11
for z in s['finite_switch']:
    assert z['success'] and np.isfinite(z['LOO_RMSE'])
for z in s['joint_shared_ladder']:
    assert all(z['optimizer_success']) and z['dense_polish_success']
    assert all(np.array(z['predicted'])>0)
assert s['asymmetric_gs']['max_curve_equivalence_error']<1e-9
assert len(np.loadtxt(OUT/'coordinate_draws.csv',delimiter=',',skiprows=1))==300
grid=np.geomspace(.02,10,301)
for name,m in s['matched'].items():
    pr=distribution(grid,p,m['J'],m['refractory'])
    assert np.min(pr)>=0 and np.max(abs(pr.sum(axis=1)-1))<1e-12
    assert np.max(abs(pr@np.arange(4)-hill(grid,p)))<1e-7
# Verify nominal upstream-profile entries include both declared endpoint scenarios.
prof=np.loadtxt(OUT/'finite_switch_profile.csv',delimiter=',',skiprows=1)
assert len(prof)==61 and np.all(prof[:,-1]==1)
checks={'passed':True,'birth_death_stationarity_max_error':max(errors),'independent_stationary_solves':len(errors),'coordinate_draws':300,'fixed_cascade_comparisons':len(s['cascade']),'joint_shared_comparisons':len(s['joint_shared_ladder']),'additional_asymmetric_comparisons':1,'notes':'Checks numerical consistency, not molecular validity, experimental fit uncertainty, or a global optimum.'}
(OUT/'final_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
