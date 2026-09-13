"""Independent identities, leave-one-coordinate-out checks, and design scenarios."""
import json,itertools,math
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit
from analyze import ROOT,OUT,SEED,Ladder,hill,fit,finite_switch,distribution,prepare_input,cascade_score,savecsv

report=json.loads((OUT/'summary.json').read_text());p=np.array(report['reference']['parameters'])
data=np.loadtxt(ROOT/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1);x,y=data.T
checks={};rng=np.random.default_rng(SEED+1)
errors=[]
# Algebraic population equivalence, including downstream nonlinearity, arbitrary n.
for n in range(1,13):
    for _ in range(100):
        r=rng.uniform(.01,.4);qa=rng.uniform(.001,.999);theta=r+(1-r)*qa
        p0b=(1-theta)**n;p0r=(1-r)*(1-qa)**n;k=np.exp(rng.uniform(-8,8))
        kr=k*(1-r)**(n-1)
        errors.extend([abs(p0r/p0b*(1-r)**(n-1)-1),abs(expit(np.log(k*p0b/theta))-expit(np.log(kr*p0r/theta)))])
checks['max_population_equivalence_error']=max(errors)
assert max(errors)<1e-10
# Inversion checked against a bracketed scalar root, not interpolator self-consistency.
err=[];fact=[];ratio=[]
for J in [-.7,0,.7]:
    lad=Ladder(3,J)
    for q in np.linspace(.02,.99,101):
        z=lad.exact_inverse(q);pr=lad.probs(z);mu=pr@np.arange(4);var=pr@(np.arange(4)-mu)**2
        err.append(abs(z-lad.inverse(q)))
        dz=1e-5;der=(lad.mean(z+dz)-lad.mean(z-dz))/(2*dz)
        fact.append(abs(der-var/3))
        ratio.append(abs(3*pr[0]*pr[2]/pr[1]**2-np.exp(2*J)))
checks.update(inverse_max_logu_error=max(err),variance_derivative_max_error=max(fact),odds_curvature_max_error=max(ratio))
assert max(err)<5e-8 and max(fact)<1e-8 and max(ratio)<1e-10

loo=[]
for model in report['finite_switch']:
    J=model['J'];lad=Ladder(3,J);res=[]
    for k in range(len(x)):
        ok=np.arange(len(x))!=k
        pp,ss,success=fit(x[ok],y[ok],lad,start=model['parameters'],multistart=3)
        assert success
        pred=float(finite_switch(x[k],pp,lad));res.append(pred-y[k]);loo.append([J,k,x[k],y[k],pred])
    model['LOO_RMSE']=float(np.sqrt(np.mean(np.array(res)**2)))
    model['point_regulatory_swing']=float(np.exp(lad.exact_inverse(model['parameters'][0]/3)-lad.exact_inverse(model['parameters'][1]/3)))
savecsv('leave_one_coordinate_out.csv','J,index,G_mM,U,predicted_U',loo)

# Maximin state-distribution separation, fixed candidates and finite measured input range.
grid=np.geomspace(.02,10,401);names=list(report['matched'])
ps=np.array([distribution(grid,p,m['J'],m['refractory']) for m in report['matched'].values()])
pairs=list(itertools.combinations(range(4),2))
tv=np.array([.5*np.sum(abs(ps[i]-ps[j]),axis=1) for i,j in pairs])
score=tv.min(axis=0);ix=int(np.argmax(score));gstar=float(grid[ix])
savecsv('state_design.csv','G_mM,min_pair_TV,'+','.join(f'TV_{i}_{j}' for i,j in pairs),np.column_stack([grid,score,tv.T]))
report['design']={'criterion':'maximize minimum pairwise total variation; candidate-specific, no measurement model','range_mM':[.02,10],'grid_points':401,'best_G_mM':gstar,'minimum_pair_TV':float(score[ix]),'pair_order':pairs}

# Prospective measurement precision is a scenario, not assumed historical data noise.
# Propagate fitted-coordinate uncertainty to the two selected assay inputs.
draws=np.loadtxt(OUT/'coordinate_draws.csv',delimiter=',',skiprows=1)
rows=[];intervals={}
for g in [gstar,10.]:
    for i,(name,m) in enumerate(report['matched'].items()):
        pr=np.array([distribution(g,pp,m['J'],m['refractory']) for pp in draws])
        for k in range(len(draws)):rows.append([g,k,i,*pr[k]])
        intervals[f'{name}_{g:.6g}']=np.percentile(pr,[2.5,50,97.5],axis=0).tolist()
report['state_prediction_percentiles']=intervals
savecsv('state_prediction_draws.csv','G_mM,draw,model,p0,p1,p2,p3',rows)

# Independent dense-grid and exact-threshold calculation for the binomial baseline.
tab=np.loadtxt(ROOT/'data/cascade_comparison.csv',delimiter=',',skiprows=1);pii,mids,obs=tab[:,:3].T
def exact_binomial(mid,ref=False):
    r=p[1]/3 if ref else 0
    def b(G):
        q=float(hill(G,p)/3)
        return (1-q)**3/((1-r)**2*q)
    a,c,t=b(0),b(1e100),b(mid)
    k=(a+c-2*t)/(t*(a+c)-2*a*c)
    def f(G):return k*b(G)/(1+k*b(G))
    low,high=f(0),f(1e100)
    roots=[brentq(lambda lg:(f(np.exp(lg))-low)/(high-low)-v,-20,20,xtol=1e-13) for v in [.1,.9]]
    return np.log(81)/(roots[1]-roots[0])
exact=np.array([exact_binomial(m) for m in mids]);dense=cascade_score(prepare_input(p,grid=20001),mids)
checks['baseline_exact_coefficients']=exact.tolist();checks['dense_grid_max_coefficient_error']=float(np.max(abs(exact-dense)))
checks['mixture_exact_cascade_error']=float(max(abs(exact_binomial(m)-exact_binomial(m,True)) for m in mids))
assert checks['dense_grid_max_coefficient_error']<2e-5 and checks['mixture_exact_cascade_error']<1e-10

# Report tolerances explicitly; summary errors are not estimated from missing replicates.
# Vary each midpoint independently by +/-10%; reported coefficients by +/-10%.
sen=[]
for name,m in report['matched'].items():
    inp=prepare_input(p,m['J'],m['refractory'],grid=6001)
    for j,mid in enumerate(mids):
        preds=[cascade_score(inp,[mid*f])[0] for f in [.9,1,1.1]]
        sen.append([names.index(name),pii[j],min(preds),max(preds),obs[j]*.9,obs[j]*1.1])
savecsv('summary_error_scenarios.csv','model,PII_uM,pred_min,pred_max,reported_minus10pct,reported_plus10pct',sen)
report['checks']=checks
report['generalized_plateau']={'independent_conditional_profile_endpoint':173,'n3_swing_endpoint':173**(1/3),'UR_ratio':7.93/1.17,'n3_forward_endpoint':173**(1/3)/(7.93/1.17),'scope':'same plateau criterion, homogeneous positive three-state-count ladder, regulatory floor; not mixture or replicate CI'}
(OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
(OUT/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps({'checks':checks,'design':report['design'],'finite_switch':report['finite_switch'],'generalized_plateau':report['generalized_plateau']},indent=2))
