"""Raw-assay parametric null calibration and independent Monte Carlo audit.

All settings are prospective sensitivity scenarios, not measured uncertainties.
Calibration simulates the composite-null grid and repeats nonlinear fitting.
It is a fixed-null parametric bootstrap, not a data-dependent plug-in bootstrap.
The supremum cutoff covers the declared finite grid only, not all nuisance values.
"""
from dataclasses import dataclass, asdict, replace
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import csv
import json
import os
import platform
import time
import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.stats import binom, t

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
P = np.array([0, .005, .01, .02, .05, .1, .2, .5, 1, 2, 5.])
G = np.array([0., 15.6, 31.2])
ET = .04
ALPHA0 = .17
N = 8
LO = np.r_[np.log([1e-4]*3), -.2, .5]
HI = np.r_[np.log([10.]*3), .2, 1.5]


@dataclass(frozen=True)
class Setting:
    kg: float = 15.6
    alpha: float = ALPHA0
    calibration: float = 1.
    signal: float = 1.
    dose_bias: float = 0.
    corrected: bool = False
    mixture: float = 0.
    differential_background: float = 0.
    differential_gain: float = 0.
    condition_affine: bool = False


def halfpoints(alpha, kg):
    g = G / kg
    return .18 * (1 + g) / (1 + g / alpha)


NULL_SCORE = float(np.mean(np.log(halfpoints(ALPHA0, 15.6)[1:] / .18)))


def occupancy(p, k, e):
    c = p + k + e
    root = np.sqrt(c*c - 4*e*p)
    v = 2*p / (c + root)
    return v, -k*v/root


def mixture_occupancy(p, k, e, w):
    # Bisection of free ligand conserves mass for both affinity classes.
    low, high = np.zeros_like(p), p.copy()
    for _ in range(55):
        free = (low + high)/2
        b = w*free/(k+free) + (1-w)*free/(5*k+free)
        above = free + e*b > p
        high = np.where(above, free, high)
        low = np.where(above, low, free)
    return w*free/(k+free) + (1-w)*free/(5*k+free)


def fit_curves(obs, e, doses, condition_affine=False):
    # With equal replicate weights, fitting replicate means is exactly the
    # same nonlinear least-squares problem, up to an additive constant.
    means = obs.mean(axis=2)
    def residual(z):
        v, _ = occupancy(doses, np.exp(z[:3,None]), e)
        bg,gain=(z[3:6,None],z[6:9,None]) if condition_affine else (z[3],z[4])
        return (bg + gain*v - means).ravel()
    def jac(z):
        v, deriv = occupancy(doses, np.exp(z[:3,None]), e)
        jj = np.zeros((3, len(P), 9 if condition_affine else 5))
        for j in range(3):
            jj[j,:,j] = (z[6+j] if condition_affine else z[4])*deriv[j]
            if condition_affine:
                jj[j,:,3+j]=1.
                jj[j,:,6+j]=v[j]
        if not condition_affine:
            jj[:,:,3] = 1.
            jj[:,:,4] = v
        return jj.reshape(len(P)*3, -1)
    start=np.r_[np.log([.2,.1,.1]),np.zeros(3),np.ones(3)] if condition_affine else np.r_[np.log([.2,.1,.1]),0.,1.]
    lower=np.r_[LO[:3],[-.2]*3,[.5]*3] if condition_affine else LO
    upper=np.r_[HI[:3],[.2]*3,[1.5]*3] if condition_affine else HI
    fit = least_squares(residual, start,
                        jac=jac, bounds=(lower,upper), max_nfev=150,
                        ftol=1e-8, xtol=1e-8, gtol=1e-8)
    score = float(np.mean(fit.x[1:3] - fit.x[0]))
    active = bool(np.any(fit.active_mask != 0))
    failed = not bool(fit.success and np.isfinite(fit.fun).all())
    # Separate numerical failure, active bound, and goodness of fit.
    pred = means + fit.fun.reshape(means.shape)
    rmse = float(np.sqrt(np.mean((pred[:,:,None]-obs)**2)))
    return score, failed, active, rmse, int(fit.nfev), fit.x


def experiment(rng, setting, raw=False):
    s = setting
    stock = np.exp(rng.normal(0, .10*s.calibration))
    em = ET*np.exp(rng.normal(0, .20*s.calibration))
    doses = np.tile(P, (3,1))
    truep = doses*stock
    truep[1:] *= np.exp(s.dose_bias)
    h = halfpoints(s.alpha,s.kg)
    if s.mixture:
        # A continuous path starts at the correctly specified one-class model.
        # At m=1, low-affinity weights are .1 at G=0 and .4 at G>0.
        weights = (1-s.mixture*np.array([.1,.4,.4]))[:,None]
        vals = mixture_occupancy(truep,h[:,None],ET,weights)
    else:
        vals = occupancy(truep,h[:,None],ET)[0]
    bg = rng.normal(0,.015*s.signal)
    gain = np.exp(rng.normal(0,.05*s.signal))
    mu = bg + gain*vals
    mu[1:] = bg + s.differential_background + gain*np.exp(s.differential_gain)*vals[1:]
    obs = mu[:,:,None] + rng.normal(size=(3,len(P),2))*np.sqrt(.025**2+(.02*vals[:,:,None])**2)
    if s.corrected:
        doses[1:] *= np.exp(s.dose_bias+rng.normal(0,.02,(2,1)))
    fit = fit_curves(obs,em,doses,s.condition_affine)
    if raw:
        return dict(observations=obs.tolist(), enzyme_measured=em,
                    doses_for_fit=doses.tolist(), setting=asdict(s),
                    fitted_parameters=fit[-1].tolist())
    return fit[:5]


def trial(seed, setting):
    rng = np.random.default_rng(seed)
    batch = [experiment(rng,setting) for _ in range(N)]
    scores = np.array([b[0] for b in batch])
    se = scores.std(ddof=1)/np.sqrt(N)
    failed = sum(b[1] for b in batch)
    active = sum(b[2] for b in batch)
    valid = not (failed or active) and np.isfinite(se) and se > 0
    stat = (scores.mean()-NULL_SCORE)/se if valid else -np.inf
    return (stat, scores.mean(), se, failed, active,
            np.mean([b[3] for b in batch]), np.mean([b[4] for b in batch]))


def cell_job(task):
    phase, index, config, trials, seed = task
    setting = Setting(**config)
    rows = np.array([trial(np.random.SeedSequence([seed,phase,index,i]),setting)
                     for i in range(trials)])
    return index, config, rows


def grid():
    return [asdict(Setting(kg=kg,calibration=c,signal=s))
            for kg in [7.8,11.7,15.6] for c in [0.,.5,1.] for s in [0.,1.,2.]]


def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(mid-half),float(mid+half)]


def run_cells(name,phase,configs,trials,workers,seed):
    folder=OUT/name;folder.mkdir(parents=True,exist_ok=True)
    tasks=[]
    for i,s in enumerate(configs):
        dest=folder/f"cell_{i:03d}.npz"
        if dest.exists():
            saved=np.load(dest)
            oldconfig=asdict(Setting(**json.loads(str(saved['config']))))
            if (oldconfig==s
                and saved['rows'].shape==(trials,7)
                and int(saved['seed'])==seed and int(saved['phase'])==phase):
                continue
            raise RuntimeError(f"Incompatible saved cell: {dest}")
        tasks.append((phase,i,s,trials,seed))
    print(json.dumps(dict(stage=name,cells=len(configs),pending=len(tasks),trials=trials)),flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(cell_job,task) for task in tasks]
        for future in as_completed(futures):
            i,s,rows=future.result()
            np.savez_compressed(folder/f"cell_{i:03d}.npz", rows=rows,
                                config=json.dumps(s,sort_keys=True),seed=seed,phase=phase)
            print(json.dumps(dict(stage=name,cell=i,config=s,failed=int(rows[:,3].sum()),
                                  active=int(rows[:,4].sum()))),flush=True)
    return [(s,np.load(folder/f"cell_{i:03d}.npz")['rows']) for i,s in enumerate(configs)]


def run(args):
    global OUT
    OUT=ROOT/args.output
    OUT.mkdir(parents=True,exist_ok=True)
    start=time.time()
    nullconfigs=grid()
    if args.pilot:
        nullconfigs=[asdict(Setting())]
    nullconfigs=[dict(s,condition_affine=args.condition_affine) for s in nullconfigs]
    calibration=run_cells('calibration',0,nullconfigs,args.calibration_trials,args.workers,args.seed)
    # A one-sided binomial tolerance order statistic gives simultaneous 95%
    # Monte Carlo confidence that every GRID null tail is at most .05.
    n=args.calibration_trials
    rank=int(binom.ppf(1-.05/len(nullconfigs),n,.95))+1
    if rank>n:
        raise ValueError('Too few calibration trials for a finite tolerance cutoff')
    thresholds=[float(np.sort(rows[:,0])[rank-1]) for s,rows in calibration]
    cutoff=max(thresholds)
    record=dict(cutoff=cutoff,empirical_cutoff=max(float(np.quantile(r[:,0],.95)) for s,r in calibration),
                rank_one_based=rank,calibration_trials=n,null_cells=len(nullconfigs),
                per_cell_cutoffs=thresholds,nominal_tail=.05,simultaneous_mc_confidence=.95,
                scope='Finite declared nuisance grid; not a uniform guarantee between grid points or under misspecification.')
    (OUT/'calibration.json').write_text(json.dumps(record,indent=2))
    configs=[dict(s,alpha=a) for s in nullconfigs for a in [ALPHA0,.19,.25,.34]]
    evalsets=[('evaluation',1,configs)]
    if not args.pilot:
        evalsets += [('mixture',2,[asdict(Setting(kg=k,mixture=float(m)))
                     for k in [7.8,11.7,15.6] for m in np.linspace(0,1,11)]),
                     ('systematic',3,[asdict(Setting(kg=k,dose_bias=b,corrected=c))
                     for k in [7.8,11.7,15.6] for b in [-.15,-.10,-.05,0.,.05,.10,.15]
                     for c in [False,True]]),
                     ('signal_bias',4,[asdict(Setting(kg=k,differential_background=b,differential_gain=g))
                     for k in [7.8,11.7,15.6] for b,g in
                     [(-.02,0.),(-.01,0.),(0.,0.),(.01,0.),(.02,0.),(0.,-.1),(0.,-.05),(0.,.05),(0.,.1)]])]
    summaries=[]
    for name,phase,configs in evalsets:
        configs=[dict(s,condition_affine=args.condition_affine) for s in configs]
        evaluation=run_cells(name,phase,configs,args.trials,args.workers,args.seed)
        for s,rows in evaluation:
            rejected=int(np.sum(rows[:,0]>cutoff));n=len(rows)
            naive=int(np.sum(rows[:,0]>t.ppf(.95,N-1)))
            lo,hi=wilson(rejected,n)
            summaries.append(dict(stage=name,**s,trials=n,rejected=rejected,
                rejection_rate=rejected/n,ci_low=lo,ci_high=hi,uncalibrated_t_rate=naive/n,
                optimizer_failed_fits=int(rows[:,3].sum()),active_bound_fits=int(rows[:,4].sum()),
                total_fits=n*N,invalid_trials=int(np.sum(~np.isfinite(rows[:,0]))),
                mean_rmse=float(rows[:,5].mean()),mean_nfev=float(rows[:,6].mean())))
    with (OUT/'summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=summaries[0]);writer.writeheader();writer.writerows(summaries)
    meta=dict(seed=args.seed,experiments_per_trial=N,technical_replicates=2,condition_affine=args.condition_affine,
              null_alpha=ALPHA0,kg_grid=[7.8,11.7,15.6],P_uM=P.tolist(),G_mM=G.tolist(),
              enzyme_uM=ET,noise_absolute=.025,noise_proportional=.02,
              stock_logsd='0.10 * calibration',enzyme_logsd='0.20 * calibration',
              background_sd='0.015 * signal',gain_logsd='0.05 * signal',
              test='Grid-supremum fixed-null parametric bootstrap; each trial includes eight nonlinear raw-curve refits.',
              independent_seeds='SeedSequence([seed,phase,cell,trial]); phase 0 calibration, phases 1-4 audit',
              failed_rule='A failed or active-bound fit invalidates its whole trial; invalid trials do not reject and remain in the denominator.',
              runtime_seconds=time.time()-start,host=platform.node(),python=platform.python_version(),
              numpy=np.__version__,scipy=scipy.__version__,workers=args.workers,
              limitations='Scenario errors are stipulated, not estimated from assay replicates. Finite-grid calibration is not native biology validation.')
    (OUT/'design.json').write_text(json.dumps(meta,indent=2))
    (OUT/'example.json').write_text(json.dumps(experiment(np.random.default_rng(829),Setting(alpha=.25,condition_affine=args.condition_affine),raw=True),indent=2))
    print(json.dumps(dict(done=True,seconds=meta['runtime_seconds'],cutoff=cutoff,cells=len(summaries))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--trials',type=int,default=2000)
    p.add_argument('--calibration-trials',type=int,default=3000)
    p.add_argument('--workers',type=int,default=16)
    p.add_argument('--seed',type=int,default=20260915)
    p.add_argument('--output',default='results')
    p.add_argument('--pilot',action='store_true')
    p.add_argument('--condition-affine',action='store_true')
    run(p.parse_args())
