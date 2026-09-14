"""Independent numerical checks for the accelerated raw-assay fitter."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares, brentq
from bootstrap_binding import Setting, experiment, fit_curves, occupancy, halfpoints, NULL_SCORE, trial


def main():
    rng=np.random.default_rng(7319)
    largest_derivative_error=0.
    largest_fit_error=0.
    for _ in range(30):
        p=np.exp(rng.uniform(-7,2,(3,11)))
        k=np.exp(rng.uniform(-6,0,(3,1)))
        e=np.exp(rng.uniform(-6,-1))
        v,d=occupancy(p,k,e)
        delta=1e-5
        fd=(occupancy(p,k*np.exp(delta),e)[0]-occupancy(p,k*np.exp(-delta),e)[0])/(2*delta)
        largest_derivative_error=max(largest_derivative_error,float(np.max(abs(fd-d))))
    assert largest_derivative_error<1e-8
    for i in range(50):
        affine=i>=25
        s=Setting(alpha=rng.choice([.17,.19,.25,.34]),kg=rng.uniform(7.8,15.6),
                  mixture=rng.choice([0.,.3,.8]))
        raw=experiment(rng,s,raw=True)
        obs=np.array(raw['observations']);p=np.array(raw['doses_for_fit']);e=raw['enzyme_measured']
        accelerated=fit_curves(obs,e,p,affine)
        def residual(z):
            k=np.exp(z[:3]);pred=np.zeros_like(p)
            for j in range(3):
                for l,pt in enumerate(p[j]):
                    pf=0. if pt==0 else brentq(lambda f:f+e*f/(k[j]+f)-pt,0,pt,xtol=1e-14)
                    bg,gain=(z[3+j],z[6+j]) if affine else (z[3],z[4])
                    pred[j,l]=bg+gain*pf/(k[j]+pf)
            return (pred[:,:,None]-obs).ravel()
        start=np.r_[np.log([.1,.05,.05]),[.01]*3,[1.02]*3] if affine else np.r_[np.log([.1,.05,.05]),.01,1.02]
        independent=least_squares(residual,start,
                                  method='lm',ftol=1e-12,gtol=1e-12,xtol=1e-12,max_nfev=1000)
        assert independent.success
        score=float(np.mean(independent.x[1:3]-independent.x[0]))
        largest_fit_error=max(largest_fit_error,abs(score-accelerated[0]))
    assert largest_fit_error<1e-4
    # Interior-grid profiling is an independent numerical check of the null
    # endpoint used to center the score; it is not a nuisance-grid size proof.
    profile=np.array([np.mean(np.log(halfpoints(.17,kg)[1:]/.18)) for kg in np.linspace(7.8,15.6,501)])
    assert abs(profile.max()-NULL_SCORE)<1e-12 and np.all(np.diff(profile)>0)
    # Recompute every stored audit decision and replay selected complete trials
    # from their seed, including nonlinear fitting rather than score resampling.
    audits=[]
    for folder in ['results','results_condition_affine']:
        base=Path(__file__).resolve().parents[1]/folder
        if not (base/'calibration.json').exists():continue
        import csv
        from scipy.stats import binom
        calibration=json.loads((base/'calibration.json').read_text())
        files=sorted((base/'calibration').glob('cell_*.npz'))
        n=calibration['calibration_trials']
        rank=int(binom.ppf(1-.05/len(files),n,.95))+1
        cutoff=max(float(np.sort(np.load(f)['rows'][:,0])[rank-1]) for f in files)
        assert abs(cutoff-calibration['cutoff'])<1e-12
        summaries=list(csv.DictReader((base/'summary.csv').open()))
        count=0;active=0;failed=0;replay_error=0.
        for stage in ['evaluation','mixture','systematic','signal_bias']:
            rows_for_stage=[r for r in summaries if r['stage']==stage]
            for index,summary in enumerate(rows_for_stage):
                saved=np.load(base/stage/f'cell_{index:03d}.npz')
                rr=saved['rows'];count+=len(rr)
                assert int(np.sum(rr[:,0]>cutoff))==int(summary['rejected'])
                assert int(rr[:,3].sum())==int(summary['optimizer_failed_fits'])
                assert int(rr[:,4].sum())==int(summary['active_bound_fits'])
                active+=int(rr[:,4].sum());failed+=int(rr[:,3].sum())
                assert np.isfinite(rr[:,1:]).all()
                if index in [0,len(rows_for_stage)-1]:
                    s=Setting(**json.loads(str(saved['config'])))
                    for tid in [0,17]:
                        replay=np.array(trial(np.random.SeedSequence([int(saved['seed']),int(saved['phase']),index,tid]),s))
                        error=float(np.max(abs(replay[1:]-rr[tid,1:])))
                        replay_error=max(replay_error,error)
                        assert error<1e-8
        audits.append(dict(folder=folder,independently_recounted_trials=count,
                           failed_fits=failed,active_bound_fits=active,max_trial_replay_error=replay_error))
    result=dict(passed=True,derivative_cases=30,independent_raw_refits=50,stored_trial_audits=audits,
                maximum_derivative_error=largest_derivative_error,
                maximum_score_difference=largest_fit_error,null_profile_points=501)
    out=Path(__file__).resolve().parents[1]/'verification.json'
    out.write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':main()
