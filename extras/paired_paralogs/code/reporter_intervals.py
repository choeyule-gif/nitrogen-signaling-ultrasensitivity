"""Sharp deterministic intervals for a calibrated paired-paralog assay.
Intervals describe admissible input sets, not automatic confidence intervals.
"""
import json
from pathlib import Path
import numpy as np

def odds(q):
    return q/(1-q)

def calibration(b, k):
    if not (0 < b[0] <= b[1] < 1 and 0 < k[0] <= k[1] < 1):
        raise ValueError('Calibration requires interior fraction intervals')
    return odds(k[0])/odds(b[1]), odds(k[1])/odds(b[0])

def sequestration(b, k, rho, free_error=0):
    """Allow |qK_free - F_rho(qB)| <= free_error as an explicit assumption.
    This term can represent bounded departures from the calibrated stationary
    relation, but cannot be inferred from sampling error or set to zero by default.
    """
    if not (0 <= b[0] <= b[1] <= 1 and 0 <= k[0] <= k[1] <= 1
            and 0 < rho[0] <= rho[1] and 0 <= free_error <= 1):
        raise ValueError('Invalid parameter interval')
    def f(q,r): return r*q/(1-q+r*q)
    lo=max(0.,f(b[0],rho[0])-free_error)
    hi=min(1.,f(b[1],rho[1])+free_error)
    if k[0] > hi: return None
    if hi == 0: return (0.,1.) # necessarily qK_total=0; all s feasible
    lower=0. if lo <= k[1] else 1-k[1]/lo
    upper=1-k[0]/hi
    return float(lower),float(upper)

def verify():
    rng=np.random.default_rng(1537)
    worst=0.
    for _ in range(10000):
        b=sorted(rng.uniform(.001,.999,2)); r=sorted(10**rng.uniform(-2,2,2))
        eps=rng.uniform(0,.15); q=rng.uniform(*b); rr=rng.uniform(*r)
        free=np.clip(rr*q/(1-q+rr*q)+rng.uniform(-eps,eps),0,1)
        s=rng.uniform(); total=(1-s)*free
        k=(max(0,total-rng.uniform(0,.1)),min(1,total+rng.uniform(0,.1)))
        bounds=sequestration(b,k,r,eps)
        assert bounds is not None
        worst=max(worst,bounds[0]-s,s-bounds[1])
        assert bounds[0]-1e-12 <= s <= bounds[1]+1e-12
    assert sequestration((0,0),(0,0),(1,2))==(0.,1.)
    assert sequestration((0,0),(.1,.2),(1,2)) is None
    assert sequestration((.1,.2),(.9,1),(1,1)) is None
    # Unknown specificity gives only the generic unmodified-fraction limit.
    scenarios=[]
    for eps in [0,.03,.1]:
        r=calibration((.48,.52),(.73,.77))
        scenarios.append(dict(free_relation_error=eps,rho_interval=r,
            bound_fraction_interval=sequestration((.48,.52),(.35,.39),r,eps)))
    result=dict(random_feasible_cases=10000,max_interval_violation=worst,
        synthetic_design_scenarios=scenarios,
        interpretation='Deterministic interval propagation. Synthetic values are not measured biology; confidence coverage requires joint input coverage.')
    out=Path(__file__).resolve().parents[1]/'results/reporter_intervals.json'
    out.parent.mkdir(exist_ok=True);out.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__': verify()
