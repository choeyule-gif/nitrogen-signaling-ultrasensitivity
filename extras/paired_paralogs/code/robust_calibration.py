"""Sharp partial identification from bounded free and bound reporter responses.
This is an observation-level calculation, not biological validation.
"""
import numpy as np
from scipy.optimize import linprog,brentq
from scipy.special import logsumexp
import json
from pathlib import Path

def mixture_interval(free,bound,total):
    for interval in [free,bound,total]:
        if not 0<=interval[0]<=interval[1]<=1:raise ValueError('fractions must be ordered in [0,1]')
    lo,hi=0.,1.
    # a*s <= c, from overlap of the attainable mixture and observed interval.
    for a,c in [(bound[0]-free[0],total[1]-free[0]),(free[1]-bound[1],free[1]-total[0])]:
        if a==0:
            if c < 0:return None
        elif a>0:hi=min(hi,c/a)
        else:lo=max(lo,c/a)
    if lo>hi:return None
    return float(lo),float(hi)

def free_reporter_interval(total,capacity):
    if not 0<=capacity<=1:raise ValueError('capacity must be in [0,1]')
    if not 0<=total[0]<=total[1]<=1:raise ValueError('fractions must be ordered in [0,1]')
    # qtot=(1-b)qfree+b*qbound, unknown qbound in [0,1], b<=capacity.
    # This sharp interval is narrower than simply qtot +/- capacity.
    if capacity==1:return (0.,1.)
    return (max(0.,(total[0]-capacity)/(1-capacity)),min(1.,total[1]/(1-capacity)))

def mean(log_u,coeff):
    idx=np.arange(len(coeff));weights=np.log(coeff)+idx*log_u
    p=np.exp(weights-logsumexp(weights));return float(p@idx/(len(coeff)-1))

def calibrated_free(qb,cb,ck):
    if qb==0:return 0.
    if qb==1:return 1.
    x=brentq(lambda x:mean(x,cb)-qb,-200,200,xtol=1e-12)
    return mean(x,ck)

def verify():
    rng=np.random.default_rng(924)
    worst=0.;cases=3000
    for _ in range(cases):
        free=np.sort(rng.uniform(0,1,2));bound=np.sort(rng.uniform(0,1,2));total=np.sort(rng.uniform(0,1,2))
        got=mixture_interval(free,bound,total)
        # Independent LP in s, free-pool modified mass x, bound modified mass y.
        A=[[free[0],-1,0],[-free[1],1,0],[-bound[0],0,-1],[-bound[1],0,1]]
        rhs=[-free[0],free[1],0,0]
        # First pair derives x >= (1-s)flo and x <= (1-s)fhi.
        A[0]=[-free[0],-1,0];A[1]=[free[1],1,0]
        # y >= s*blo.
        A[2]=[bound[0],0,-1]
        A += [[0,1,1],[0,-1,-1]];rhs += [total[1],-total[0]]
        low=linprog([1,0,0],A_ub=A,b_ub=rhs,bounds=[(0,1),(0,None),(0,None)],method='highs')
        high=linprog([-1,0,0],A_ub=A,b_ub=rhs,bounds=[(0,1),(0,None),(0,None)],method='highs')
        assert low.success==high.success
        if not low.success:assert got is None;continue
        assert got is not None
        worst=max(worst,abs(got[0]-low.fun),abs(got[1]+high.fun))
    assert worst<1e-8
    err=0.
    for _ in range(1000):
        cb=np.exp(rng.uniform(-6,6,rng.integers(2,8)));ck=np.exp(rng.uniform(-6,6,rng.integers(2,14)))
        x=rng.uniform(-5,5);qb=mean(x,cb);fk=mean(x,ck);s=rng.uniform();bk=rng.uniform(0,.2)
        pred=calibrated_free(qb,cb,ck);err=max(err,abs(pred-fk))
        k=(1-s)*fk+s*bk
        interval=mixture_interval((max(0,pred-1e-9),min(1,pred+1e-9)),(bk,bk),(k,k))
        assert interval is not None and interval[0]-1e-7<=s<=interval[1]+1e-7
    assert err<1e-8
    reporter_error=0.
    for _ in range(500):
        qlo,qhi=np.sort(rng.uniform(0,1,2));cap=rng.uniform(0,.999)
        got=free_reporter_interval((qlo,qhi),cap)
        # Feasible free fractions satisfy (1-cap)f <= qhi and
        # (1-cap)f+cap >= qlo; independently optimize the one-variable LP.
        A=[[1-cap],[-(1-cap)]];rhs=[qhi,cap-qlo]
        low=linprog([1],A_ub=A,b_ub=rhs,bounds=[(0,1)],method='highs')
        high=linprog([-1],A_ub=A,b_ub=rhs,bounds=[(0,1)],method='highs')
        assert low.success and high.success
        reporter_error=max(reporter_error,abs(got[0]-low.fun),abs(got[1]+high.fun))
    assert reporter_error<1e-10
    assert free_reporter_interval((.2,.8),0)==(.2,.8)
    assert free_reporter_interval((.2,.8),1)==(0.,1.)
    assert mixture_interval((0,0),(0,0),(0,0))==(0.,1.)
    assert mixture_interval((1,1),(1,1),(1,1))==(0.,1.)
    assert mixture_interval((0,0),(1,1),(.2,.8))==(.2,.8)
    assert mixture_interval((.75,.75),(0,0),(.375,.375))==(.5,.5)
    assert mixture_interval((.5,.5),(.5,.5),(.5,.5))==(0.,1.)
    assert mixture_interval((.5,.5),(.5,.5),(.6,.6)) is None
    cases_out=[]
    for bmax in [0,.05,.15,.5,1]:
        cases_out.append({'bound_modification_max':bmax,'s_interval':mixture_interval((.70,.80),(0,bmax),(.35,.39))})
    result={'reporter_lp_cases':500,'max_reporter_lp_discrepancy':reporter_error,'mixture_lp_cases':cases,'max_lp_discrepancy':worst,'arbitrary_ladder_pairs':1000,'max_curve_inversion_error':err,'synthetic_selectivity_scenarios':cases_out,'scope':'Observation-level identifiability under externally justified calibration/interval coverage. No native occupancy validation.'}
    out=Path(__file__).resolve().parents[1]/'results'/'robust_calibration_verification.json';out.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':verify()
