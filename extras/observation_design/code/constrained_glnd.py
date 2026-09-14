"""Cross-assay stress fit retaining two endpoint specificity ratios.
Fixed effective C, the earlier middle-state constraint and the transported
0.080 mM UR half-range are released together.
Transported ratios are sensitivity hypotheses, not matched-assay constraints.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
x,y=np.loadtxt(ROOT/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1).T
r0=(137/3)/(2.7/2.3);b2=(6.5/.82)/(2.7/2.3)

def decode(z):
 k1,k2=np.exp(z[:2]);a1=r0*expit(z[2]);a2=a1*expit(z[3]);b1=1+(b2-1)*expit(z[4]);return k1,k2,np.array([r0,a1,a2]),np.array([1,b1,b2])
def response(g,z):
 k1,k2,a,b=decode(z);w=np.stack([np.ones_like(g),g/k1,g*g/(k1*k2)],axis=-1);return 3*(w@a)/(w@(a+b))
rng=np.random.default_rng(59);best=None
for _ in range(50):
 z=np.r_[rng.uniform(-4,2,2),rng.uniform(-5,3,3)];f=least_squares(lambda z:response(x,z)-y,z,bounds=([-12,-12,-20,-20,-20],[5,5,20,20,20]),max_nfev=3000,ftol=1e-12,gtol=1e-12,xtol=1e-12)
 if best is None or np.sum(f.fun**2)<np.sum(best.fun**2):best=f
k1,k2,a,b=decode(best.x)
lo=max(a[2]/a[1],b[0]/b[1]);hi=min(a[0]/a[1],b[2]/b[1])
g=np.geomspace(1e-5,1e3,500);base=response(g,best.x);worst=0;family=[]
for lam in np.linspace(lo,hi,5):
 w=np.stack([np.ones_like(g),g/(k1*lam),g*g/(k1*k2)],axis=-1);aa=a.copy();bb=b.copy();aa[1]*=lam;bb[1]*=lam
 pred=3*(w@aa)/(w@(aa+bb));worst=max(worst,float(np.max(abs(pred-base))))
 family.append(dict(lam=float(lam),a=aa.tolist(),b=bb.tolist(),K1_mM=float(k1*lam),K2_mM=float(k2/lam)))
assert worst<1e-12
report=dict(ratio_unliganded=r0,UR_fold=b2,RMSE=float(np.sqrt(np.mean(best.fun**2))),K1_mM=k1,K2_mM=k2,a=a.tolist(),b=b.tolist(),monotone_capacity_lambda_interval=[lo,hi],family=family,max_equivalence_error=worst,zero_input_max_UMP=3*r0/(1+r0),predictions=response(x,best.x).tolist(),scope='Transfer of two 1998 endpoint apparent-specificity ratios with monotone UT/UR sequences; fixed effective C, the earlier intermediate-state constraint and the 0.080 mM UR half-range are not imposed. Cross-assay transport is not validated; fit to 12 rounded coordinates is descriptive.')
(OUT/'constrained_glnd.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
