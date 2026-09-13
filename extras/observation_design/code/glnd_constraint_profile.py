from pathlib import Path
import numpy as np,json,csv
from scipy.optimize import least_squares
from scipy.special import expit
R=Path(__file__).resolve().parents[1];OUT=R/'results';x,y=np.loadtxt(R/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1).T
r0=(137/3)/(2.7/2.3);b2=(6.5/.82)/(2.7/2.3);rng=np.random.default_rng(355);rows=[]
for k1 in [.026,.1,.3,1.,3.,10.,100.]:
 def decode(z):
  k2=np.exp(z[0]);a1=r0*expit(z[1]);a2=a1*expit(z[2]);b1=1+(b2-1)*expit(z[3]);return k2,np.array([r0,a1,a2]),np.array([1,b1,b2])
 def pred(z):
  k2,a,b=decode(z);w=np.array([np.ones_like(x),x/k1,x*x/(k1*k2)]).T;return 3*(w@a)/(w@(a+b))
 fits=[least_squares(lambda z:pred(z)-y,np.r_[rng.uniform(-8,1),rng.uniform(-4,3,3)],bounds=([-20]*4,[10,20,20,20]),max_nfev=2000,ftol=1e-11,xtol=1e-11,gtol=1e-11)for _ in range(12)]
 f=min(fits,key=lambda f:sum(f.fun**2));k2,a,b=decode(f.x);G=np.sqrt(k1*k2);p1=(G/k1)/(2+G/k1);lo=max(a[2]/a[1],1/b[1]);hi=min(r0/a[1],b2/b[1]);p1lo=(G/(k1*lo))/(2+G/(k1*lo));p1hi=(G/(k1*hi))/(2+G/(k1*hi))
 rows.append(dict(K1_mM=k1,K2_mM=float(k2),RMSE=float(np.sqrt(np.mean(f.fun**2))),midpoint_mM=float(G),p1_at_midpoint=float(p1),lambda_low=float(lo),lambda_high=float(hi),max_p1_separation=float(p1lo-p1hi)))
with (OUT/'glnd_constraint_profile.csv').open('w')as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))
