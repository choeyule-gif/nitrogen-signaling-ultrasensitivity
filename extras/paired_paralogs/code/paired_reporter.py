"""Paired paralog inference, conditional on independent free-site ladders.
No claim that the 2017 time points are steady. Data application is a sensitivity
map over unknown rho; it is not an estimate of native sequestration.
"""
from pathlib import Path
import json,csv
import numpy as np
from scipy.special import comb
R=Path(__file__).resolve().parents[1]
def free_k(qb,rho):return rho*qb/(1-qb+rho*qb)
def infer(qb,qk,rho):return 1-qk/free_k(qb,rho)
def stationary(qb,rho,beta):
 u=qb/(1-qb);v=rho*u
 b=comb(3,np.arange(4))*u**np.arange(4);b/=b.sum()
 weights=np.r_[comb(3,np.arange(4))*v**np.arange(4),beta]
 k=weights/weights.sum();return b,k

def independent_solve(ut,ur,beta,off=1):
 q=np.zeros((5,5))
 for i in range(3):q[i+1,i]=(3-i)*ut;q[i,i+1]=(i+1)*ur
 q[4,0]=beta*off;q[0,4]=off;q[np.diag_indices(5)]=-q.sum(0)
 a=q.copy();a[-1]=1;y=np.zeros(5);y[-1]=1;return np.linalg.solve(a,y),q
rng=np.random.default_rng(9214);err=0.;mass=0.;n=2000
for _ in range(n):
 qb=rng.uniform(.02,.98);rho=10**rng.uniform(-1,1);beta=10**rng.uniform(-2,3)
 b,k=stationary(qb,rho,beta);qk=k[:4]@np.arange(4)/3;s=k[4]
 z,Q=independent_solve(rho*qb/(1-qb),1,beta,10**rng.uniform(-2,2))
 err=max(err,abs(infer(qb,qk,rho)-s),np.max(abs(z-k)));mass=max(mass,np.max(abs(Q@k)))
 # An independent no-AmtB condition calibrates rho even when its input weight differs.
 qb0=rng.uniform(.03,.97);qk0=free_k(qb0,rho);cal=qk0/(1-qk0)/(qb0/(1-qb0));assert abs(cal/rho-1)<1e-12
assert err<1e-9 and mass<1e-8
with (R/'results/observations.csv').open() as f:raw=list(csv.DictReader(f))
rows=[]
for row in sorted({z['source_row'] for z in raw if z['strain']=='WT'},key=int):
 d={z['observable']:z for z in raw if z['source_row']==row}
 try:
  def frac(p):
   u,v=[float(d[p+x]['mean_copy_cell']) for x in ['_UMP','_unmodified']];return u/(u+v)
  qb,qk=frac('GlnB'),frac('GlnK')
 except (KeyError,ZeroDivisionError):continue
 for rho in [1,1.5,3,9,27]:
  ss=infer(qb,qk,rho);rows.append(dict(time_min=d['GlnB_UMP']['time_min'],rho_scenario=rho,theta_B=qb,theta_K=qk,inferred_bound_fraction=ss,compatible_point=0<=ss<=1,source_row=row))
with (R/'results/paired_reporter_scenarios.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(R/'results/paired_reporter_verification.json').write_text(json.dumps(dict(random_cases=n,max_reconstruction_error=err,max_stationary_residual=mass,calibration='No-AmtB control identifies rho conditional on unchanged relative catalytic specificities and independent-site free ladders.',data_interpretation='Scenario map only. Neither stationarity nor rho is established by the 2017 data; negative inferred fractions reject the chosen combination of assumptions, not the biological data.'),indent=2));print('Paired reporter verified',n,err,mass)
