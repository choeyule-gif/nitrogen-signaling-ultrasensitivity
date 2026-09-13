"""A positive count-state GlnB/GlnK model with explicit K0 sequestration.
Nine fractions: B0..B3, K0..K3, AmtB-bound K0. Each paralog's fractions sum to 1.
Buffered metabolite drivers and measured protein totals are external inputs.
Only fully unmodified K0 can bind. Binding is pseudo-first-order: free AmtB is
buffered, NOT an estimated conserved transporter pool. Protein abundance rise
adds unmodified trimers; decline removes all states proportionately. Neither
choice identifies underlying synthesis/degradation mechanisms.
"""
from pathlib import Path
import json,csv
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.special import comb
from dynamic_pilot import inputs,params,obs
R=Path(__file__).resolve().parents[1];idx=np.arange(4);un=(3-idx)/3;mod=idx/3

def data(strain):
 rows=[]
 for ri in sorted({int(x['source_row']) for x in obs if x['strain']==strain}):
  d={x['observable']:x for x in obs if int(x['source_row'])==ri}
  for j,p in enumerate(['GlnB','GlnK']):
   if (strain=='delta_glnB' and j==0) or (strain=='delta_glnK' and j==1):continue
   try:
    a,b=d[p+'_UMP'],d[p+'_unmodified'];u,v=float(a['mean_copy_cell']),float(b['mean_copy_cell']);su,sv=float(a['SE_copy_cell']),float(b['SE_copy_cell'])
   except (KeyError,ValueError):continue
   if u+v<=0:continue
   se=(v*su+u*sv)/(u+v)**2
   if se<=0:continue
   rows.append([float(a['time_min']),j,u/(u+v),se])
 return np.array(rows)

def interpolate(t,xy):
 x,y=xy;j=min(max(np.searchsorted(x,t,side='right')-1,0),len(x)-2)
 return float(np.interp(t,x,y)),float((y[j+1]-y[j])/(x[j+1]-x[j])) if x[0]<=t<=x[-1] else 0.

def simulate(strain,lam,logs,method='LSODA',rtol=2e-6):
 ut_gain,ur_gain,beta,koff=np.exp(logs[:4]);specificity=np.exp(logs[4]) if len(logs)>4 else 1.;d=inputs(strain);dt=data(strain)
 runout={'WT':-40.,'delta_glnB':-47.333333333333336,'delta_glnK':-52.06666666666667}[strain]
 end=max(dt[:,0]);dt=dt[dt[:,0]>=runout-1e-7]
 y=np.zeros(9)
 for j in [0,1]:
  z=dt[dt[:,1]==j];q=float(z[0,2]) if len(z) else 0.;y[j*4:j*4+4]=comb(3,idx)*q**idx*(1-q)**(3-idx)
 times=sorted(set([runout,end]+[float(t) for x,v in d.values() for t in x if runout<t<end]));sols=[]
 def rhs(t,z):
  bt,db=interpolate(t,d['p2t']);kt,dk=interpolate(t,d['gt']);g=np.interp(t,*d['q'])/params['kgln1'];ak=np.interp(t,*d['akg'])
  bu,ku=bt*(un@z[:4]),kt*(un@z[4:8]);bm,km=bt*(mod@z[:4]),kt*(mod@z[4:8])
  f=(1+2*g+g*g)/(1+2*g/lam+g*g)
  a=f*params['vmaxut']*ut_gain/((1+bu/params['kmp2']+ku/params['kmg'])*(1+g))
  b=f*params['vmaxur']*ur_gain/((1+bm/params['kmp2u']+km/params['kmgu'])*(1+1/g))
  dz=np.zeros(9)
  for offset,kp,km_,T,Tdot in [(0,params['kmp2'],params['kmp2u'],bt,db),(4,params['kmg'],params['kmgu'],kt,dk)]:
   for i in range(3):
    flux=(3-i)*a/kp*z[offset+i]-(i+1)*b/km_*(specificity if offset==4 else 1.)*z[offset+i+1];dz[offset+i]-=flux;dz[offset+i+1]+=flux
   if T>0 and Tdot>0:
    # Positive production of unmodified protein, normalized by the imposed total.
    n=4 if offset==0 else 5;dz[offset:offset+n]-=Tdot/T*z[offset:offset+n];dz[offset]+=Tdot/T
  bind=koff*(beta/(1+(ak/params['kakg1'])**params['nut']))*z[4]-koff*z[8]
  dz[4]-=bind;dz[8]+=bind
  return dz
 for lo,hi in zip(times[:-1],times[1:]):
  sol=solve_ivp(rhs,[lo,hi],y,method=method,rtol=rtol,atol=rtol*1e-3,dense_output=True)
  if not sol.success:raise RuntimeError(sol.message)
  y=sol.y[:,-1];sols.append(sol)
 predictions=[]
 for t,j,q,se in dt:
  k=min(max(np.searchsorted(times,t,side='right')-1,0),len(sols)-1);z=sols[k].sol(t);offset=int(j)*4
  predictions.append([t,j,q,se,mod@z[offset:offset+4],z[8]])
 allz=np.concatenate([z.sol(np.linspace(z.t[0],z.t[-1],31)).T for z in sols]);mass=max(np.max(abs(allz[:,:4].sum(1)-1)),np.max(abs(allz[:,4:].sum(1)-1)))
 assert allz.min()>-1e-6 and mass<1e-6
 return np.array(predictions),dict(min_state=float(allz.min()),mass_error=float(mass),max_bound=float(allz[:,8].max()))

def main():
 report=[];out=[]
 for lam in [1.,2.,4.]:
  def residual(v):
   p,_=simulate('WT',lam,v);return (p[2:,4]-p[2:,2])/p[2:,3] # initial two observations are conditioning values
  starts=[np.log([1,1,1.5,1]),np.log([.5,.5,10,.1])]
  fits=[least_squares(residual,v,bounds=([-4,-4,-4,-5],[4,4,6,5]),diff_step=.005,ftol=3e-5,xtol=3e-5,gtol=3e-5,max_nfev=70) for v in starts]
  fit=min(fits,key=lambda f:f.fun@f.fun)
  for strain in ['WT','delta_glnK','delta_glnB']:
   z,qa=simulate(strain,lam,fit.x);ninit=2 if strain=='WT' else 1;score=z[ninit:]
   rec=dict(lam=lam,strain=strain,parameters=np.exp(fit.x).tolist(),n_scored=len(score),weighted_SSE=float(np.sum(((score[:,4]-score[:,2])/score[:,3])**2)),RMSE_fraction=float(np.sqrt(np.mean((score[:,4]-score[:,2])**2))),optimizer_success=bool(fit.success),**qa);report.append(rec);print(rec,flush=True)
   out.extend(dict(lam=lam,strain=strain,time_min=t,protein='GlnB' if j==0 else 'GlnK',observed_fraction=q,SE_upper=se,predicted_fraction=p,bound_fraction=s) for t,j,q,se,p,s in z)
 (R/'results/native_trimer_summary.json').write_text(json.dumps(report,indent=2))
 with (R/'results/native_trimer_predictions.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
if __name__=='__main__':main()
