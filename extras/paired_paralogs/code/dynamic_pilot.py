"""Independent in-vivo data: pilot test of stationary-equivalent UT/UR timing.
Original rate model: Figshare 4880003 model.txt (CC BY 4.0).
The lift w=(1,2g,g²), UT=(1,1/2,0), UR=(0,1/2,1) realizes original
1/(1+g), g/(1+g) regulatory factors exactly. Transform w1 -> w1/lam,
a1,b1 -> lam*a1,lam*b1 multiplies both rates by the same F_lam(g).
Only a common speed is calibrated to WT. Mutant data are held out.
Input means treated as fixed; summary SE is a working residual scale.
"""
import re,json,csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
R=Path(__file__).resolve().parents[1];src=(R/'sources/model.txt').read_text()
with (R/'results/observations.csv').open() as f:obs=list(csv.DictReader(f))
params={k:float(v) for k,v in re.findall(r'^([a-zA-Z0-9_]+)\s*=\s*([0-9.]+)\s*$',src,re.M)}
drives={}
for name,x,y in re.findall(r'^(p2t|gt|q|akg|gst)=interp1SB\(\[([^]]+)\],\[([^]]+)\],time\)',src,re.M):
 drives[name]=(np.array([float(z) for z in x.split(',')])-101,np.array([float(z) for z in y.split(',')]))
def inputs(strain):
 if strain=='WT':return drives
 out={}
 for name,obsname,mult in [('p2t','GlnB_total',.00161),('gt','GlnK_total',.00161),('q','glutamine',.00000161),('akg','aKG',.00000161)]:
  z=sorted([(float(x['time_min']),float(x['mean_copy_cell'])*mult) for x in obs if x['strain']==strain and x['observable']==obsname]);out[name]=(np.array([x for x,y in z]),np.array([y for x,y in z]))
 return out
def observed(strain):
 z=[]
 for x in obs:
  if x['strain']!=strain or x['observable'] not in ['GlnB_UMP','GlnK_UMP']:continue
  if strain=='delta_glnK' and x['observable']=='GlnK_UMP':continue
  if strain=='delta_glnB' and x['observable']=='GlnB_UMP':continue
  try:se=float(x['SE_copy_cell'])*.00161
  except ValueError:continue
  if se<=0:continue
  z.append((float(x['time_min']),0 if x['observable']=='GlnB_UMP' else 1,float(x['mean_copy_cell'])*.00161,se))
 return np.array(z)
def simulate(strain,lam,speed):
 d=inputs(strain);ts=observed(strain);t0=max(v[0][0] for v in d.values()) if strain!='WT' else -101
 # Start at earliest jointly driven observation; initialise each measured modified pool by interpolation.
 if strain=='WT':y0=[1.4,1.31]
 else:
  y0=[]
  for j in [0,1]:
   z=ts[ts[:,1]==j];y0.append(float(np.interp(t0,z[:,0],z[:,2])) if len(z) else 0.)
 end=max(x[-1] for x,y in d.values());times=sorted(set(np.r_[t0,end,*[x for x,y in d.values()]]));times=[x for x in times if t0<=x<=end]
 minrate=0.;minpool=0.
 def rhs(t,z):
  nonlocal minrate,minpool
  b,k,q,ak=[np.interp(t,*d[x]) for x in ['p2t','gt','q','akg']]
  gts=k*(params['x']+(1-params['x'])*ak**params['nut']/(params['kakg1']**params['nut']+ak**params['nut']))
  b0,k0=b-z[0],gts-z[1];g=q/params['kgln1'];fac=(1+2*g+g*g)/(1+2*g/lam+g*g)*speed
  ut=fac*params['vmaxut']*np.array([b0/params['kmp2'],k0/params['kmg']])/((1+b0/params['kmp2']+k0/params['kmg'])*(1+g))
  ur=fac*params['vmaxur']*np.array([z[0]/params['kmp2u'],z[1]/params['kmgu']])/((1+z[0]/params['kmp2u']+z[1]/params['kmgu'])*(1+1/g))
  minrate=min(minrate,float(ut.min()),float(ur.min()));minpool=min(minpool,b0,k0,float(np.min(z)))
  return ut-ur
 sols=[]
 for lo,hi in zip(times[:-1],times[1:]):
  sol=solve_ivp(rhs,[lo,hi],y0,method='LSODA',rtol=2e-7,atol=2e-9,dense_output=True);assert sol.success;y0=sol.y[:,-1];sols.append(sol)
 accepted_min=0.
 for sol in sols:
  for tt,zz in zip(np.linspace(sol.t[0],sol.t[-1],31),sol.sol(np.linspace(sol.t[0],sol.t[-1],31)).T):
   b,k,ak=[np.interp(tt,*d[x]) for x in ["p2t","gt","akg"]]
   available=k*(params["x"]+(1-params["x"])*ak**params["nut"]/(params["kakg1"]**params["nut"]+ak**params["nut"]))
   accepted_min=min(accepted_min,b-zz[0],available-zz[1],float(zz.min()))
 pred=[]
 for t,j,y,se in ts:
  if t<t0 or t>end:continue
  ix=min(max(np.searchsorted(times,t,side='right')-1,0),len(sols)-1);pred.append([t,j,y,se,float(sols[ix].sol(t)[int(j)])])
 return np.array(pred),dict(min_accepted_available_pool=accepted_min,min_evaluated_rate=minrate,min_evaluated_pool=minpool,initial_time=t0,interpretation='Negative trial RHS values may include solver stages; inspect accepted trajectories before biochemical interpretation.')
if __name__ == '__main__':
 rows=[];report=[]
 for lam in [1.,2.,4.]:
  def loss(logs):
   z,_=simulate('WT',lam,np.exp(logs));return np.sum(((z[:,4]-z[:,2])/z[:,3])**2)
  fit=minimize_scalar(loss,bounds=(-3,3),method='bounded',options={'xatol':.002});speed=np.exp(fit.x)
  for strain in ['WT','delta_glnK','delta_glnB']:
   z,qa=simulate(strain,lam,speed);record=dict(lam=lam,strain=strain,speed=speed,n=len(z),weighted_SSE=float(np.sum(((z[:,4]-z[:,2])/z[:,3])**2)),RMSE_uM=float(np.sqrt(np.mean((z[:,4]-z[:,2])**2))),**qa);report.append(record);print(record,flush=True)
   rows.extend(dict(lam=lam,strain=strain,speed=speed,time_min=t,observable='GlnB_UMP' if j==0 else 'GlnK_UMP',observed_uM=y,SE_uM=se,predicted_uM=p) for t,j,y,se,p in z)
 with (R/'results/dynamic_pilot_predictions.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 (R/'results/dynamic_pilot_summary.json').write_text(json.dumps(report,indent=2))
