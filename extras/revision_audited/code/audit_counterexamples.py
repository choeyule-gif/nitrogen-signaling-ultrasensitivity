from pathlib import Path
import json,sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'data/audit'; mp.mp.dps=80
src=json.loads((OUT/'counterexamples_parameters.json').read_text())
def prepare(d,kind):
 n=d['n']; A=mp.mpf(str(d['A'])); Bh=list(map(lambda x:mp.mpf(str(x)),d['Bhat'][:n]));
 if kind=='i':
  Gh=list(map(lambda x:mp.mpf(str(x)),d['Ghat'])); rho=list(map(lambda x:mp.mpf(str(x)),d['rho'])); K=[rho[i]*Bh[i]/(A*Gh[i]) for i in range(n)]
 else:
  rho=[mp.mpf(str(d['rho']))]*n;K=list(map(lambda x:mp.mpf(str(x)),d['K']));Gh=[rho[i]*Bh[i]/(A*K[i]) for i in range(n)]
 c=[mp.mpf(1)]
 for k in K:c.append(c[-1]*k)
 TE,TS,TL=[mp.mpf(str(d[k])) for k in ['TE','TS','TL']]
 def state(l):
  rr=[c[i]/l**i for i in range(n+1)]; M=sum(rr); W=sum(Bh[i]*rr[i] for i in range(n)); V=sum(Gh[i]*rr[i+1]*l for i in range(n)); Q=W+A*V
  aa=(1+A*l)*Q;bb=(1+A*l)*M+(TS-TE)*Q
  # Rationalized positive root avoids catastrophic cancellation.
  e=2*TE*M/(bb+mp.sqrt(bb**2+4*aa*TE*M));s0=TS/(M+e*Q);s=[s0*x for x in rr];eps=A*e*l
  beta=[Bh[i]*e*s[i] for i in range(n)];gamma=[Gh[i]*eps*s[i+1] for i in range(n)]
  return l+A*e*l+sum(gamma)-TL,e,eps,s,beta,gamma
 roots=[]
 for old in d['roots_double']:
  rt=mp.findroot(lambda l:state(l)[0],(mp.mpf(str(old))*.99,mp.mpf(str(old))*1.01),tol=mp.mpf('1e-70'))
  v,e,eps,s,b,g=state(rt)
  flux=[rho[i]*b[i]-g[i] for i in range(n)]
  residuals=[A*e*rt-eps,e+eps+sum(b)+sum(g)-TE,rt+eps+sum(g)-TL,sum(s)+sum(b)+sum(g)-TS,*flux]
  # Raw mass action: binding, dissociation, and catalysis with positive rates.
  residuals += [Bh[i]*(1+rho[i])*e*s[i]-(1+rho[i])*b[i] for i in range(n)]
  residuals += [2*Gh[i]*eps*s[i+1]-2*g[i] for i in range(n)]
  p=[z/sum(s) for z in s];mu=sum(i*p[i] for i in range(n+1));Wmean=sum(Bh[i]*p[i] for i in range(n));cov=sum(Bh[i]*i*p[i] for i in range(n))-Wmean*mu
  roots.append(dict(l=mp.nstr(rt,35),e=mp.nstr(e,20),min_species=mp.nstr(min([rt,e,eps,*s,*b,*g]),5),max_balance_residual=mp.nstr(max(map(abs,residuals)),5),omega_elasticity=float(-cov/Wmean)))
 # Global extrema of elasticity on a broad domain, endpoints tend to 1 and 0.
 cf=np.array(list(map(float,c)));bh=np.r_[list(map(float,Bh)),0];ii=np.arange(n+1)
 def elast(z):
  w=np.log(cf)-ii*z; pp=np.exp(w-logsumexp(w)); mm=pp@ii; ww=pp@bh;return -((pp*ii)@bh-ww*mm)/ww
 zz=np.linspace(-30,30,3001);vv=np.array([elast(z) for z in zz]);mn=vv.min();mx=vv.max()
 for i in np.where((vv[1:-1]>vv[:-2])&(vv[1:-1]>vv[2:]))[0]+1:mx=max(mx,-minimize_scalar(lambda z:-elast(z),bounds=(zz[i-1],zz[i+1]),method='bounded').fun)
 for i in np.where((vv[1:-1]<vv[:-2])&(vv[1:-1]<vv[2:]))[0]+1:mn=min(mn,minimize_scalar(elast,bounds=(zz[i-1],zz[i+1]),method='bounded').fun)
 return dict(n=n,rho_spread=float(max(rho)/min(rho)),omega_elasticity_sample_and_refined=[mn,mx],roots=roots)
out={'condition_i':prepare(src['condition_i'],'i'),'condition_ii':[prepare(x,'ii') for x in src['condition_ii']]}
# Audit the original S5A caption against the actual constants.
k=np.array([3.,1.,1/3]);rho=k*.8/1.5
out['original_S5A']={'rho':rho.tolist(),'rho_spread':float(rho.max()/rho.min()),'condition_i_satisfied':False,'correction':'Bhat_i=(3-i)*1.5, Ghat_i=(i+1)*0.8 gives constant rho=0.8/1.5 for Khat=A=1.'}
(OUT/'omission_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

replacement=json.loads((OUT/'replacement_condition_i.json').read_text())
(OUT/'replacement_condition_i_checks.json').write_text(json.dumps(prepare(replacement,'i'),indent=2))
