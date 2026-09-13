"""Fresh calculations for the v1.1 research extensions (no stored pass flags)."""
from pathlib import Path
import importlib.util
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix
R=Path(__file__).resolve().parents[1]
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,R/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def test_regulatory_equivalences():
 rc=module('regulatory_core','extras/robust_design/code/regulatory_core.py')
 x,y=np.loadtxt(R/'extras/regulatory_states/data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1).T
 v,sse=rc.fit(x,y);assert np.sqrt(sse/len(x))<.037
 grid=np.geomspace(.002,20,71);a,b,o=rc.glnd(grid,v)
 for lam in [2,4]:
  c,d,p=rc.glnd(grid,v,lam);np.testing.assert_allclose(a/(a+b),c/(c+d),atol=2e-14);assert np.max(np.abs(o-p))>.1
 par=dict(KG=15.6,KP=.18,KU=.35,a1=.17,a2=4,beta=1,betap=10,k=1)
 rng=np.random.default_rng(6102)
 for _ in range(40):
  g,p,u=10**rng.uniform(-2,1,3);at,ar,occ=rc.glne(g,p,u,par);z,err=rc.receptor_stationary(g,p,u,par);np.testing.assert_allclose(z,occ,atol=1e-11);assert err<1e-10
  alt={**par,'a1':par['a1']*4,'betap':par['betap']*4};bt,br,_=rc.glne(g,p,u,alt);np.testing.assert_allclose(at/ar,bt/br,rtol=1e-12)

def test_closed_cascade_and_productive_boundary():
 m=module('closed_model','extras/closed_cascade/code/model.py')
 for productive in [0,1]:
  for g in [20,530,3000]:
   for j in [0,.4]:
    kw=dict(G=g,J=j,Dt=.1,Et=.05,productive=productive)
    a,b=m.Cascade(**kw),m.Cascade(**kw,lam=4,sigma=4)
    ya,yb=a.steady(),b.steady();oa,ob=a.obs(ya),b.obs(yb)
    for model,y in [(a,ya),(b,yb)]:
     assert np.min(y)>=0;assert np.max(np.abs(model.rhs(0,y)))<1e-8;np.testing.assert_allclose(model.cons@y,[5,10,.1,.05],atol=1e-9);np.testing.assert_allclose(model.cons@model.N,0,atol=1e-14)
    np.testing.assert_allclose(oa['pfree'],ob['pfree'],atol=1e-12)
    eps=.04;nu=1 if j==0 else 12;bound=np.tanh(nu/4*(-np.log1p(-eps)))+.005
    assert abs(oa['GStotal']-ob['GStotal'])<=bound+1e-12
 # Integrate full reactions from an independently assigned initial state.
 a=m.Cascade(G=530,Dt=.1,Et=.05);z=np.zeros(len(a.names));z[a.P[0]]=5;z[a.S[0]]=10;z[a.D[0]]=.1;z[a.E[0]]=.05
 dep=np.zeros((len(a.ks),len(z)))
 for k,ix in enumerate(a.rs):dep[k,ix]=1
 sol=solve_ivp(a.rhs,[0,2e6],z,method='BDF',rtol=2e-8,atol=2e-11,jac_sparsity=csr_matrix(np.abs(a.N)@dep>0))
 assert sol.success;np.testing.assert_allclose(sol.y[:,-1],a.steady(),atol=3e-7)

def test_ladder_and_design_bounds():
 from scipy.stats import norm
 rng=np.random.default_rng(931)
 for _ in range(500):
  n=int(rng.integers(1,14));i=np.arange(n+1);c=rng.normal(0,8,n+1);x=rng.uniform(-8,8);delta=rng.uniform(0,3)
  def mean(t):
   z=c+i*t;v=np.exp(z-z.max());return v@i/(v.sum()*n)
  assert abs(mean(x+delta)-mean(x))<=np.tanh(n*delta/4)+1e-12
 zcrit=norm.ppf(.95);power=norm.cdf(np.sqrt(155)*.01/.05-zcrit);assert power>=.8
 # Occupancy separation with independently varying K nuisance ranges.
 g=np.geomspace(.5,2,1001)
 null=2*g/(1+2*g+g*g);alt=.5*g/(1+.5*g+g*g)
 assert abs(null.min()-alt.max()-11/45)<1e-8

if __name__=='__main__':
 for name,fn in sorted(list(globals().items())):
  if name.startswith('test_'):
   fn();print('PASS',name,flush=True)
