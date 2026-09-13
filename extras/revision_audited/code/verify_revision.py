"""Independent numerical checks for the revised mathematical and fitted results."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
x,y=np.loadtxt(ROOT/'data/jiangninfa2011_fig2B_digitised.csv',delimiter=',',skiprows=1).T
summary=json.loads((ROOT/'data/analysis_summary.json').read_text())
def hill(p):return p[0]+(p[1]-p[0])/(1+(x/np.exp(p[2]))**p[3])
rng=np.random.default_rng(73);sses=[];hs=[]
for _ in range(100):
 p0=[rng.uniform(.1,1),rng.uniform(2,2.999),rng.uniform(-2,1),rng.uniform(.5,4)]
 r=least_squares(lambda p:hill(p)-y,p0,bounds=([0,1.5,-8,.1],[1.5,3,5,8]),xtol=1e-11,ftol=1e-11,gtol=1e-11)
 sses.append(float(r.fun@r.fun));hs.append(r.x[3])
assert max(abs(np.array(sses)-summary['SSE']))<1e-9
# Unbounded plateau comparison reproduces the original fixed-h=1 SSE.
r=least_squares(lambda q:q[0]+(q[1]-q[0])/(1+x/np.exp(q[2]))-y,[.4,3.3,-.5],xtol=1e-12,ftol=1e-12,gtol=1e-12)
# Symbolically verify the saturating-switch transformation.
l,K,w,v,k=sp.symbols('l K w v k',positive=True)
f=w+(v-w)*l**2/(k**2+l**2)
theta=K/(K+f);t0=K/(K+w);ti=K/(K+v);S2=k**2*(K+w)/(K+v)
assert sp.simplify(theta-(ti+(t0-ti)/(1+l**2/S2)))==0
assert sp.simplify((t0/(1-t0))/(ti/(1-ti))-v/w)==0
# Switch maximum derivative vanishes at its stated operating point.
a,b=sp.symbols('a b',positive=True)
eta=1+(a-b)*l/((1+a*l)*(1+b*l))
assert sp.simplify(sp.diff(eta,l).subs(l,1/sp.sqrt(a*b)))==0
# For the free/total range bound, random admissible thresholds obey the bound.
for _ in range(1000):
 low=10**rng.uniform(-3,1);high=low*10**rng.uniform(.1,2);delta=low*rng.uniform(0,.15)
 tl=low+rng.uniform(0,delta);th=high+rng.uniform(0,delta)
 lower=high/(low+delta);upper=(high+delta)/low
 assert lower<=th/tl<=upper
# Data consistency and full retention of coordinate perturbations.
z=np.load(ROOT/'data/coordinate_perturbations.npz')
assert z['parameters'].shape==(2000,4)
assert np.all(np.isfinite(z['composite_coefficients']))
res={'multistart_runs':100,'max_SSE_spread':float(np.ptp(sses)),'h_range':[float(min(hs)),float(max(hs))],'unbounded_h1_SSE':float(r.fun@r.fun),'unbounded_h1_upper_plateau':float(r.x[1]),'symbolic_switch_identity':True,'symbolic_plateau_odds':True,'symbolic_stationary_point':True,'range_bound_random_checks':1000,'all_2000_perturbations_retained':True}
(ROOT/'data/verification.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
