from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
import numpy as np,json
from scipy.special import logsumexp
from scipy.optimize import least_squares,brentq
from pathlib import Path
data=json.loads((ROOT/'data/audit/conversion_factor_author.json').read_text());xx=np.logspace(-3,3,20000);ii=np.arange(4);target=data['ratio'];out={}
def theta(c,x):
 lw=np.log(c)[:,None]-2*ii[:,None]*np.log(np.atleast_1d(x));p=np.exp(lw-logsumexp(lw,axis=0));return (p*ii[:,None]).sum(0)/3
for kind in ['fit','range']:
 errors=[]
 for pt in data['locus'][kind+'_points']:
  c=[1,3*pt['a'],3*pt['b'],1];yy=theta(c,xx)
  if kind=='fit':
   def f(q):return q[0]+(q[1]-q[0])/(1+(xx/np.exp(q[2]))**q[3])-yy
   p=least_squares(f,[0,1,0,2],gtol=1e-11,xtol=1e-11,ftol=1e-11).x;coef=p[3]
  else:
   low,high=yy[-1],yy[0];roots=[brentq(lambda z:theta(c,np.exp(z))[0]-(low+q*(high-low)),-15,15) for q in [.1,.9]];coef=np.log(81)/(roots[0]-roots[1])
  errors.append(abs(coef/2-target))
 out[kind]={'number_checked':len(errors),'max_ratio_error':max(errors)}
(ROOT/'data/audit/locus_checks.json').write_text(json.dumps(out,indent=2));print(out)
