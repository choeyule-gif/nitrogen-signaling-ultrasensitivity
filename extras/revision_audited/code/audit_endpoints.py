from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
from pathlib import Path
import numpy as np,json
from scipy.optimize import least_squares
p=ROOT/'data/jiangninfa2011_fig2B_digitised.csv';x,y=np.loadtxt(p,delimiter=',',skiprows=1).T
out=[]
for kind in ['supplied','zero_input','zero_input_y3','drop_first']:
 xx=x.copy();yy=y.copy()
 if kind.startswith('zero'):xx[0]=0
 if kind=='zero_input_y3':yy[0]=3
 if kind=='drop_first':xx=xx[1:];yy=yy[1:]
 def f(p):return p[0]+(p[1]-p[0])/(1+(xx/np.exp(p[2]))**p[3])-yy
 r=least_squares(f,[.45,2.97,np.log(.56),2.],bounds=([0,1.5,-8,.1],[1.5,3,5,8]),gtol=1e-12,xtol=1e-12,ftol=1e-12)
 out.append(dict(case=kind,Umin=r.x[0],Umax=r.x[1],S=np.exp(r.x[2]),h=r.x[3],SSE=float(r.fun@r.fun)))
(ROOT/'data/audit/endpoint_checks.json').write_text(json.dumps(out,indent=2));print(out)
