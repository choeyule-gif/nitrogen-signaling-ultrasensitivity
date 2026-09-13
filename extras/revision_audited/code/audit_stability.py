from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
from pathlib import Path
import ast,json
import numpy as np
p=ROOT/'code/ref_model.py'
tree=ast.parse(p.read_text());ns={'np':np}
exec(compile(ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef)],type_ignores=[]),str(p),'exec'),ns)
rows=[];maxerr=0;maxred=0
for n in [1,2,3,4,6]:
 L,q,pi,phi,m=ns['chain'](n);d=pi@phi
 vals=[]
 for om in np.logspace(-3,3,300):
  x=np.linalg.solve(1j*om*np.eye(m)-L,phi);F=phi@(pi*(-L@x));old=ns['ReF'](L,q,pi,phi,m,om);maxerr=max(maxerr,abs(F.real-old));vals.append(F.real)
  # Q0 contains departures from all complexes; verify rank-one loop reduction at e=A=l=1.
  Q=L.T; v=-.5*(phi<0);z=pi;Q0=Q.copy();Q0[:,phi>0]=0;u=-Q0@z;g=v@np.linalg.solve(1j*om*np.eye(m)-Q,u)
  maxred=max(maxred,abs(g+F/8))
 rows.append(dict(n=n,F0=float(1-d*d),F_at_min_frequency=vals[0],min_ReF=min(vals),circulation_threshold=float(ns['circ_threshold'](L,q,pi,m))))
(ROOT/'data/audit/ref_checks.json').write_text(json.dumps(dict(max_direct_resolvent_error=maxerr,max_rank1_reduction_error=maxred,rows=rows),indent=2));print(maxerr,maxred,rows)
