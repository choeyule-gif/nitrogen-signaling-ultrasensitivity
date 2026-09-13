from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
import sympy as s,json
from pathlib import Path
l,e=s.symbols('l e');n=5;A=s.Rational(1,50);bh=[5,4,3,2,1];gh=[s.Rational(1,50),20,s.Rational(21,2),25600,s.Rational(65,2)];c=[1,5,10,10,5,1];TE,TS,TL=550,960,196
M=sum(c[i]*l**(n-i) for i in range(n+1));W=sum(bh[i]*c[i]*l**(n-i) for i in range(n));V=sum(gh[i]*c[i+1]*l**(n-i) for i in range(n));Q=W+A*V
P=(1+A*l)*Q*e**2+((1+A*l)*M+(TS-TE)*Q)*e-TE*M
G=A*l*Q*e**2+((l-TL)*Q+A*l*M+A*TS*V)*e+(l-TL)*M
res=s.factor(s.resultant(P,G,e));print(res)
factors=s.factor_list(res)[1];out=[]
for f,m in factors:
 p=s.Poly(f,l);cnt=p.count_roots(0,s.oo);out.append(dict(degree=int(p.degree()),multiplicity=int(m),nonnegative_roots=int(cnt),polynomial=str(f)))
(ROOT/'data/audit/exact_root_count.json').write_text(json.dumps(out,indent=2));print(out)
p12=[f for f,m in factors if s.degree(f,l)==12][0];rr=s.polys.polytools.intervals(p12,eps=s.Rational(1,10**25));physical=[]
for (a,b),mult in rr:
 if a<=0:continue
 lv=(a+b)/2; ev=s.N(-s.resultant(P,G,e),5) if False else None
 # Solve the linear combination that cancels e^2.
 polyE=s.Poly(P,e);polyG=s.Poly(G,e);linear=s.expand(polyG.coeff_monomial(e**2)*P-polyE.coeff_monomial(e**2)*G);lin=s.Poly(linear,e)
 er=-lin.coeff_monomial(1).subs(l,lv)/lin.coeff_monomial(e).subs(l,lv)
 physical.append({'l':str(s.N(lv,20)),'e':str(s.N(er,20)),'positive_e':bool(er>0),'l_below_TL':bool(lv<TL)})
(ROOT/'data/audit/exact_physical_roots.json').write_text(json.dumps(physical,indent=2));print(physical)
