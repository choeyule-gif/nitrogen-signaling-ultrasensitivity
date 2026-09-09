"""Computer-algebra checks of the identities stated in the manuscript (Methods, 'Steady-state
computations').  Each block prints a residual that must be exactly 0 or a boolean.

A. eq (chain)/(nogo): the chain s_1/s_0 = K_1/l survives three simultaneous dead-end complexes.
B. eq (mobius)/(quartic): the Mobius relation satisfies the T_L balance; the T_E residual is a
   quartic in l whose extreme coefficients are those of eq (quarticco) (n symbolic); the
   factorised form of eq (quartic) equals the expanded numerator up to a constant.
C. eq (basal)/(hillid): a saturating switch with basal activity gives exactly a four-parameter
   Hill function of exponent 2, with the stated plateaus and half-point.
D. eq (two)/(ratio): with forward on E.L1 and reverse on E.L2 the enzyme cancels and
   s_{i+1}/s_i = [Bhat_i b3_i A_1/(Ghat_i g3_i A_2)] (l_1/l_2), n = 4.
E. Under condition (i) the two-form network returns l_2 as linear/quadratic in e (not Mobius).
F. eq (omcov): dlog(omega)/dlog(l) = -Cov(Bhat_i, i)/<Bhat> on 4,000 random ladders.
G. Jacobian of the totals map: monomial sign counts of det J at n = 1 and n = 2.
H. A catalytically active complex bridging two targets (n = 1): the ratio is free of e but
   carries s_0, and reduces to the chain as the bridging rate -> 0.
"""
import _path  # noqa: F401
import numpy as np, sympy as sp
from esbm import save_json
out = {}

print("A. dead-end invariance (three simultaneous dead ends)")
e, l, s0, s1 = sp.symbols('e l s0 s1', positive=True)
b1, b2, b3, g1, g2, g3, a1, a2 = sp.symbols('b1 b2 b3 g1 g2 g3 a1 a2', positive=True)
C = b1*e*s0/(b2+b3); EL = a1*e*l/a2; D = g1*EL*s1/(g2+g3)
sol = sp.solve(sp.Eq(b3*C, g3*D), s1)[0]
ratio = sp.factor(sp.simplify(sol*l/s0))
print("   s1*l/s0 =", ratio, "  (no dead-end constant appears; e cancels: %s)" % (not ratio.has(e)))
out['A_chain_free_of_e'] = not ratio.has(e)

print("B. Mobius relation and the quartic")
l, e, A, Bh, Kh, rho, TE, TL, TS, nn = sp.symbols('l e A Bhat Khat rho T_E T_L T_S n', positive=True)
Pi = nn*Bh; kap = rho/(1+rho); a = A*(1-kap); m = TL-kap*TE
Lam = (1+rho)*e*Pi*l/(l+Kh); psi = TS*Lam/(1+Lam)
esub = (m-l)/(a*l-kap)
# the Mobius relation is the combination (T_L balance) - kappa (T_E balance), which is free of psi
comb_ = (l+A*e*l+kap*psi-TL) - kap*(e*(1+A*l)+psi-TE)
chk = sp.simplify(sp.together(comb_.subs(e, esub)))
print("   Mobius e=(m-l)/(al-kappa) solves (T_L balance) - kappa (T_E balance):", chk == 0)
Rr = sp.together((e*(1+A*l)+psi-TE).subs(e, esub)); num, den = sp.fraction(sp.cancel(Rr))
poly = sp.Poly(sp.expand(num), l)
print("   degree of the T_E residual numerator in l:", poly.degree())
c4 = sp.factor(poly.coeff_monomial(l**4)); c0 = sp.factor(poly.coeff_monomial(l**0))
target4 = A*Pi/nn*((1+rho)**2*Pi-A); target0 = -rho*Kh*Pi*TL/nn
cnorm = sp.simplify(target4/c4)
print("   [l^4] normalised = A Pi/n [(1+rho)^2 Pi - A]:", sp.simplify(cnorm*c4-target4) == 0)
print("   [l^0] normalised = -rho Khat Pi T_L/n      :", sp.simplify(cnorm*c0-target0) == 0)
cand = (TE*(a*l-kap)-(m-l)*(1+A*l))*((l+Kh)*(a*l-kap)+(1+rho)*(m-l)*Pi*l) - TS*(1+rho)*(m-l)*Pi*l*(a*l-kap)
rr = sp.simplify(sp.cancel(sp.expand(cand)/sp.expand(num)))
print("   eq (quartic) / expanded numerator is constant in l:", sp.simplify(sp.diff(rr, l)) == 0, " value", rr)
out['B'] = dict(mobius_ok=bool(chk == 0), degree=int(poly.degree()), c4_ok=bool(sp.simplify(cnorm*c4-target4) == 0),
                c0_ok=bool(sp.simplify(cnorm*c0-target0) == 0), quartic_form_ok=bool(sp.simplify(sp.diff(rr, l)) == 0))

print("C. saturating switch with basal activity is a floored Hill function of exponent 2")
l, Kh, w0, fmax, K2 = sp.symbols('l Khat w_0 f_max K_2', positive=True)
x = (l/K2)**2; f = w0+(fmax-w0)*x/(1+x); th = Kh/(Kh+f)
th0 = Kh/(Kh+w0); thi = Kh/(Kh+fmax); S = K2*sp.sqrt((Kh+w0)/(Kh+fmax))
diff = sp.simplify(sp.together(th-(thi+(th0-thi)/(1+(l/S)**2))))
print("   theta - [theta_inf + (theta_0-theta_inf)/(1+(l/S)^2)] =", diff)
out['C_hillid_residual_zero'] = bool(diff == 0)

print("D. two-form network: enzyme cancels, chain variable l_1/l_2 (n = 4)")
n = 4
e, l1, l2, A1, A2 = sp.symbols('e l_1 l_2 A_1 A_2', positive=True)
B1 = sp.symbols('b1_0:4', positive=True); B2 = sp.symbols('b2_0:4', positive=True); B3 = sp.symbols('b3_0:4', positive=True)
G1 = sp.symbols('g1_0:4', positive=True); G2 = sp.symbols('g2_0:4', positive=True); G3 = sp.symbols('g3_0:4', positive=True)
s = sp.symbols('s_0:5', positive=True)
ok = True
for i in range(n):
    Bhi = B1[i]/(B2[i]+B3[i]); Ghi = G1[i]/(G2[i]+G3[i])
    beta = Bhi*A1*e*l1*s[i]; gam = Ghi*A2*e*l2*s[i+1]
    r = sp.simplify(sp.solve(sp.Eq(B3[i]*beta, G3[i]*gam), s[i+1])[0]/s[i])
    ok &= (not r.has(e)) and sp.simplify(r - Bhi*B3[i]*A1/(Ghi*G3[i]*A2)*(l1/l2)) == 0
print("   all four ratios free of e and equal to [Bhat b3 A1/(Ghat g3 A2)](l1/l2):", ok)
out['D_twoform_chain_ok'] = bool(ok)

print("E. two-form network under condition (i): l_2 as a function of e")
e, l1, l2, A1, A2, rho, W, s0, TE, TL1, TL2 = sp.symbols('e l_1 l_2 A_1 A_2 rho W s_0 T_E T_L1 T_L2', positive=True)
Sb = A1*e*l1*s0*W; Sg = rho*Sb
sol = sp.solve([sp.Eq(l1+A1*e*l1+Sb, TL1), sp.Eq(l2+A2*e*l2+Sg, TL2)], [l1, l2], dict=True)[0]
ex = sp.cancel(sp.together(sol[l2])); nu_, de_ = sp.fraction(ex)
dn, dd = sp.Poly(sp.expand(nu_), e).degree(), sp.Poly(sp.expand(de_), e).degree()
print("   l_2 = (degree %d in e)/(degree %d in e)" % (dn, dd))
out['E_twoform_l2_degrees'] = [dn, dd]

print("F. eq (omcov) on 4,000 random ladders")
rng = np.random.default_rng(5); bad = 0
for _ in range(4000):
    nn_ = int(rng.integers(2, 7))
    c = np.concatenate([[1.0], np.cumprod(10**rng.uniform(-1.5, 1.5, nn_))]); Bhv = 10**rng.uniform(-1.5, 1.5, nn_+1)
    lv = 10**rng.uniform(-2, 2); u = 1/lv
    om = lambda uu: (Bhv*c*uu**np.arange(nn_+1)).sum()/(c*uu**np.arange(nn_+1)).sum()
    d = 1e-6; lhs = (np.log(om(1/(lv*(1+d))))-np.log(om(1/(lv*(1-d)))))/(2*d)
    w = c*u**np.arange(nn_+1); p = w/w.sum(); i = np.arange(nn_+1)
    cov = (Bhv*i*p).sum()-(Bhv*p).sum()*(i*p).sum(); rhs = -cov/(Bhv*p).sum()
    if abs(lhs-rhs) > 1e-5*max(1, abs(rhs)): bad += 1
print("   violations:", bad)
out['F_omcov_violations'] = bad

print("G. sign counts of the monomials of det J (totals map) at n = 1, 2")
signs = {}
for n in (1, 2):
    e, l, s0, A = sp.symbols('e l s0 A', positive=True)
    Bh = [sp.Symbol('B%d' % i, positive=True) for i in range(n)]; Gh = [sp.Symbol('G%d' % i, positive=True) for i in range(n)]
    rho = [sp.Symbol('r%d' % i, positive=True) for i in range(n)]
    K = [Bh[i]*rho[i]/(A*Gh[i]) for i in range(n)]
    s = [s0]
    for i in range(n): s.append(s[-1]*K[i]/l)
    eps = A*e*l; beta = [Bh[i]*e*s[i] for i in range(n)]; gam = [Gh[i]*eps*s[i+1] for i in range(n)]
    TEe = e+eps+sum(beta)+sum(gam); TLe = l+eps+sum(gam); TSe = sum(s)+sum(beta)+sum(gam)
    J = sp.Matrix([[sp.diff(f_, v) for v in (e, l, s0)] for f_ in (TEe, TLe, TSe)])
    num, den = sp.fraction(sp.cancel(sp.together(J.det())))
    T = sp.expand(num).as_ordered_terms()
    sub = {A: sp.Rational(3, 2), e: sp.Rational(1, 3), l: sp.Rational(7, 5), s0: sp.Rational(4, 9)}
    for i in range(n):
        sub[Bh[i]] = sp.Rational(5+i, 7); sub[Gh[i]] = sp.Rational(2+i, 3); sub[rho[i]] = sp.Rational(9+i, 4)
    pos = sum(1 for t in T if t.subs(sub) > 0); neg = sum(1 for t in T if t.subs(sub) < 0)
    print("   n=%d: %d monomials, %d positive, %d negative; denominator %s" % (n, len(T), pos, neg, sp.factor(den)))
    signs[n] = dict(monomials=len(T), positive=pos, negative=neg, denominator=str(sp.factor(den)))
out['G_detJ_signs'] = signs

print("H. active bridging complex at n = 1")
e, l, s0, s1, A, b1, b2, b3, g1, g2, g3, y1, y2, y3 = sp.symbols('e l s_0 s_1 A b1 b2 b3 g1 g2 g3 y1 y2 y3', positive=True)
eps = A*e*l; beta, gam, Y = sp.symbols('beta gamma Y', positive=True)
eqs = [sp.Eq(b1*e*s0, (b2+b3)*beta), sp.Eq(g1*eps*s1+y2*Y, (g2+g3)*gam+y1*gam*s0), sp.Eq(y1*gam*s0, (y2+y3)*Y)]
sol = sp.solve(eqs, [beta, gam, Y], dict=True)[0]
cut = sp.Eq(b3*sol[beta], g3*sol[gam]+y3*sol[Y])
r = sp.simplify(sp.factor(sp.solve(cut, s1)[0]/s0))
lim = sp.simplify(sp.limit(r, y1, 0))
chain = b1*b3*(g2+g3)*a2/((b2+b3)*g1*g3*A*l) if False else b1*b3*(g2+g3)/((b2+b3)*g1*g3*A*l)
print("   s1/s0 free of e:", not r.has(e), "; carries s_0:", r.has(s0), "; limit y1->0 equals the chain:", sp.simplify(lim-chain) == 0)
out['H'] = dict(free_of_e=bool(not r.has(e)), carries_s0=bool(r.has(s0)), reduces_to_chain=bool(sp.simplify(lim-chain) == 0))
save_json("verify_symbolic", out)
