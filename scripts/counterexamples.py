"""Counterexamples to the two conditions of eq (cond), each with three positive steady states.

(i)  A ladder whose catalytic ratios rho_i alternate by a factor 2.89 (condition (i) broken,
     (ii) intact), found by a seeded random search at n in {2,3,4}; this is Fig 6b.
(ii) Strongly cooperative ladders (condition (ii) broken, (i) intact), found by a seeded
     random search at n in {3..6}; three are kept (S1 Appendix).
Every root is refined and re-verified in 60-digit arithmetic: residual, positivity of every
species, and the elasticity dlog(omega)/dlog(l) at the root.
Writes results/counterexamples_parameters.json (rate constants, totals, roots) and results/counterexamples.json.
"""
import _path  # noqa: F401
import os, json, numpy as np
from scipy.optimize import brentq
from mpmath import mp
from esbm import RESULTS, save_json

# ============================================================ (i): spread 2.89 in rho_i
def resid_i(n, A, c, Bh, Gh, TE, TS, TL, ls):
    """Residual of the T_L balance along the variety parametrised by l (vectorised in l)."""
    u = 1.0/ls; ar = np.arange(n+1)[:, None]; ar2 = np.arange(n)[:, None]
    Mp = (c[:, None]*u**ar).sum(0); Vt = ((Gh*c[1:])[:, None]*u**ar2).sum(0)
    Q = ((Bh*c[:n])[:, None]*u**ar2).sum(0)+A*Vt
    a = (1+A*ls)*Q; b = (1+A*ls)*Mp+Q*TS-TE*Q; cc = -TE*Mp
    e = (-b+np.sqrt(np.maximum(b*b-4*a*cc, 0)))/(2*a); s0 = TS/(Mp+e*Q)
    return ls+A*e*ls+A*e*s0*Vt-TL, e, s0
rng = np.random.default_rng(4242)
ls = np.logspace(-4, 4, 3000); hit = None; trials = 0
while hit is None:
    trials += 1
    n = int(rng.integers(2, 5)); A = 10**rng.uniform(-.8, .8); base = 10**rng.uniform(-1.5, 1.5)
    rho = base*np.array([2.89**(i % 2) for i in range(n)])
    Bh = 10**rng.uniform(-1.5, 1.5, n); Gh = 10**rng.uniform(-1.5, 1.5, n)
    c = np.concatenate([[1.0], np.cumprod(rho*Bh/(A*Gh))])
    TE = 10**rng.uniform(-1.2, 1.2); TS = 10**rng.uniform(-1, 2); TL = 10**rng.uniform(-1.2, 1.6)
    Rv, e, s0 = resid_i(n, A, c, Bh, Gh, TE, TS, TL, ls)
    ok = np.isfinite(Rv) & (e > 0) & (s0 > 0)
    if ok.sum() > 10 and (np.diff(np.sign(Rv[ok])) != 0).sum() >= 3:
        hit = (n, A, c, Bh, Gh, TE, TS, TL, rho)
n, A, c, Bh, Gh, TE, TS, TL, rho = hit
print("(i) spread-2.89 example after %d trials: n=%d rho=%s (max/min %.2f)" % (trials, n, np.array2string(rho, precision=5), rho.max()/rho.min()))
Rv, e, s0 = resid_i(n, A, c, Bh, Gh, TE, TS, TL, ls); ok = np.isfinite(Rv) & (e > 0) & (s0 > 0)
xs, vs = ls[ok], Rv[ok]; idx = np.where(np.diff(np.sign(vs)) != 0)[0]
roots_i = []
for i in idx:
    a_, b_ = xs[i], xs[i+1]
    for _ in range(200):
        m = np.sqrt(a_*b_)
        if np.sign(resid_i(n, A, c, Bh, Gh, TE, TS, TL, np.array([m]))[0][0]) == np.sign(vs[i]): a_ = m
        else: b_ = m
    roots_i.append(np.sqrt(a_*b_))
mp.dps = 60
def R_i_mp(L):
    L = mp.mpf(L); u = 1/L
    Mp = sum(mp.mpf(float(c[k]))*u**k for k in range(n+1))
    Vt = sum(mp.mpf(float(Gh[k]*c[k+1]))*u**k for k in range(n))
    Q = sum(mp.mpf(float(Bh[k]*c[k]))*u**k for k in range(n))+mp.mpf(float(A))*Vt
    aa = (1+mp.mpf(float(A))*L)*Q; bb = (1+mp.mpf(float(A))*L)*Mp+Q*mp.mpf(float(TS))-mp.mpf(float(TE))*Q
    ccc = -mp.mpf(float(TE))*Mp
    ee = (-bb+mp.sqrt(bb*bb-4*aa*ccc))/(2*aa); ss = mp.mpf(float(TS))/(Mp+ee*Q)
    return L+mp.mpf(float(A))*ee*L+mp.mpf(float(A))*ee*ss*Vt-mp.mpf(float(TL)), ee, ss
def bisect_mp(F, r, rel=1e-3, iters=230):
    """Bisection in mp arithmetic on a bracket [r(1-rel), r(1+rel)] widened until F changes sign."""
    a, b = mp.mpf(r)*(1-rel), mp.mpf(r)*(1+rel); fa, fb = F(a), F(b)
    while mp.sign(fa) == mp.sign(fb):
        a /= 1+rel; b *= 1+rel; fa, fb = F(a), F(b)
    for _ in range(iters):
        m = (a+b)/2; fm = F(m)
        if mp.sign(fm) == mp.sign(fa): a, fa = m, fm
        else: b, fb = m, fm
    return (a+b)/2
ver_i = []
for r in roots_i:
    r60 = bisect_mp(lambda L: R_i_mp(L)[0], r)
    v, ee, ss = R_i_mp(r60)
    rr = [mp.mpf(1)]
    for k in range(n): rr.append(rr[-1]*mp.mpf(float(c[k+1]/c[k]))/r60)
    species = [ss*x for x in rr]+[ee]
    ver_i.append(dict(l=mp.nstr(r60, 25), residual=mp.nstr(abs(v), 3), e=mp.nstr(ee, 10), s0=mp.nstr(ss, 10), all_positive=bool(min(species) > 0)))
    print("    root l = %s  |R| = %s  e = %s  s0 = %s  positive %s" % (mp.nstr(r60, 12), mp.nstr(abs(v), 3), mp.nstr(ee, 6), mp.nstr(ss, 6), min(species) > 0))
ex_i = dict(condition_broken="(i)", n=n, A=A, rho=rho.tolist(), Bhat=Bh.tolist(), Ghat=Gh.tolist(), c=c.tolist(),
            K=(c[1:]/c[:-1]).tolist(), TE=TE, TS=TS, TL=TL, roots_double=roots_i, roots_60digit=ver_i, spread=rho.max()/rho.min(), trials=trials)

# ============================================================ (ii): strongly cooperative ladders
def make(n, A, Bh, K, rho, TE, TL, TS):
    kap = rho/(1+rho)
    def MW(l):
        r = np.empty(n+1); r[0] = 1.0
        for i in range(1, n+1): r[i] = r[i-1]*K[i-1]/l
        return r.sum(), (Bh*r).sum()
    def psi(e, l):
        M, W = MW(l); Lam = (1+rho)*e*W/M; return TS*Lam/(1+Lam)
    def e_of(l):
        f = lambda e: l+A*e*l+kap*psi(e, l)-TL
        if f(0.0) > 0: return None
        hi = 1.0
        while f(hi) < 0:
            hi *= 10
            if hi > 1e14: return None
        return brentq(f, 0.0, hi, xtol=1e-300, rtol=8.9e-16)
    def R(l):
        e = e_of(l)
        return np.nan if e is None else TE-e*(1+A*l)-psi(e, l)
    def elast(l):
        h = 1e-6; M1, W1 = MW(l*(1-h)); M2, W2 = MW(l*(1+h))
        return (np.log(W2/M2)-np.log(W1/M1))/(2*h)
    return R, elast
rng = np.random.default_rng(7)
lgrid = np.logspace(-9, 9, 1400); found = []; trials2 = 0
while len(found) < 3 and trials2 < 300000:
    trials2 += 1
    n = int(rng.integers(3, 7))
    A = 10**rng.uniform(-1, 3); rho = 10**rng.uniform(-2, 2)
    Bh = np.zeros(n+1); Bh[:n] = 10**rng.uniform(-2, 1, size=n)
    s = 10**rng.uniform(1, 4)
    K = np.array([10**rng.uniform(-4, -1)]+[10**rng.uniform(-1, 1) for _ in range(n-2)]+[s])
    TE = 10**rng.uniform(-2, 3); TL = 10**rng.uniform(-2, 3); TS = 10**rng.uniform(-2, 3)
    R, el = make(n, A, Bh, K, rho, TE, TL, TS)
    try: es = np.array([el(x) for x in np.logspace(-6, 6, 80)])
    except Exception: continue
    mx = np.nanmax(es)
    if not (mx > 1.0 or np.nanmin(es) < 0.0): continue
    try: v = np.array([R(x) for x in lgrid])
    except Exception: continue
    ok = np.isfinite(v); v2, x2 = v[ok], lgrid[ok]
    if int(np.sum(v2[:-1]*v2[1:] < 0)) >= 3:
        roots = [brentq(R, x2[i], x2[i+1], xtol=1e-300, rtol=8.9e-16) for i in range(len(v2)-1) if v2[i]*v2[i+1] < 0]
        found.append(dict(n=n, A=A, rho=rho, Bh=Bh, K=K, TE=TE, TL=TL, TS=TS, roots=roots, maxel=mx))
        print("(ii) found after %d trials: n=%d roots=%s max elasticity %.3f" % (trials2, n, ["%.6g" % r for r in roots], mx))
ex_ii = []
for d in found:
    n = d['n']; A = mp.mpf(d['A']); rho = mp.mpf(d['rho']); Bh = [mp.mpf(x) for x in d['Bh']]
    K = [mp.mpf(x) for x in d['K']]; TE = mp.mpf(d['TE']); TL = mp.mpf(d['TL']); TS = mp.mpf(d['TS'])
    kap = rho/(1+rho)
    def MW(l):
        r = [mp.mpf(1)]
        for i in range(1, n+1): r.append(r[-1]*K[i-1]/l)
        return sum(r), sum(Bh[i]*r[i] for i in range(n+1)), r
    def psi(e, l):
        M, W, _ = MW(l); Lam = (1+rho)*e*W/M; return TS*Lam/(1+Lam)
    def e_of(l):
        f = lambda e: l+A*e*l+kap*psi(e, l)-TL
        lo, hi = mp.mpf(0), mp.mpf(1)
        while f(hi) < 0: hi *= 10
        for _ in range(230):
            mid = (lo+hi)/2
            if f(mid) < 0: lo = mid
            else: hi = mid
        return (lo+hi)/2
    def R(l):
        e = e_of(l); return TE-e*(1+A*l)-psi(e, l)
    def el(l):
        hh = mp.mpf(10)**-20
        M1, W1, _ = MW(l*(1-hh)); M2, W2, _ = MW(l*(1+hh))
        return (mp.log(W2/M2)-mp.log(W1/M1))/(2*hh)
    ver = []
    for r0 in d['roots']:
        l = bisect_mp(R, r0)
        e = e_of(l); M, W, rr = MW(l); s0 = TS/(M+(1+rho)*e*W); sp = [s0*x for x in rr]
        pos = bool(e > 0 and l > 0 and all(x > 0 for x in sp))
        ver.append(dict(l=mp.nstr(l, 25), residual=mp.nstr(abs(R(l)), 3), elasticity_at_root=mp.nstr(el(l), 8), all_positive=pos))
        print("    n=%d root l=%s |R|=%s elast=%s positive %s" % (n, mp.nstr(l, 12), mp.nstr(abs(R(l)), 3), mp.nstr(el(l), 6), pos))
    ee = [el(mp.mpf(10)**(mp.mpf(k)/4)) for k in range(-24, 25)]
    print("    max dlog omega/dlog l over 12 decades: %s" % mp.nstr(max(ee), 8))
    ex_ii.append(dict(condition_broken="(ii)", n=n, A=d['A'], rho=d['rho'], Bhat=d['Bh'].tolist(), K=d['K'].tolist(), TE=d['TE'], TL=d['TL'], TS=d['TS'],
                      roots_double=d['roots'], roots_60digit=ver, max_elasticity=float(max(ee))))
with open(os.path.join(RESULTS, "counterexamples_parameters.json"), "w") as f:
    json.dump(dict(condition_i=ex_i, condition_ii=ex_ii), f, indent=1, default=float)
save_json("counterexamples", dict(cond_i=dict(n=ex_i['n'], spread=ex_i['spread'], n_roots=len(roots_i), all_positive=all(v['all_positive'] for v in ver_i)),
                                  cond_ii=[dict(n=x['n'], n_roots=len(x['roots_double']), max_elasticity=x['max_elasticity'], all_positive=all(v['all_positive'] for v in x['roots_60digit'])) for x in ex_ii],
                                  trials_i=trials, trials_ii=trials2))
print("written results/counterexamples_parameters.json")
