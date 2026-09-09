"""Adversarial search for a violation of monotonicity of T_L in the free effector.

Multistationarity of the effector-switched network is equivalent to T_L failing to be
monotone in l at fixed (T_E, T_S).  Conditions (i)-(ii) are proved only for a 1:1
effector complex; this searches the cooperative extension eps = A e l^h, over the
identical-site family (A, Bhat, Ghat, rho, T_E, T_S), each spanning 10^{+-3}, for each
(n, h).  We MINIMISE dlogT_L/dlogl: a negative value would refute the conjecture.
"""
import _path  # noqa: F401
import os, numpy as np, sys
from esbm import RESULTS, save_json
from scipy.optimize import brentq, minimize
np.seterr(all='ignore')

def dlogTL(n, h, A, Bh, Gh, rho, TE, TS, l):
    Kh = Bh*rho/(A*Gh); kap = rho/(1+rho)
    def TL(l):
        u = Kh/l**h
        om = n*Bh/(1+u)
        if not np.isfinite(om) or om <= 0: return np.nan
        def f(e):
            v = e*(1+A*l**h) + TS*((1+rho)*e*om)/(1+(1+rho)*e*om) - TE
            return v if np.isfinite(v) else np.inf
        lo, hi = 0.0, 1.0
        if not np.isfinite(f(lo)): return np.nan
        while f(hi) < 0:
            hi *= 10
            if hi > 1e18: return np.nan
        e = brentq(f, lo, hi, xtol=1e-300, rtol=8.9e-16)
        Lam = (1+rho)*e*om
        return l + A*e*l**h + kap*TS*Lam/(1+Lam)
    d = 1e-6
    a, b = TL(l*(1-d)), TL(l*(1+d))
    if not (np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0): return np.nan
    return (np.log(b)-np.log(a))/(2*d)

def obj(th, n, h):
    A, Bh, Gh, rho, TE, TS, l = [10.0**x for x in th]
    v = dlogTL(n, h, A, Bh, Gh, rho, TE, TS, l)
    return 1e9 if not np.isfinite(v) else v

RESTARTS, MAXFEV = 60, 2500
rng = np.random.default_rng(20260903)
print("n   h   restarts  evaluations   min dlogT_L/dlogl      violations")
worst = {}
for n in (1, 2, 3, 4, 6, 8, 12):
    for h in (1, 2, 3):
        best, nev, viol = np.inf, 0, 0
        for _ in range(RESTARTS):
            th0 = np.concatenate([rng.uniform(-3, 3, 6), rng.uniform(-2, 2, 1)])
            r = minimize(obj, th0, args=(n, h), method='Nelder-Mead',
                         options=dict(maxiter=MAXFEV, maxfev=MAXFEV, xatol=1e-11, fatol=1e-15))
            nev += r.nfev
            if r.fun < best: best, bx = r.fun, r.x
            if r.fun <= 0: viol += 1
        worst[(n, h)] = (best, bx)
        print(f"{n:<3} {h:<3} {RESTARTS:>8} {nev:>12} {best:>21.6e} {viol:>15}")
        sys.stdout.flush()
np.save(os.path.join(RESULTS, 'mono_adv.npy'), np.array([[k[0], k[1], v[0]] for k, v in worst.items()]))
print("\nsmallest over all (n,h): %.6e" % min(v[0] for v in worst.values()))
for h in (1, 2, 3):
    print("  h=%d: smallest over n = %.3e" % (h, min(v[0] for k, v in worst.items() if k[1] == h)))

# ---- re-evaluate every optimum in 50-digit arithmetic -----------------------
import mpmath as mp
mp.mp.dps = 50
def dlogTL_mp(n, h, A, Bh, Gh, rho, TE, TS, l):
    A,Bh,Gh,rho,TE,TS,l = map(mp.mpf,(A,Bh,Gh,rho,TE,TS,l))
    h = mp.mpf(h); Kh = Bh*rho/(A*Gh); kap = rho/(1+rho)
    def TL(l):
        om = n*Bh/(1+Kh/l**h)
        f = lambda e: e*(1+A*l**h) + TS*((1+rho)*e*om)/(1+(1+rho)*e*om) - TE
        hi = mp.mpf(1)
        while f(hi) < 0: hi *= 10
        lo = mp.mpf(0)
        for _ in range(400):
            mid = (lo+hi)/2
            if f(mid) < 0: lo = mid
            else: hi = mid
        e = (lo+hi)/2
        Lam = (1+rho)*e*om
        return l + A*e*l**h + kap*TS*Lam/(1+Lam)
    d = mp.mpf(10)**-25
    return (mp.log(TL(l*(1+d)))-mp.log(TL(l*(1-d))))/(2*d)

print("\n50-digit re-evaluation of each optimum")
print("n   h    d log T_L / d log l (50 dps)")
sm = {}
for (n, h), (v, x) in sorted(worst.items()):
    pars = [10.0**t for t in x]
    r = dlogTL_mp(n, h, *pars)
    sm.setdefault(h, []).append(r)
    print(f"{n:<3} {h:<3} {mp.nstr(r, 10):>26}   {'VIOLATION' if r <= 0 else ''}")
print("\nsmallest verified value, by h:")
for h in sorted(sm): print("  h=%d : %s" % (h, mp.nstr(min(sm[h]), 4)))
print("smallest overall: %s" % mp.nstr(min(min(v) for v in sm.values()), 4))

save_json("mono_adversarial", dict(double_precision_min=float(min(v[0] for v in worst.values())),
                                   verified_min_by_h={str(h): float(min(sm[h])) for h in sm},
                                   verified_min=float(min(min(v) for v in sm.values())),
                                   verified_max=float(max(max(v) for v in sm.values())),
                                   n_pairs=len(worst), violations=int(sum(1 for v in sm.values() for r in v if r <= 0))))
