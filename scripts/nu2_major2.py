"""Major 2: let the second layer's ladder factor nu_2 float.

The composite of R3/R4 puts the ATase-GS ladder at nu_2 = 1 (identical, independent),
the very assumption the Discussion doubts on Hart's avidity result.  Here nu_2 is a free
parameter: the second layer is given a non-binomial Adair ladder on n = 12 sites with
cooperativity w, nu_2 is measured from it, and the composite and residual recomputed.
"""
import _path  # noqa: F401
from esbm import save_json

import numpy as np
from scipy.optimize import brentq
from math import comb

U0, Um, S, h = 2.96771, 0.45350, 0.55954, 2.01789
th0, thi = U0/3, Um/3
theta = lambda g: thi + (th0-thi)/(1+(g/S)**h)
HALF = [0.225, 0.530, 0.600]; MEAS = np.array([4.48, 6.46, 5.23]); PII = [36, 5, 0.5]
n2 = 12

def ladder(w, Kh2):
    """K_i = (n-i+1)/i * Kh2 * w^(i-1);  w = 1 is binomial (nu_2 = 1)."""
    return np.array([(n2-i+1)/i*Kh2*w**(i-1) for i in range(1, n2+1)])

def theta2_and_nu(w, Kh2, u):
    K = ladder(w, Kh2)
    c = np.concatenate([[1.0], np.cumprod(K)])
    p = c*u**np.arange(n2+1)
    p = p/p.sum()
    i = np.arange(n2+1); m = (p*i).sum(); v = (p*(i-m)**2).sum()
    nu = n2*v/(m*(n2-m)) if 0 < m < n2 else np.nan
    return m/n2, nu

def n1090(F, lo=1e-14, hi=1e14):
    a = brentq(lambda x: F(x)-0.1, lo, hi); b = brentq(lambda x: F(x)-0.9, lo, hi)
    return np.log(81)/np.log(b/a)

f1 = lambda g: (1-theta(g))**3
fchain = lambda g: theta(g)/f1(g)                 # l_2 / l_1, linear reading
nH1 = n1090(lambda g: (th0-theta(g))/(th0-thi))

def composite(w, Kh2):
    u = lambda g: 1.0/fchain(g)                   # u = 1/f
    F = lambda g: theta2_and_nu(w, Kh2, u(g))[0]
    a, b = F(1e-12), F(1e12)
    return lambda g: (F(g)-a)/(b-a)

print("%6s %7s   %-28s %-28s" % ("w", "nu_2", "composite (36,5,0.5)", "residual (36,5,0.5)"))
rows = []
for w in (1.0, 1.10, 1.25, 1.5, 2.0, 3.0):
    comps, nus = [], []
    for hp in HALF:
        lk = brentq(lambda t: np.log(brentq(lambda g: composite(w, np.exp(t))(g)-0.5, 1e-10, 1e10)/hp), -40, 40)
        Kh2 = np.exp(lk); C = composite(w, Kh2)
        comps.append(n1090(C, 1e-10, 1e10))
        g50 = brentq(lambda g: C(g)-0.5, 1e-10, 1e10)
        nus.append(theta2_and_nu(w, Kh2, 1.0/fchain(g50))[1])
    comps = np.array(comps); res = MEAS/comps
    rows.append((w, np.mean(nus), comps, res))
    print("%6.2f %7.3f   %-28s %-28s" % (w, np.mean(nus),
          np.array2string(comps, precision=2), np.array2string(res, precision=3)))
save_json("nu2_major2", dict(rows=[dict(w=float(r[0]), nu2=float(r[1]), composite=r[2].tolist(), residual=r[3].tolist()) for r in rows],
                             peak_always=bool(all(r[3][1] > max(r[3][0], r[3][2]) for r in rows))))
print("\nresidual peaks at the middle PII in every case:",
      all(r[3][1] > max(r[3][0], r[3][2]) for r in rows))
w_abs = [r[0] for r in rows if r[3].max() < 1.05]
print("nu_2 that would absorb the residual entirely (max residual < 1.05):", w_abs or "none in this range")
