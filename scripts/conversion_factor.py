"""Methods, 'Converting a reported coefficient to a local slope'.

Every three-site ladder is fixed, up to the scale symmetry c_i -> lambda^i c_i, by two shape
coordinates.  Three one-parameter families are used: an MWC-type ladder c_i = C(3,i) w^{i(i-1)/2},
a ladder with only the last step strengthened (c_3 -> w^3 c_3), and one with the middle states
suppressed (c_1, c_2 -> c_1/w, c_2/w).  For each family the parameter w is solved so that the
ratio of the trimer's coefficient to the monomer's equals the observed 1.1106, with the
coefficient taken either as the fitted exponent n_H^fit or as the range measure n_H^90, and
the local nu is evaluated at the trimer's own half-point.  The conversion is nu / 1.1106.
Manuscript: 1.026-1.037 (n_H^fit), 1.036-1.059 (n_H^90); quoted 1.04 +- 0.02.
Also checks h-independence of the conversion for h in [0.5, 7].
"""
import _path  # noqa: F401
import numpy as np
from math import comb
from scipy.optimize import curve_fit, brentq
from esbm import save_json
from esbm.model import hill4, n1090
from esbm.data import constants
v = constants()["ventura2010"]
RATIO = v["hill_wildtype_trimer"][0]/v["hill_monovalent_construct"][0]
def lgrid(h): return np.logspace(-6/h, 6/h, 20000)      # +-3 decades at h = 2, scaled so that a pure power is a reparametrisation
l = lgrid(2.0)
def theta(c, l, h=2.0):
    u = l**(-h); n = len(c)-1; w = np.array([c[k]*u**k for k in range(n+1)])
    return (np.arange(n+1)*w).sum()/(n*w.sum())
def nu_local(c, l, h=2.0):
    u = l**(-h); n = len(c)-1; w = np.array([c[k]*u**k for k in range(n+1)]); p = w/w.sum(); i = np.arange(n+1)
    m = (i*p).sum(); vv = ((i-m)**2*p).sum(); return n*vv/(m*(n-m))
def coeff(c, kind, h=2.0):
    l = lgrid(h); th = np.array([theta(c, x, h) for x in l])
    if kind == 'fit':
        p, _ = curve_fit(hill4, l, th, p0=[th[0], th[-1], 1.0, h], maxfev=400000); return p[3]
    F = lambda x: (theta(c, x, h)-th[0])/(th[-1]-th[0])                      # increasing in x
    return n1090(F, 10**(-12/h), 10**(12/h))
def nu_at_half(c, h=2.0):
    l = lgrid(h); th = np.array([theta(c, x, h) for x in l]); ih = np.argmin(abs(th-(th[0]+th[-1])/2))
    return nu_local(c, l[ih], h)
FAM = {"MWC c_i=C(3,i)w^{i(i-1)/2}": lambda w: np.array([comb(3, k)*w**(k*(k-1)/2) for k in range(4)], float),
       "last step only c_3*=w^3": lambda w: np.array([1, 3, 3, 1*w**3], float),
       "middle suppressed c_1,c_2/=w": lambda w: np.array([1, 3/w, 3/w, 1], float)}
out = {}
for kind in ('fit', 'range'):
    mono = coeff(np.array([1.0, 1.0]), kind)
    print("estimator n_H^%s  (monomer coefficient %.4f)" % ('fit' if kind == 'fit' else '90', mono))
    res = {}
    for lbl, f in FAM.items():
        w = brentq(lambda w: coeff(f(w), kind)/mono-RATIO, 1.0001, 8.0, xtol=1e-9)
        nu = nu_at_half(f(w)); res[lbl] = dict(w=w, nu=nu, conversion=nu/RATIO)
        print("   %-32s w = %.4f  nu(half-point) = %.4f  conversion = %.4f" % (lbl, w, nu, nu/RATIO))
    out[kind] = res
    out[kind+'_range'] = [min(r['conversion'] for r in res.values()), max(r['conversion'] for r in res.values())]
print("locus of ladders c = (1, 3a, 3b, 1) reproducing the observed ratio (a on a grid, b solved):")
locus = {}
for kind in ('fit', 'range'):
    mono = coeff(np.array([1.0, 1.0]), kind); conv = []; pts = []
    for a in np.logspace(np.log10(0.25), np.log10(4.0), 25):
        F = lambda b: coeff(np.array([1.0, 3*a, 3*b, 1.0]), kind)/mono-RATIO
        bs = np.logspace(-1.5, 1.5, 31); vals = np.array([F(b) for b in bs])
        for i in range(len(bs)-1):
            if vals[i]*vals[i+1] < 0:
                b = brentq(F, bs[i], bs[i+1], xtol=1e-9); c = np.array([1.0, 3*a, 3*b, 1.0])
                th = np.array([theta(c, x) for x in l]); p, _ = curve_fit(hill4, l, th, p0=[th[0], th[-1], 1.0, 2.0], maxfev=400000)
                rmse = float(np.sqrt(np.mean((hill4(l, *p)-th)**2)))
                conv.append(nu_at_half(c)/RATIO); pts.append(dict(a=float(a), b=float(b), conversion=conv[-1], hill_fit_rmse=rmse))
    locus[kind] = [min(conv), max(conv), len(conv)]
    for q in sorted(pts, key=lambda q: q['conversion'])[::max(1, len(pts)//6)]:
        print("      a=%.3f b=%.3f conversion %.4f  (RMSE of a Hill fit to the trimer curve %.4f)" % (q['a'], q['b'], q['conversion'], q['hill_fit_rmse']))
    locus[kind+'_points'] = pts
    print("   n_H^%s: conversion along the locus in [%.4f, %.4f] (%d ladders)" % ('fit' if kind == 'fit' else '90', min(conv), max(conv), len(conv)))
out['locus'] = locus
print("h-independence (MWC family, n_H^fit estimator):")
f = FAM["MWC c_i=C(3,i)w^{i(i-1)/2}"]; hind = {}
for h in (0.5, 1.0, 2.0, 4.0, 7.0):
    mono = coeff(np.array([1.0, 1.0]), 'fit', h)
    w = brentq(lambda w: coeff(f(w), 'fit', h)/mono-RATIO, 1.0001, 8.0, xtol=1e-9)
    hind[h] = nu_at_half(f(w), h)/RATIO; print("   h=%.1f conversion %.6f" % (h, hind[h]))
out['h_independence'] = hind; out['ratio'] = RATIO
save_json("conversion_factor", out)
