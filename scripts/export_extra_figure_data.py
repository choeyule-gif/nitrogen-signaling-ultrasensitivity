"""Data for S2 Fig (the supremum of the ladder factor over every three-site ladder shape, with
the locus of ladders reproducing the observed ratio 1.1106), S4 Fig (reported/calculated cascade
coefficient against the downstream ladder factor at each PII), the measured hyperbola of Fig 2A,
and the three fixed-exponent fits of Fig 3A drawn as curves.  Run after export_figure_data.py."""
import os, sys, json, csv
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
import numpy as np
from math import comb
from scipy.optimize import brentq
from esbm import ROOT
from esbm.model import n1090
from esbm.data import digitised_titration
from scipy.optimize import least_squares
FD = os.path.join(ROOT, "figdata")
os.makedirs(FD, exist_ok=True)
def csvw(name, cols, comment):
    keys = list(cols); rows = zip(*[np.asarray(cols[k]).ravel() for k in keys])
    with open(os.path.join(FD, name), "w", newline="") as f:
        f.write("# " + comment + "\n"); w = csv.writer(f); w.writerow(keys); w.writerows(rows)
# ---- Fig 2A: the hyperbola nu*eta = 2.018 of the fitted first-layer coefficient
nums = json.load(open(os.path.join(ROOT, "results", "revision_numbers.json"))); hfit = nums['fit']['h']
nu = np.linspace(.55, 3, 300)
csvw("Fig2A_measured_hyperbola.csv", dict(nu=nu, eta=hfit/nu), "Fig 3A: nu*eta = %.4f, the fitted first-layer coefficient; constrains the product only" % hfit)
# ---- S2 Fig: sup_u nu over three-site ladder shapes c = (1, 3a, 3b, 1)
def nu_of(c, u):
    i = np.arange(len(c)); w = c*u**i; p = w/w.sum(); m = (i*p).sum(); v = ((i-m)**2*p).sum()
    n = len(c)-1; return n*v/(m*(n-m))
la = np.linspace(-3, 3, 61); lb = np.linspace(-3, 3, 61); lu = np.linspace(-12, 12, 481)
Z = np.empty((len(lb), len(la)))
for j, b in enumerate(10**lb):
    for i, a in enumerate(10**la):
        c = np.array([1.0, 3*a, 3*b, 1.0]); Z[j, i] = max(nu_of(c, 10**x) for x in lu)
LA, LB = np.meshgrid(la, lb)
csvw("FigS2_nu_sup_surface.csv", dict(log10_a=LA.ravel(), log10_b=LB.ravel(), nu_sup=Z.ravel()), "S2 Fig: sup over the input of the ladder factor nu for the three-site ladder c = (1, 3a, 3b, 1); a = b = 1 is the binomial ladder (nu = 1); the bound nu <= 3 is approached as a, b -> 0")
conv = json.load(open(os.path.join(ROOT, "results", "conversion_factor.json")))
for kind in ("fit", "range"):
    pts = sorted(conv['locus'][kind+'_points'], key=lambda q: q['a'])
    zs = [max(nu_of(np.array([1.0, 3*q['a'], 3*q['b'], 1.0]), 10**x) for x in lu) for q in pts]
    csvw("FigS2_locus_%s.csv" % kind, dict(log10_a=[np.log10(q['a']) for q in pts], log10_b=[np.log10(q['b']) for q in pts], nu_sup=zs, conversion=[q['conversion'] for q in pts]),
         "S2 Fig: three-site ladders reproducing the observed ratio 1.1106 with the %s estimator" % ("fitted-exponent" if kind == "fit" else "range"))
print("S2: nu_sup range %.3f .. %.3f ; binomial point %.4f" % (Z.min(), Z.max(), Z[np.argmin(abs(lb)), np.argmin(abs(la))]))
# ---- S4 Fig: reported/calculated against the downstream ladder factor nu_2 (same construction as scripts/nu2_major2.py)
U0, Um, S, h = nums['fit']['Umax'], nums['fit']['Umin'], nums['fit']['S_mM'], hfit
th0, thi = U0/3, Um/3
theta = lambda g: thi + (th0-thi)/(1+(g/S)**h)
HALF = [0.225, 0.530, 0.600]; MEAS = np.array([4.48, 6.46, 5.23]); PII = [36, 5, 0.5]; n2 = 12
def ladder(w, Kh2): return np.array([(n2-i+1)/i*Kh2*w**(i-1) for i in range(1, n2+1)])
def theta2_and_nu(w, Kh2, u):
    K = ladder(w, Kh2); c = np.concatenate([[1.0], np.cumprod(K)]); p = c*u**np.arange(n2+1); p = p/p.sum()
    i = np.arange(n2+1); m = (p*i).sum(); v = (p*(i-m)**2).sum()
    return m/n2, n2*v/(m*(n2-m))
fchain = lambda g: theta(g)/(1-theta(g))**3
def composite(w, Kh2):
    F = lambda g: theta2_and_nu(w, Kh2, 1.0/fchain(g))[0]; a, b = F(1e-12), F(1e12)
    return lambda g: (F(g)-a)/(b-a)
ws = np.concatenate([np.linspace(1.0, 1.5, 21), np.linspace(1.55, 3.0, 15)])
rows = {"w": [], "nu2": []}
for p_ in PII: rows["ratio_PII_%g" % p_] = []
for w in ws:
    nus = []
    for hp, m_ in zip(HALF, MEAS):
        def g_half(t):                       # half-point of the composite, NaN where the ladder is saturated at both ends
            C = composite(w, np.exp(t))
            try: return brentq(lambda g: C(g)-0.5, 1e-10, 1e10)
            except ValueError: return np.nan
        tg = np.linspace(-30, 30, 121); vals = np.array([np.log(g_half(t)/hp) for t in tg]); ok = np.isfinite(vals)
        idx = [i for i in range(len(tg)-1) if ok[i] and ok[i+1] and vals[i]*vals[i+1] < 0][0]
        lk = brentq(lambda t: np.log(g_half(t)/hp), tg[idx], tg[idx+1])
        Kh2 = np.exp(lk); C = composite(w, Kh2); comp = n1090(C, 1e-10, 1e10)
        g50 = brentq(lambda g: C(g)-0.5, 1e-10, 1e10); nus.append(theta2_and_nu(w, Kh2, 1.0/fchain(g50))[1])
        rows["ratio_PII_%g" % ([36, 5, 0.5][HALF.index(hp)])].append(m_/comp)
    rows["w"].append(w); rows["nu2"].append(float(np.mean(nus)))
csvw("FigS4_ratio_vs_nu2.csv", rows, "S4 Fig: reported/calculated cascade coefficient against the downstream ladder factor nu_2 (twelve-site ladder K_i = (n-i+1)/i Khat2 w^(i-1)), Khat2 pinned to each reported midpoint")
print("S4: nu2 from %.3f to %.3f; ratios at nu2~1.34: %s" % (rows['nu2'][0], rows['nu2'][-1], [round(rows['ratio_PII_%g' % p_][np.argmin(abs(np.array(rows['nu2'])-1.34))], 3) for p_ in PII]))

# ---- Fig 4A: the bounded fixed-exponent fits as curves, on the grid of Fig3A_fit.csv
xt, yt = digitised_titration()
BOUNDS = ([0, 1.5, -8, .1], [1.5, 3, 5, 8])
def hill(xx, q):
    lo, hi, ls, h = q
    return lo + (hi-lo)/(1+(np.asarray(xx)/np.exp(ls))**h)
def fit_fixed(h):
    r = least_squares(lambda q: hill(xt, [*q, h])-yt, [.45, 2.97, np.log(.56)],
                      bounds=(np.array(BOUNDS[0])[:3], np.array(BOUNDS[1])[:3]),
                      xtol=1e-10, ftol=1e-10, gtol=1e-10)
    return np.r_[r.x, h]
xx = np.logspace(-2, 1.25, 350)
csvw("Fig3A_fixed_h_fits.csv", dict(glutamine_mM=xx, **{"fit_h%d" % h: hill(xx, fit_fixed(h)) for h in (1, 2, 3)}),
     "Fig 4A: the bounded fixed-exponent fits h = 1, 2, 3 as curves, on the grid of Fig3A_fit.csv")
