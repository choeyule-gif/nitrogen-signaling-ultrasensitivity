"""The cascade accounting of R4/R5 (Fig 8), with the digitisation bootstrap.

Point estimate: fit the four-parameter Hill function to the digitised titration, build the
composite of the first layer with the all-or-none readout, pin the second layer's
half-point Khat_2 to each measured cascade half-point, and score everything with the
range measure n_H^90 = log81/log(S10/S90) that the experiments use.  Then repeat in each
bootstrap replicate of the digitisation (+-0.06 ordinate, +-6% abscissa, uniform).
Also reports the supremum of the readout factor over 24 decades of Khat_2 (both readings
of l_2) and the matched factor (composite = (1-theta)^3 exactly), for h = 1..4.
Writes results/cascade_budget.json and results/boot_residuals.npy (used by sensitivity_major1.py).
"""
import _path  # noqa: F401
import os, numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from esbm import RESULTS, save_json
from esbm.data import digitised_titration, constants
from esbm.model import n1090

g, U = digitised_titration()
cc = constants()["jiangninfa2011"]
HALF = cc["cascade_halfpoint_mM"]; MEAS = np.array(cc["cascade_hill_90"]); PII = cc["cascade_PII_uM"]

def fit4(gg, UU):
    def sse(p):
        U0, Um, K, h = p[0], p[1], np.exp(p[2]), p[3]
        return np.sum((Um + (U0-Um)*K/(K+gg**h) - UU)**2)
    best = None
    for s0 in ([2.97, 0.45, np.log(0.31), 2.0], [3.0, 0.4, -1.0, 1.5], [2.9, 0.5, -1.5, 2.5]):
        r = minimize(sse, s0, method='Nelder-Mead', options=dict(xatol=1e-11, fatol=1e-14, maxiter=20000, maxfev=20000))
        if best is None or r.fun < best.fun: best = r
    U0, Um, K, h = best.x[0], best.x[1], np.exp(best.x[2]), best.x[3]
    return U0/3, Um/3, K**(1.0/h), h            # theta_0, theta_inf, S, h

def budget(th0, thi, S, h, reading='lin'):
    theta = lambda x: thi + (th0-thi)/(1+(x/S)**h)
    n1 = n1090(lambda x: (th0-theta(x))/(th0-thi), 1e-12, 1e12)
    def f(x):
        t = theta(x); l1 = (1-t)**3
        return ((1-l1) if reading == 'compl' else t)/l1        # l_2/l_1
    f0, fi = f(1e-12), f(1e12)
    def z(K2): return lambda x: ((f0-f(x))/(f0-fi))*((K2+fi)/(K2+f(x)))
    out = []
    for hp in HALF:
        K2 = np.exp(brentq(lambda lk: np.log(brentq(lambda x: z(np.exp(lk))(x)-0.5, 1e-12, 1e12)/hp), -45, 45))
        nc = n1090(z(K2), 1e-12, 1e12); out.append((nc/n1, nc, K2))
    return n1, out, f, z

th0, thi, S, h = fit4(g, U)
res = dict(fit=dict(theta0=th0, theta_inf=thi, S_mM=S, h=h, U0=3*th0, Umin=3*thi))
print("fit: U0 %.5f Umin %.5f S %.5f h %.5f" % (3*th0, 3*thi, S, h))
for reading in ('lin', 'compl'):
    n1, base, f, z = budget(th0, thi, S, h, reading)
    res[reading] = dict(first_layer=n1, readout=[b[0] for b in base], composite=[b[1] for b in base],
                        K2_pinned=[b[2] for b in base], residual=[m/b[1] for m, b in zip(MEAS, base)])
    print("reading %-5s first layer %.4f" % (reading, n1))
    for p, m, (fac, comp, K2) in zip(PII, MEAS, base):
        print("  PII=%5g uM  Khat2 %.3g  readout %.4f  composite %.4f  residual %.4f" % (p, K2, fac, comp, m/comp))
    # supremum over 24 decades of Khat_2
    fac_of = lambda lk: n1090(z(np.exp(lk)), 1e-12, 1e12)/n1
    grid = np.linspace(-12*np.log(10), 12*np.log(10), 241)
    vals = np.array([fac_of(lk) for lk in grid])
    r = minimize_scalar(lambda lk: -fac_of(lk), bounds=(grid[vals.argmax()]-2, grid[vals.argmax()]+2), method='bounded', options=dict(xatol=1e-10))
    res[reading].update(factor_min_24dec=vals.min(), factor_sup=-r.fun, K2_at_sup=np.exp(r.x))
    print("  factor over 24 decades of Khat2: %.4f .. sup %.5f at Khat2 = %.4g" % (vals.min(), -r.fun, np.exp(r.x)))
# matched factor (composite = (1-theta)^3 exactly) and predicted half-point; h-independence of both
def matched(hh):
    theta = lambda x: thi + (th0-thi)/(1+(x/S)**hh)
    n1 = n1090(lambda x: (th0-theta(x))/(th0-thi), 1e-12, 1e12)
    l1 = lambda x: (1-theta(x))**3; a, b = l1(1e-12), l1(1e12)
    Fm = lambda x: (l1(x)-a)/(b-a)
    return n1090(Fm, 1e-12, 1e12)/n1, brentq(lambda x: Fm(x)-0.5, 1e-12, 1e12)
mf, half_m = matched(h)
res['matched'] = dict(factor=mf, predicted_halfpoint_mM=half_m)
print("matched factor %.5f ; predicted cascade half-point %.3f mM" % (mf, half_m))
sup_h, mat_h = {}, {}
for hh in (1.0, 2.0, 3.0, 4.0):
    theta = lambda x, hh=hh: thi + (th0-thi)/(1+(x/S)**hh)
    n1 = n1090(lambda x: (th0-theta(x))/(th0-thi), 1e-12, 1e12)
    def fl(x, hh=hh):
        t = theta(x); return t/(1-t)**3
    f0, fi = fl(1e-12), fl(1e12)
    fac_of = lambda lk: n1090(lambda x: ((f0-fl(x))/(f0-fi))*((np.exp(lk)+fi)/(np.exp(lk)+fl(x))), 1e-12, 1e12)/n1
    r = minimize_scalar(lambda lk: -fac_of(lk), bounds=(-12*np.log(10), 12*np.log(10)), method='bounded', options=dict(xatol=1e-10))
    sup_h[hh] = -r.fun; mat_h[hh] = matched(hh)[0]
    print("h=%.0f  sup factor %.5f  matched factor %.5f" % (hh, -r.fun, mat_h[hh]))
res['sup_vs_h'] = sup_h; res['matched_vs_h'] = mat_h

# ---------------- bootstrap of the digitisation (600 replicates, seed 20260903)
rng = np.random.default_rng(20260903)
REP = 600; acc = []
for _ in range(REP):
    gg = g*(1+0.06*rng.uniform(-1, 1, g.size))
    UU = np.clip(U + 0.06*rng.uniform(-1, 1, U.size), 1e-3, 2.999)
    try:
        p = fit4(gg, UU)
        if not (0 < p[1] < p[0] < 1 and p[2] > 0 and 0.2 < p[3] < 6): continue
        nn, o, _, _ = budget(*p)
        acc.append([p[3], nn] + [x for fc in o for x in fc[:2]] + [m/fc[1] for m, fc in zip(MEAS, o)])
    except Exception:
        continue
A = np.array(acc); print("\n%d/%d replicates usable" % (len(A), REP))
q = lambda v: (float(np.median(v)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
resid = A[:, 8:11]
mono = ((resid[:, 0] <= resid[:, 1]) & (resid[:, 1] <= resid[:, 2])) | ((resid[:, 0] >= resid[:, 1]) & (resid[:, 1] >= resid[:, 2]))
peak = (resid[:, 1] > resid[:, 0]) & (resid[:, 1] > resid[:, 2])
print("residual peaks at the middle PII in %d/%d replicates; non-monotone in %d" % (peak.sum(), len(resid), (~mono).sum()))
print("h              %.3f [%.3f, %.3f]" % q(A[:, 0]))
print("first layer    %.3f [%.3f, %.3f]" % q(A[:, 1]))
boot = dict(n_usable=len(A), n_total=REP, peak_middle=int(peak.sum()), h=q(A[:, 0]), first_layer=q(A[:, 1]))
for i, p in enumerate(PII):
    print("PII=%5g  readout %.3f [%.3f, %.3f]   composite %.3f [%.3f, %.3f]   residual %.3f [%.3f, %.3f]"
          % ((p,) + q(A[:, 2+2*i]) + q(A[:, 3+2*i]) + q(A[:, 8+i])))
    boot['PII_%g' % p] = dict(readout=q(A[:, 2+2*i]), composite=q(A[:, 3+2*i]), residual=q(A[:, 8+i]))
res['bootstrap'] = boot
np.save(os.path.join(RESULTS, "boot_residuals.npy"), resid)
save_json("cascade_budget", res)
