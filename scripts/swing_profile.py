"""R3: the two plateaus measure the product of the two regulations.

Reparametrise the fit by (theta_0, swing s = f_max/w_0, S, h) and profile the likelihood in s.
Manuscript: theta_0 = 0.989, theta_inf = 0.151, w_0/Khat = 0.011, maximum-likelihood swing
about 5e2, 95% interval [175, inf), swing = R_rev = 6.8 rejected at 2dL = 54, hence a
transferase inhibition of at least 26-fold (= 175/6.8)."""
import _path  # noqa: F401
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
from esbm import save_json
from esbm.data import digitised_titration, constants
g, U = digitised_titration(); n = len(U)
Rrev = constants()["jiangninfa2011"]["removing_activity_fold_activation_R_rev"]
def model(p, s):
    th0, S, h = p; w = 1/th0-1; thinf = 1/(1+s*w)
    return 3*(thinf+(th0-thinf)/(1+(g/S)**h))
def sse(p, s):
    if not (0 < p[0] < 1 and p[1] > 0 and p[2] > 0): return 1e9
    r = U-model(p, s); return r@r
def nll(S_): return n/2*np.log(2*np.pi*S_/n)+n/2
def prof(s):
    b = 1e9
    for _ in range(150):
        x0 = [np.random.uniform(0.95, 0.999), np.random.uniform(0.3, 1.0), np.random.uniform(1.5, 2.5)]
        r = minimize(sse, x0, args=(s,), method='Nelder-Mead', options=dict(maxiter=20000, fatol=1e-14))
        b = min(b, r.fun)
    return b
# plateaus of the free (h free) fit, which R3 quotes
from scipy.optimize import curve_fit
from esbm.model import hill4
p4, _ = curve_fit(hill4, g, U, p0=[3, 0.45, 0.55, 2], maxfev=400000)
th0f, thif = p4[0]/3, p4[1]/3; w0f = 1/th0f-1; fmf = 1/thif-1
print("free fit: theta_0 = %.3f, theta_inf = %.3f, w_0/Khat = %.4f, f_max/Khat = %.3f, swing = %.0f" % (th0f, thif, w0f, fmf, fmf/w0f))
np.random.seed(1)
grid = np.logspace(0.6, 4.5, 40)
L = np.array([nll(prof(s)) for s in grid]); d = 2*(L-L.min()); thr = chi2.ppf(0.95, 1)
ins = grid[d <= thr]
smin = grid[np.argmin(d)]
print("minimum at swing = %.4g ; 95%% region [%.4g, %.4g] (grid top %.3g) ; 2dL at the top of the grid %.3f" % (smin, ins.min(), ins.max(), grid[-1], d[-1]))
lower = np.exp(np.interp(thr, d[:np.argmin(d)+1][::-1], np.log(grid[:np.argmin(d)+1])[::-1]))
print("interpolated lower 95%% bound on the swing: %.1f  -> transferase inhibition at least %.1f-fold (R_rev = %.1f)" % (lower, lower/Rrev, Rrev))
pts = {}
for s in (Rrev, 100, 175, 516):
    pts[str(s)] = 2*(nll(prof(s))-L.min()); print("  swing=%6.1f: 2dL = %.2f" % (s, pts[str(s)]))
save_json("swing_profile", dict(theta0_free=float(th0f), theta_inf_free=float(thif), w0_over_Khat=float(w0f), swing_point=float(fmf/w0f), swing_ML=float(smin), lower95=float(lower), stat_at_top=float(d[-1]), closed_above=bool(d[-1] > thr),
                                min_inhibition_fold=float(lower/Rrev), stat=pts))
