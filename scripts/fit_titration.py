"""Reanalysis of the digitised titration (R2 and Methods, 'Data reanalysis and model selection').

1. Fits with h fixed at 1..4 and with h free, each from 400 (600) random Nelder-Mead starts;
   SSE, RMSE, AIC, BIC with the error variance counted as a parameter (Table 'tab:ic').
2. Shapiro-Wilk on the free-h residuals.
3. Profile likelihood of h on a grid of step 0.02; 95% interval at chi2(1,0.95) = 3.841.
4. Standardised residuals r/sigma_hat of the h = 1, 2, 3 fits (sigma_hat from the h = 2 fit,
   nine degrees of freedom): maximum and sign runs.
5. The floorless model U = U0 K/(K + g^h): 2*dL, SSE ratio, fitted h.
6. The refractory-fraction model U = 3q + 3(1-q) K/(K + g^h): q, h, 2*dL, AIC comparison.
7. End-point drops: h with the highest, the lowest, and the last two points removed.
8. The local logit slope eta = dlog f/dlog l implied by the h = 2 fit (eq. hillid), at
   theta = 1/2 and at l = S, and the two plateaus theta_0, theta_inf.
9. U_min with its standard error from the covariance of the free fit.
Writes results/fit_titration.json and results/profile_h.npy (used by fig5_integer.py).
"""
import _path  # noqa: F401
import os, numpy as np
from scipy.optimize import minimize, curve_fit, brentq
from scipy.stats import shapiro, chi2
from esbm import RESULTS, save_json
from esbm.data import digitised_titration
from esbm.model import hill4
rng = np.random.default_rng(2026)
g, U = digitised_titration(); n = len(U)

def pred(th, h, x=g):
    U0, Um, lK = th; K = np.exp(lK)
    return Um + (U0-Um)*K/(K+x**h)
def sse(th, h): r = U-pred(th, h); return r@r
def fit(h, starts=400):
    best = None
    for _ in range(starts):
        x0 = np.array([rng.uniform(2.5, 3.5), rng.uniform(0.0, 0.9), rng.uniform(-6, 4)])
        r = minimize(sse, x0, args=(h,), method='Nelder-Mead', options=dict(maxiter=20000, xatol=1e-12, fatol=1e-14))
        if best is None or r.fun < best.fun: best = r
    return best
def nllh(S): return 0.5*n*np.log(2*np.pi*S/n)+0.5*n
def ic(S, k): L = nllh(S); return 2*k+2*L, k*np.log(n)+2*L

out = {}
print("1) model selection (k counts sigma)")
rows = []; fits = {}
for h in (1, 2, 3, 4):
    b = fit(h); fits[h] = b; A_, B_ = ic(b.fun, 4); rows.append(("h=%d" % h, b.fun, A_, B_, 4))
best = None
for _ in range(600):
    x0 = np.array([rng.uniform(2.5, 3.5), rng.uniform(0, 0.9), rng.uniform(-6, 4), rng.uniform(0.5, 4)])
    r = minimize(lambda p: sse(p[:3], p[3]), x0, method='Nelder-Mead', options=dict(maxiter=40000, xatol=1e-12, fatol=1e-14))
    if best is None or r.fun < best.fun: best = r
Sf = best.fun; Af, Bf = ic(Sf, 5); rows.append(("h free", Sf, Af, Bf, 5))
amin = min(r[2] for r in rows); bmin = min(r[3] for r in rows)
print("%-8s%11s%9s%4s%10s%10s%9s" % ("model", "SSE", "RMSE", "k", "AIC", "BIC", "dAIC"))
tab = {}
for lbl, S, A_, B_, k in rows:
    print("%-8s%11.4f%9.4f%4d%10.2f%10.2f%9.2f" % (lbl, S, np.sqrt(S/n), k, A_, B_, A_-amin))
    tab[lbl] = dict(SSE=S, RMSE=np.sqrt(S/n), k=k, AIC=A_, BIC=B_, dAIC=A_-amin, dBIC=B_-bmin)
hfree = best.x[3]; S50 = np.exp(best.x[2])**(1/hfree)
print("   free-h optimum: h = %.4f, S0.5 = %.4f mM" % (hfree, S50))
r = U-pred(best.x[:3], hfree); sw = shapiro(r)
print("   Shapiro-Wilk W = %.4f, p = %.3f" % (sw.statistic, sw.pvalue))
out['table'] = tab; out['h_free'] = hfree; out['S50_free_mM'] = S50; out['shapiro_W'] = sw.statistic; out['shapiro_p'] = sw.pvalue
out['dAIC_h1_vs_h2'] = tab['h=1']['AIC']-tab['h=2']['AIC']; out['dAIC_h3_vs_h2'] = tab['h=3']['AIC']-tab['h=2']['AIC']
out['dAIC_h2_vs_free'] = tab['h free']['AIC']-tab['h=2']['AIC']

print("\n2) profile likelihood of h")
Lmin = nllh(Sf); grid = np.arange(0.80, 4.01, 0.02); prof = []
for h in grid:
    b = fit(h, starts=120); prof.append(2*(nllh(b.fun)-Lmin))
prof = np.array(prof); thr = chi2.ppf(0.95, 1); inside = grid[prof <= thr]
np.save(os.path.join(RESULTS, "profile_h.npy"), np.vstack([grid, prof]))
print("   minimum at h = %.2f ; 95%% profile interval [%.2f, %.2f]" % (grid[np.argmin(prof)], inside.min(), inside.max()))
print("   2dL at h=1: %.1f ; h=2: %.2f ; h=3: %.1f" % (np.interp(1, grid, prof), np.interp(2, grid, prof), np.interp(3, grid, prof)))
out['profile'] = dict(h_min=grid[np.argmin(prof)], lo=inside.min(), hi=inside.max(),
                      stat_h1=np.interp(1, grid, prof), stat_h2=np.interp(2, grid, prof), stat_h3=np.interp(3, grid, prof))

print("\n3) standardised residuals (sigma_hat of the h=2 fit, %d dof)" % (n-3))
s2 = np.sqrt(fits[2].fun/(n-3)); out['sigma_hat_h2'] = s2
print("   sigma_hat = %.4f" % s2)
resid = {}
for h in (1, 2, 3):
    rr = (U-pred(fits[h].x, h))/s2
    runs = 1+int(np.sum(np.sign(rr[1:]) != np.sign(rr[:-1])))
    resid[h] = dict(max_abs=float(np.abs(rr).max()), signs=''.join('+' if v > 0 else '-' for v in rr), runs=runs)
    print("   h=%d  max|r/sigma| = %.2f  signs %s  (%d runs)" % (h, np.abs(rr).max(), resid[h]['signs'], runs))
out['standardised_residuals'] = resid

print("\n4) floorless model U = U0 K/(K+g^h)")
def hill0(x, U0, S, h): return U0/(1+(x/S)**h)
q0, _ = curve_fit(hill0, g, U, p0=[3, 0.55, 2], maxfev=400000)
S0 = ((U-hill0(g, *q0))**2).sum()
out['floorless'] = dict(h=q0[2], SSE=S0, SSE_ratio=S0/Sf, two_dL=n*np.log(S0/Sf))
print("   h = %.2f ; SSE ratio %.1f ; 2dL = %.1f (threshold 3.84)" % (q0[2], S0/Sf, n*np.log(S0/Sf)))

print("\n5) refractory-fraction model U = 3q + 3(1-q) K/(K+g^h)")
def hillq(x, q, S, h): return 3*q + 3*(1-q)/(1+(x/S)**h)
qq, _ = curve_fit(hillq, g, U, p0=[0.15, 0.55, 2], maxfev=400000)
Sq = ((U-hillq(g, *qq))**2).sum()
out['refractory'] = dict(q=qq[0], h=qq[2], SSE=Sq, two_dL=n*np.log(Sq/Sf), dAIC_vs_mechanistic=ic(Sq, 4)[0]-ic(Sf, 5)[0])
print("   q = %.3f, h = %.2f ; 2dL = %.2f against the free fit ; dAIC (artefact - mechanistic) = %.2f" % (qq[0], qq[2], n*np.log(Sq/Sf), ic(Sq, 4)[0]-ic(Sf, 5)[0]))

print("\n6) end-point drops")
def fit_h(gg, UU):
    f = lambda p: np.sum((p[1]+(p[0]-p[1])*np.exp(p[2])/(np.exp(p[2])+gg**p[3])-UU)**2)
    b = None
    for s0 in ([2.97, .45, -1.2, 2.], [3., .4, -1., 1.5], [2.9, .5, -1.5, 2.5]):
        r = minimize(f, s0, method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-14, maxiter=20000, maxfev=20000))
        if b is None or r.fun < b.fun: b = r
    return b.x[3]
ends = dict(drop_highest=fit_h(g[:-1], U[:-1]), drop_lowest=fit_h(g[1:], U[1:]), drop_last_two=fit_h(g[:-2], U[:-2]))
print("   drop highest %.3f ; drop lowest %.3f ; drop last two %.3f" % (ends['drop_highest'], ends['drop_lowest'], ends['drop_last_two']))
out['endpoints'] = ends

print("\n7) eq. hillid: plateaus and local logit slope of the h=2 fit")
p4, cov4 = curve_fit(hill4, g, U, p0=[3, 0.45, 0.55, 2], maxfev=400000); se4 = np.sqrt(np.diag(cov4))
out['free_fit'] = dict(U0=p4[0], Umin=p4[1], S=p4[2], h=p4[3], se=list(se4))
print("   free fit: U0 %.3f+-%.3f  Umin %.3f+-%.3f  S %.3f+-%.3f  h %.3f+-%.3f" % tuple(np.ravel(list(zip(p4, se4)))))
U0, Um, S = fits[2].x[0], fits[2].x[1], np.exp(fits[2].x[2])**0.5   # h=2 fit: K = S^2
th0, thi = U0/3, Um/3
w0 = 1/th0-1; fmax = 1/thi-1                    # in units of Khat
KL = S*np.sqrt((1+fmax)/(1+w0))
def eta_at(l):
    x = (l/KL)**2; f = w0+(fmax-w0)*x/(1+x)
    return (fmax-w0)*2*x/(1+x)**2/f
l_half = brentq(lambda l: w0+(fmax-w0)*(l/KL)**2/(1+(l/KL)**2)-1.0, 1e-6, 1e6)   # f = Khat  <=>  theta = 1/2
out['hillid'] = dict(theta0=th0, theta_inf=thi, w0_over_Khat=w0, fmax_over_Khat=fmax, swing=fmax/w0,
                     eta_at_theta_half=eta_at(l_half), eta_at_S=eta_at(S), S_h2=S, K_L=KL)
print("   theta0 %.3f theta_inf %.3f  w0/Khat %.4f  fmax/Khat %.2f  swing %.0f" % (th0, thi, w0, fmax, fmax/w0))
print("   eta at theta=1/2: %.3f ; eta at l=S: %.3f" % (eta_at(l_half), eta_at(S)))
save_json("fit_titration", out)
