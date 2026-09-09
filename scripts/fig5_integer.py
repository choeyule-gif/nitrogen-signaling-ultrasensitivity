"""Fig 5: cooperative effector binding sets the sensor factor.  Run fit_titration.py and mcmc_h.py first.

(a) theta/(1-theta) against l for h = 1, 2, 3: slope -h;
(b) the digitised titration with fits at fixed h = 1, 2, 3 and the +-sigma_hat band of the h = 2 fit;
(c) standardised residuals r/sigma_hat of the three fits, band +-2 sigma_hat;
(d) likelihood ratio L/L_max for h (profile) with the marginal posterior from MCMC; dashed line,
    the chi2(1, 0.95) threshold; bars, the 95% profile interval and credible interval.
"""
import _path  # noqa: F401
import os, numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import gaussian_kde
from esbm import FIGURES, RESULTS, load_json, save_csv
from esbm.data import digitised_titration
from esbm.style import MC, GRAY, use, tidy, panel
use()
rng = np.random.default_rng(4)
g, U = digitised_titration(); n = len(U)
def pred(t, h, x): U0, Um, lK = t; K = np.exp(lK); return Um+(U0-Um)*K/(K+x**h)
def sse(t, h): r = U-pred(t, h, g); return r@r
def fit(h, st=300):
    best = None
    for _ in range(st):
        x0 = np.array([rng.uniform(2.5, 3.5), rng.uniform(0, .9), rng.uniform(-6, 4)])
        r = minimize(sse, x0, args=(h,), method='Nelder-Mead', options=dict(maxiter=20000, xatol=1e-12, fatol=1e-14))
        if best is None or r.fun < best.fun: best = r
    return best
ft = load_json("fit_titration"); mc = load_json("mcmc_h")
pg, pv = np.load(os.path.join(RESULTS, "profile_h.npy")); hh = np.load(os.path.join(RESULTS, "mcmc_h.npy"))

fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.4)); ax = ax.ravel()
l = np.logspace(-2, 2, 600); Kh = 1.0
for i, h in enumerate((1, 2, 3)):
    th = Kh/(Kh+l**h); ax[0].plot(l, th/(1-th), color=MC[i], label=r'$h=%d$' % h)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlim(1e-2, 1e2); ax[0].set_ylim(1e-6, 1e6)
ax[0].set_xlabel(r'$l$'); ax[0].set_ylabel(r'$\theta/(1-\theta)$')
ax[0].legend(loc='upper right'); tidy(ax[0], 1, 1, nx=3, ny=4); panel(ax[0], 'a')

fits = {h: fit(h) for h in (1, 2, 3)}
xx = np.logspace(np.log10(0.014), np.log10(14), 500)
s2 = np.sqrt(fits[2].fun/(n-3))
ax[1].fill_between(xx, pred(fits[2].x, 2, xx)-s2, pred(fits[2].x, 2, xx)+s2, color=GRAY, alpha=0.55, lw=0)
for i, h in enumerate((1, 2, 3)):
    ax[1].plot(xx, pred(fits[h].x, h, xx), color=MC[i], label=r'$h=%d$' % h)
ax[1].plot(g, U, 'o', mfc='white', mec='k', mew=1.2, ms=5.2, zorder=5)
ax[1].set_xscale('log'); ax[1].set_xlim(0.014, 14); ax[1].set_ylim(0, 3.3)
ax[1].set_xlabel(r'$[\mathrm{Gln}]\ /\ \mathrm{mM}$'); ax[1].set_ylabel(r'$\overline{n}_{\mathrm{UMP}}$')
ax[1].legend(loc='lower left'); tidy(ax[1], 1, 0, nx=3, ny=4); panel(ax[1], 'b')

for i, h in enumerate((1, 2, 3)):
    r = (U-pred(fits[h].x, h, g))/s2
    ax[2].plot(g, r, 'o-', color=MC[i], ms=4.5, lw=1.0, mfc='white', mew=1.2)
ax[2].axhline(0, color='k', lw=1.0)
ax[2].axhspan(-2, 2, color=GRAY, alpha=0.35, lw=0, zorder=0)
ax[2].set_xscale('log'); ax[2].set_xlim(0.014, 14); ax[2].set_ylim(-8, 8); ax[2].set_yticks([-8, -4, 0, 4, 8])
ax[2].set_xlabel(r'$[\mathrm{Gln}]\ /\ \mathrm{mM}$'); ax[2].set_ylabel(r'$r/\hat{\sigma}$')
tidy(ax[2], 1, 0, nx=3, ny=4); ax[2].set_yticks([-8, -4, 0, 4, 8]); panel(ax[2], 'c')

LR = np.exp(-pv/2)
kd = gaussian_kde(hh, bw_method=0.18); xs = np.linspace(1.55, 2.62, 700); dens = kd(xs); dens /= dens.max()
ax[3].fill_between(xs, 0, dens, color=GRAY, alpha=0.85, lw=0, label=r'MCMC')
ax[3].plot(pg, LR, '-', color='k', lw=1.7, label=r'$L/L_{\max}$')
ax[3].axhline(np.exp(-3.841/2), color='k', ls='--', lw=1.0)
ax[3].plot([ft['profile']['lo'], ft['profile']['hi']], [0.055, 0.055], '-', color='k', lw=3.2, solid_capstyle='butt')
ax[3].plot(mc['ci95'], [0.020, 0.020], '-', color='0.55', lw=3.2, solid_capstyle='butt')
ax[3].set_xlim(1.55, 2.62); ax[3].set_ylim(0, 1.06)
ax[3].set_xlabel(r'$h$'); ax[3].set_ylabel(r'$L/L_{\max}$')
ax[3].legend(loc='upper left'); tidy(ax[3], 0, 0)
ax[3].set_xticks([1.6, 2.0, 2.4]); ax[3].set_yticks([0, 0.5, 1.0]); panel(ax[3], 'd')
plt.tight_layout(w_pad=1.7, h_pad=1.3)
plt.savefig(os.path.join(FIGURES, 'Fig5.pdf'), metadata={"CreationDate": None})
save_csv("Fig5b_fits", dict(glutamine_mM=xx, **{"fit_h%d" % h: pred(fits[h].x, h, xx) for h in (1, 2, 3)}), "Fig 5b; the data points are data/jiangninfa2011_fig2B_digitised.csv")
save_csv("Fig5c_standardised_residuals", dict(glutamine_mM=g, **{"r_over_sigma_h%d" % h: (U-pred(fits[h].x, h, g))/s2 for h in (1, 2, 3)}), "Fig 5c, sigma_hat = %.4f" % s2)
save_csv("Fig5d_profile_and_posterior", dict(h=pg, likelihood_ratio=LR, **{"posterior_density_at_h": np.interp(pg, xs, dens)}), "Fig 5d")
print("Fig5 written; max |r/sigma|:", {h: round(float(np.abs((U-pred(fits[h].x, h, g))/s2).max()), 2) for h in (1, 2, 3)})
