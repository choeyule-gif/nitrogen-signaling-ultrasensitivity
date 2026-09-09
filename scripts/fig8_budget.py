"""Fig 8: what the two layers account for.  Run scripts/cascade_budget.py first.

(a) the budget at each PII concentration: first layer, composite with the readout, and
    the measured cascade coefficient;
(b) the residual, with bootstrap intervals, against PII -- and the band a mechanism
    monotone in PII could occupy;
(c) the readout factor against the second layer's half-point, over twenty-four decades
    and both readings of l_2.
"""
import _path  # noqa: F401
import os, numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from esbm import FIGURES, load_json, save_csv
from esbm.data import constants
from esbm.model import n1090
from esbm.style import MC, use_compact, no_minor
use_compact()

cb = load_json("cascade_budget"); cc = constants()["jiangninfa2011"]
PII = np.array(cc["cascade_PII_uM"], float); HALF = cc["cascade_halfpoint_mM"]; MEAS = np.array(cc["cascade_hill_90"])
L1 = cb['lin']['first_layer']; COMP = np.array(cb['lin']['composite'])
RES = np.array([cb['bootstrap']['PII_%g' % p]['residual'][0] for p in cc["cascade_PII_uM"]])
RLO = np.array([cb['bootstrap']['PII_%g' % p]['residual'][1] for p in cc["cascade_PII_uM"]])
RHI = np.array([cb['bootstrap']['PII_%g' % p]['residual'][2] for p in cc["cascade_PII_uM"]])
th0, thi, S, h = cb['fit']['theta0'], cb['fit']['theta_inf'], cb['fit']['S_mM'], cb['fit']['h']
theta = lambda g: thi + (th0-thi)/(1+(g/S)**h)
def fchain(g, r):
    t = theta(g); l1 = (1-t)**3
    return ((1-l1) if r == 'compl' else t)/l1
n1 = n1090(lambda g: (th0-theta(g))/(th0-thi))
def factor(K2, r):
    f0, fi = fchain(1e-14, r), fchain(1e14, r)
    z = lambda g: ((f0-fchain(g, r))/(f0-fi))*((K2+fi)/(K2+fchain(g, r)))
    return n1090(z)/n1
K2g = np.logspace(-12, 12, 420)
fac_lin = np.array([factor(k, 'lin') for k in K2g]); fac_cmp = np.array([factor(k, 'compl') for k in K2g])
K2star = cb['lin']['K2_pinned']

fig = plt.figure(figsize=(7.15, 2.5))
def ax_at(l, w): return fig.add_axes([l, 0.20, w, 0.755])
ax = ax_at(0.062, 0.235)
x = np.arange(3); w = 0.26
ax.bar(x-w, [L1]*3, w, color=MC[0], ec='k', lw=0.6)
ax.bar(x, COMP, w, color=MC[2], ec='k', lw=0.6)
ax.bar(x+w, MEAS, w, color='w', ec='k', lw=0.9, hatch='///')
ax.set_xticks(x); ax.set_xticklabels(['36', '5', '0.5'])
ax.set_xlim(-0.55, 2.55); ax.set_ylim(0, 7.2); ax.set_yticks([0, 2, 4, 6])
ax.set_xlabel(r'$T_{\mathrm{PII}}\ /\ \mu$M'); ax.set_ylabel(r'$n_H$')
ax = ax_at(0.395, 0.235)
lo, hi = min(RES[0], RES[2]), max(RES[0], RES[2])
ax.fill_between([0.3, 60], lo, hi, color='0.85', lw=0, zorder=0)
ax.axhline(1.0, color='k', ls=(0, (3.7, 1.6)), lw=0.8, zorder=1)
ax.errorbar(PII, RES, yerr=[RES-RLO, RHI-RES], fmt='o', color=MC[1], ms=5.4, capsize=2.8, lw=1.3, zorder=3)
ax.semilogx(PII, RES, color=MC[1], lw=1.0, zorder=2)
ax.set_xlim(60, 0.3); ax.set_ylim(0.6, 2.0); ax.set_xticks([10, 1]); ax.set_yticks([0.8, 1.2, 1.6, 2.0])
ax.set_xlabel(r'$T_{\mathrm{PII}}\ /\ \mu$M'); ax.set_ylabel('residual')
ax = ax_at(0.728, 0.245)
ax.axhline(3.0, color='k', ls=(0, (1.5, 1.5)), lw=0.9)
ax.semilogx(K2g, fac_lin, color=MC[0], label=r'$l_2\propto\theta$')
ax.semilogx(K2g, fac_cmp, color=MC[1], label=r'$l_2=1-(1-\theta)^3$')
for k in K2star: ax.axvline(k, color='0.55', lw=0.7, ls=(0, (2.5, 2.0)))
ax.set_xlim(1e-8, 1e8); ax.set_ylim(0.9, 3.25); ax.set_xticks([1e-6, 1e-2, 1e2, 1e6]); ax.set_yticks([1, 2, 3])
ax.set_xlabel(r'$\hat K_2$'); ax.set_ylabel('readout factor')
ax.legend(loc='upper left', handlelength=1.3, borderpad=0.3, labelspacing=0.25)
no_minor(fig)
for xx, s in ((0.006, 'a'), (0.339, 'b'), (0.672, 'c')):
    fig.text(xx, 0.945, s, fontsize=11, fontweight='bold', va='baseline')
fig.savefig(os.path.join(FIGURES, "Fig8.pdf"), metadata={"CreationDate": None})
save_csv("Fig8ab_budget", dict(PII_uM=PII, first_layer=[L1]*3, composite=COMP, measured=MEAS, residual=RES, residual_lo=RLO, residual_hi=RHI), "Fig 8a,b")
save_csv("Fig8c_readout_factor", dict(Khat2=K2g, factor_linear=fac_lin, factor_complementary=fac_cmp), "Fig 8c")
print("Fig8 written; readout sup linear %.5f complementary %.5f" % (fac_lin.max(), fac_cmp.max()))
