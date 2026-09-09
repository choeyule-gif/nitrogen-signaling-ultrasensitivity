"""Fig 2: identical, independent sites do not sharpen the response.

(a) theta(l) for n = 1, 2, 3, 12 -- the four curves coincide;
(b) the local Hill coefficient, equal to 1, computed by numerical differentiation of the
    full steady state rather than from the closed form;
(c) Hill coefficient of the TOTAL modified fraction against T_E/T_S, at n = 3,
    A*Khat = 1, Bhat = Ghat = 0.4/Khat: sequestration only lowers it.
Reports: max |n_H - 1| in (b); n_H at T_E/T_S = 0.1 and 0.2 in (c).
"""
import _path  # noqa: F401  (makes esbm importable when run as a script)
import os, numpy as np
from math import comb
import matplotlib.pyplot as plt
from esbm import FIGURES, save_json, save_csv
from esbm.model import nH_local_tot
from esbm.style import MC, use_compact, no_minor
use_compact(**{'legend.fontsize': 9})

Kh = 1.0
def theta_exact(l, n):
    u = Kh/l
    num = sum(i*comb(n, i)*u**i for i in range(n+1)); den = sum(comb(n, i)*u**i for i in range(n+1))
    return num/(n*den)
l = np.logspace(-2.5, 2.5, 400); NS = (1, 2, 3, 12)
def nH(l, n):
    d = 1e-6
    g = lambda x: np.log(theta_exact(x, n)/(1-theta_exact(x, n)))
    return -(g(l*(1+d))-g(l*(1-d)))/(2*d)
lb = np.logspace(-2.5, 2.5, 200)
dev = max(abs(nH(x, n)-1) for n in NS for x in lb)
print("panel (b): max |n_H - 1| over four n and five decades = %.2e" % dev)
TS_c = 1.0
xc = np.logspace(-3, 1.4, 200)
yc = np.array([nH_local_tot(t*TS_c, TS_c) for t in xc])
nh01, nh02 = nH_local_tot(0.1, TS_c), nH_local_tot(0.2, TS_c)
print("panel (c): nH at T_E/T_S = 0.1 is %.5f, at 0.2 is %.5f ; range over the sweep [%.4f, %.4f]" % (nh01, nh02, yc.min(), yc.max()))

fig = plt.figure(figsize=(7.373, 2.370))
W, H = 530.857, 170.623
def ax_at(x0, x1): return fig.add_axes([x0/W, (H-141.45)/H, (x1-x0)/W, (141.45-7.47)/H])
ax = ax_at(35.33, 160.87)
for i, n in enumerate(NS):
    ax.semilogx(l, [theta_exact(x, n) for x in l], color=MC[i], label=r'$n=%d$' % n)
ax.set_xlim(10**-2.5, 10**2.5); ax.set_ylim(-0.02, 1.02)
ax.set_xticks([1e-2, 1e0, 1e2]); ax.set_yticks([0.0, 0.3, 0.6, 0.9])
ax.set_xlabel(r'$l$'); ax.set_ylabel(r'$\theta$')
ax.legend(loc='upper right', handlelength=1.4, borderpad=0.35, labelspacing=0.35)
ax = ax_at(219.61, 345.14)
for i, n in enumerate(NS):
    ax.semilogx(lb, [nH(x, n) for x in lb], color=MC[i])
ax.set_xlim(10**-2.5, 10**2.5); ax.set_ylim(0.978, 1.022)
ax.set_xticks([1e-2, 1e0, 1e2]); ax.set_yticks([0.98, 1.00, 1.02])
ax.set_xlabel(r'$l$'); ax.set_ylabel(r'$n_H$')
ax = ax_at(403.88, 529.42)
ax.axhline(1.0, color='k', ls=(0, (3.7, 1.6)), lw=0.8)
ax.semilogx(xc, yc, color=MC[0])
ax.set_xlim(1e-3, 10**1.4); ax.set_ylim(0.62, 1.06)
ax.set_xticks([1e-3, 1e-1, 1e1]); ax.set_yticks([0.75, 0.90, 1.05])
ax.set_xlabel(r'$T_E/T_S$'); ax.set_ylabel(r'$n_H$')
no_minor(fig)
for x, s in ((1.4, 'a'), (185.7, 'b'), (370.0, 'c')):
    fig.text(x/W, 1-9.7/H, s, fontsize=11, fontweight='bold', va='baseline')
fig.savefig(os.path.join(FIGURES, "Fig2.pdf"))
save_csv("Fig2a_theta", dict(l=l, **{"theta_n%d" % n: [theta_exact(x, n) for x in l] for n in NS}), "Fig 2a")
save_csv("Fig2b_local_hill", dict(l=lb, **{"nH_n%d" % n: [nH(x, n) for x in lb] for n in NS}), "Fig 2b")
save_csv("Fig2c_total_fraction_hill", dict(TE_over_TS=xc, nH=yc), "Fig 2c: n=3, A*Khat=1, Bhat=Ghat=0.4/Khat, T_S=Khat")
save_json("fig2_sites", dict(max_abs_nH_minus_1=dev, nH_TE_over_TS_0p1=nh01, nH_TE_over_TS_0p2=nh02,
                             nH_sweep_min=yc.min(), nH_sweep_max=yc.max()))
