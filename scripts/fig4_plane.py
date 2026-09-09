"""Fig 4: the (nu, eta) plane.

Grey: iso-n_H hyperbolae nu*eta = const.  Blue: the measured n_H = 2.02 of the first
cycle.  Orange: the ceiling eq:central puts on eta, with the Mg2+ and Mn2+ affinity
ratios marked.  Purple: the ceiling eq:nubound puts on nu at n = 3.  Black: the
measured (nu, eta) = (1.11 +- 0.11, 1.82 +- 0.18).  Right-hand ticks: the integers eta = h.
Reports the two eta ceilings and the product nu*eta.
"""
import _path  # noqa: F401  (makes esbm importable when run as a script)
import os, numpy as np
import matplotlib.pyplot as plt
from esbm import FIGURES, save_json, save_csv
from esbm.data import constants
from esbm.style import MC, use_compact, no_minor
use_compact()
c = constants()["jiang1998a"]
eta_of = lambda r: 2*np.sqrt(r)/(1+np.sqrt(r))
r_mg = c["Kact_glutamine_UR_uM_Mg"]/c["Ki_glutamine_UT_uM_Mg"][0]     # 80/70
r_mn = c["Ki_glutamine_UT_uM_Mn"]/c["Kact_glutamine_UR_uM_Mn"]         # 700/150
eta_mg, eta_mn = eta_of(r_mg), eta_of(r_mn)
NU, NU_E, ETA, ETA_E, NH = 1.11, 0.11, 1.82, 0.18, 2.02

fig = plt.figure(figsize=(3.55, 3.2))
ax = fig.add_axes([0.16, 0.145, 0.79, 0.835])
x0, x1, y0, y1 = 0.6, 3.2, 0.6, 3.2
xx = np.linspace(x0, x1, 500)
for cc in (1, 1.5, 3, 4, 6):
    yy = cc/xx; m = (yy >= y0) & (yy <= y1)
    if m.sum() < 3: continue
    ax.plot(xx[m], yy[m], color='0.78', lw=0.7, zorder=1)
    xm, ym = xx[m][len(xx[m])//2], yy[m][len(yy[m])//2]
    ax.text(xm, ym, f'{cc:g}', color='0.45', fontsize=7, ha='center', va='center', zorder=1, bbox=dict(fc='w', ec='none', pad=0.6))
yy = NH/xx; m = (yy >= y0) & (yy <= y1)
ax.plot(xx[m], yy[m], color=MC[0], lw=1.6, zorder=3)
ax.text(2.45, NH/2.45+0.055, f'{NH:g}', color=MC[0], fontsize=8, ha='center', va='bottom', zorder=3, bbox=dict(fc='w', ec='none', pad=0.6))
ax.fill_between([x0, x1], 2.0, y1, color=MC[1], alpha=0.10, lw=0, zorder=0)
ax.axhline(2.0, color=MC[1], lw=1.4, zorder=2)
for e, mk in ((eta_mg, 'o'), (eta_mn, 's')):
    ax.plot([x0, x1], [e, e], color=MC[1], lw=0.8, ls=(0, (1, 1.8)), zorder=2)
    ax.plot([x0+0.09], [e], mk, ms=4.8, mec=MC[1], mfc='w', mew=1.3, zorder=4)
ax.fill_betweenx([y0, y1], 3.0, x1, color=MC[3], alpha=0.10, lw=0, zorder=0)
ax.axvline(3.0, color=MC[3], lw=1.4, zorder=2)
for h in (1, 2, 3):
    ax.plot([x1, x1-0.075], [h, h], color=MC[0], lw=1.6, solid_capstyle='butt', zorder=4)
ax.errorbar([NU], [ETA], xerr=[NU_E], yerr=[ETA_E], fmt='o', color='k', ms=5.4, capsize=2.8, lw=1.3, zorder=5)
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_xticks([1, 2, 3]); ax.set_yticks([1, 2, 3])
ax.set_xlabel(r'$\nu$'); ax.set_ylabel(r'$\eta$')
no_minor(fig)
fig.savefig(os.path.join(FIGURES, "Fig4.pdf"))
print("eta ceilings: Mg2+ %.4f (ratio %.3f)  Mn2+ %.4f (ratio %.3f);  point (%.2f,%.2f) product %.3f" % (eta_mg, r_mg, eta_mn, r_mn, NU, ETA, NU*ETA))
save_csv("Fig4_plane_markers", dict(quantity=["nu", "nu_err", "eta", "eta_err", "nH_first_cycle", "eta_ceiling_Mg", "eta_ceiling_Mn", "nu_ceiling_n3"], value=[NU, NU_E, ETA, ETA_E, NH, eta_mg, eta_mn, 3.0]), "Fig 4")
save_json("fig4_plane", dict(eta_ceiling_Mg=eta_mg, eta_ceiling_Mn=eta_mn, ratio_Mg=r_mg, ratio_Mn=r_mn, nu_eta_product=NU*ETA))
