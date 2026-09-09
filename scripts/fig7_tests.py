"""Fig 7: the two measurements the analysis asks for.

(a) equilibrium binding of glutamine to GlnD, drawn for one site and for two
    cooperative sites at the same half-point;
(b)-(d) the glutamine titration predicted at two enzyme and two target levels under
    three accounts, all solved in the same mass-action framework:
    (b) eqs (1)-(3) as they stand -- the four curves coincide;
    (c) with a catalytically active complex bridging two targets -- the target levels
        separate, the enzyme levels do not;
    (d) with the enzyme raised until it sequesters the effector -- the enzyme levels separate.
Reports the local coefficient at the half-point in each condition (caption: 1.000; 1.08 vs 1.26).
"""
import _path  # noqa: F401
import os, numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from esbm import FIGURES, save_json, save_csv
from esbm.style import MC, use_compact, no_minor
use_compact(**{'lines.linewidth': 1.5, 'legend.fontsize': 7.5})

A, Bh, Gh, Kh = 1.0, 0.4, 0.4, 1.0
C = Bh + A*Gh*Kh
LEV = [(1.0, 1.0), (2.0, 1.0), (1.0, 4.0), (2.0, 4.0)]      # (enzyme x, target x)
LAB = [r'$1\times,1\times$', r'$2\times,1\times$', r'$1\times,4\times$', r'$2\times,4\times$']

def theta_bridge(l, TS, c):
    """n = 1 with a catalytically active complex bridging two targets: s_1 = s_0 Khat/(l - c s_0)."""
    hi = min(TS, l/c)*(1-1e-12)
    s0 = brentq(lambda s0: s0 + s0*Kh/(l - c*s0) - TS, 1e-14, hi, xtol=1e-18, rtol=8.9e-16)
    s1 = s0*Kh/(l - c*s0)
    return s1/(s0+s1)

def theta_tot(l, TE, TS):
    n = 3; u = Kh/l; P = (1+u)**n; Pm = (1+u)**(n-1); Pm2 = (1+u)**(n-2)
    f = lambda e: e*(1+A*l) + n*e*(TS/(P+n*e*C*Pm))*C*Pm - TE
    hi = 1.0
    while f(hi) < 0: hi *= 10
    e = brentq(f, 1e-300, hi, xtol=1e-18, rtol=8.9e-16)
    s0 = TS/(P + n*e*C*Pm)
    M = s0*(n*u*Pm + Bh*e*n*(n-1)*u*Pm2 + A*Gh*Kh*e*n*Pm2*(1+n*u))
    TL = l + A*e*l + A*Gh*Kh*n*e*s0*Pm
    return M/(n*TS), TL

ll = np.logspace(-1.2, 1.2, 260)
fig = plt.figure(figsize=(7.15, 2.05))
def ax_at(l, w): return fig.add_axes([l, 0.245, w, 0.705])
ax = ax_at(0.072, 0.172)
for hh, col in ((1, MC[0]), (2, MC[1])):
    ax.semilogx(ll, ll**hh/(1.0 + ll**hh), color=col, label=r'$h=%d$' % hh)
ax.set_xlim(ll[0], ll[-1]); ax.set_ylim(-0.02, 1.02); ax.set_xticks([0.1, 1, 10]); ax.set_yticks([0, 0.5, 1])
ax.set_xlabel(r'$[\mathrm{Gln}]$'); ax.set_ylabel('fraction bound')
ax.legend(loc='upper left', handlelength=1.2, borderpad=0.3, labelspacing=0.22)
ax = ax_at(0.328, 0.172)
for (fe, ft), col, lab in zip(LEV, MC, LAB):
    ax.semilogx(ll, [theta_tot(x, 1e-4*fe, 1.0*ft)[0] for x in ll], color=col, label=lab)
ax.set_xlim(ll[0], ll[-1]); ax.set_ylim(-0.02, 1.02); ax.set_xticks([0.1, 1, 10]); ax.set_yticks([0, 0.5, 1])
ax.set_xlabel(r'$l$'); ax.set_ylabel(r'$\theta$')
ax.legend(loc='upper right', handlelength=1.1, borderpad=0.25, labelspacing=0.18)
ax = ax_at(0.584, 0.172)
for (fe, ft), col in zip(LEV, MC):
    ax.semilogx(ll, [theta_bridge(x, 1.0*ft, c=0.35) for x in ll], color=col)
ax.set_xlim(ll[0], ll[-1]); ax.set_ylim(-0.02, 1.02); ax.set_xticks([0.1, 1, 10]); ax.set_yticks([0, 0.5, 1])
ax.set_xlabel(r'$l$'); ax.set_ylabel(r'$\theta$')
ax = ax_at(0.840, 0.172)
for (fe, ft), col in zip(LEV, MC):
    out = [theta_tot(x, 0.6*fe, 1.0*ft) for x in ll]
    ax.semilogx([o[1] for o in out], [o[0] for o in out], color=col)
ax.set_xlim(ll[0], ll[-1]); ax.set_ylim(-0.02, 1.02); ax.set_xticks([0.1, 1, 10]); ax.set_yticks([0, 0.5, 1])
ax.set_xlabel(r'$T_L$'); ax.set_ylabel(r'$\theta$')
no_minor(fig)
for x, s in ((0.006, 'a'), (0.262, 'b'), (0.518, 'c'), (0.774, 'd')):
    fig.text(x, 0.94, s, fontsize=11, fontweight='bold', va='baseline')
fig.savefig(os.path.join(FIGURES, "Fig7.pdf"))

def nH(y, x):
    yy = np.clip(np.array(y), 1e-9, 1-1e-9); gg = np.log(yy/(1-yy))
    i = np.argmin(np.abs(yy-0.5))
    return abs(np.gradient(gg, np.log(x))[i])
out = {}
for lab, fn in (("b_ours", lambda fe, ft: [theta_tot(x, 1e-4*fe, ft)[0] for x in ll]),
                ("c_bridge", lambda fe, ft: [theta_bridge(x, ft, c=0.35) for x in ll]),
                ("d_sequestering", lambda fe, ft: [theta_tot(x, 0.6*fe, ft)[0] for x in ll])):
    out[lab] = [nH(fn(fe, ft), ll) for fe, ft in LEV]
    print("%-15s" % lab, "  ".join("%.3f" % v for v in out[lab]))
cols = dict(l=ll)
for (fe, ft), lab in zip(LEV, ("1x1x", "2x1x", "1x4x", "2x4x")):
    cols["b_theta_%s" % lab] = [theta_tot(x, 1e-4*fe, ft)[0] for x in ll]
    cols["c_theta_bridge_%s" % lab] = [theta_bridge(x, ft, c=0.35) for x in ll]
    o = [theta_tot(x, 0.6*fe, ft) for x in ll]; cols["d_theta_%s" % lab] = [v[0] for v in o]; cols["d_TL_%s" % lab] = [v[1] for v in o]
save_csv("Fig7_titrations", cols, "Fig 7b-d: n=3, A=1, Bhat=Ghat=0.4, Khat=1; (enzyme x, target x)")
save_json("fig7_tests", out)
