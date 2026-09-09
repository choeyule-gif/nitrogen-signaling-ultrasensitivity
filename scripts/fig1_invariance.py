"""Fig 1: two invariance checks, both solved in the raw species variables.

(a) the ratios s_i/s_{i-1} at n = 4 while the free enzyme is swept over 5,000-fold at
    fixed free effector;
(b) the same ratios at n = 3 while a dead-end complex Z = E.S_0 has its association
    constant swept over eleven decades.
Neither solve imposes the chain relation; it is what is being tested.
Reports: the enzyme range, the drift of the ratios, and |ratio - K_i/l| (Fig 1 caption).
"""
import _path  # noqa: F401  (makes esbm importable when run as a script)
import os, numpy as np
from scipy.optimize import fsolve, brentq
import matplotlib.pyplot as plt
from esbm import FIGURES, save_json, save_csv
from esbm.style import MC, use_compact, no_minor
use_compact()

# ---------------- panel (a): non-identical ladder, n = 4, buffered effector
n_a, A_a, l_a = 4, 1.3, 0.9
Bh_a = np.array([0.80, 0.35, 1.20, 0.55])
Gh_a = np.array([0.60, 1.10, 0.45, 0.90])
rho_a = np.array([1.4, 0.7, 2.1, 1.0])
K_a = Bh_a*rho_a/(A_a*Gh_a)                       # K_{i+1}, i = 0..n-1
TS_a = 1.0
def ratios_a(TE, guess):
    """Raw steady state at fixed l and T_S: n cut balances plus the two totals."""
    def eqs(v):
        e = np.exp(v[0]); s = np.exp(v[1:])
        eps = A_a*e*l_a
        b = Bh_a*e*s[:-1]; g = Gh_a*eps*s[1:]
        out = list(b*rho_a - g)
        out.append(s.sum() + b.sum() + g.sum() - TS_a)
        out.append(e + eps + b.sum() + g.sum() - TE)
        return out
    v = fsolve(eqs, guess, full_output=False, xtol=1e-14)
    e = np.exp(v[0]); s = np.exp(v[1:])
    return e, s[1:]/s[:-1], v
TE_grid = np.logspace(-4, -0.3, 90)
E, Rr = [], []
guess = np.log([1e-4] + [TS_a/(n_a+1)/(1.6**i) for i in range(n_a+1)])
for TE in TE_grid:
    e, r, guess = ratios_a(TE, guess); E.append(e); Rr.append(r)
E = np.array(E); Rr = np.array(Rr)
o = np.argsort(E); E, Rr = E[o], Rr[o]
fold_a = E.max()/E.min(); drift_a = np.max(np.abs(Rr/Rr[0]-1)); dev_a = np.max(np.abs(Rr - K_a/l_a))
print("(a) free enzyme spans %.1f-fold; ratios drift %.2e ; |ratio - K_i/l| max %.2e" % (fold_a, drift_a, dev_a))

# ---------------- panel (b): dead-end complex, n = 3, identical independent sites
n_b, A_b, Bh_b, Gh_b, Kh_b, TS_b = 3, 1.0, 0.4, 0.4, 1.0, 1.0
C_b = Bh_b + A_b*Gh_b*Kh_b
TE_b, TL_b = 1e-9, 0.9 + 1e-9
def ratios_b(KZ):
    def e_of(l):
        u = Kh_b/l; P = (1+u)**n_b; Pm = (1+u)**(n_b-1); D = n_b*C_b*Pm + KZ
        f = lambda e: e*(1+A_b*l) + e*D*TS_b/(P+e*D) - TE_b
        hi = 1.0
        while f(hi) < 0: hi *= 10
        e = brentq(f, 1e-300, hi, xtol=1e-300, rtol=8.9e-16)
        return e, TS_b/(P+e*D), Pm
    def res(l):
        e, s0, Pm = e_of(l)
        return l + A_b*e*l + A_b*Gh_b*Kh_b*n_b*e*s0*Pm - TL_b
    l = brentq(res, 1e-6, 1e3, xtol=1e-300, rtol=8.9e-16)
    r = np.array([1.0]+[0.0]*n_b)
    for i in range(1, n_b+1): r[i] = r[i-1]*((n_b-i+1)/i*Kh_b)/l
    return r[1:]/r[:-1]
xb = np.logspace(-6, 5, 120)
yb = np.array([ratios_b(k) for k in xb])
drift_b = np.max(np.abs(yb/yb[0]-1))
print("(b) ratios drift %.2e over 11 decades at T_E/T_L = %.1e" % (drift_b, TE_b/TL_b))

# ---------------- draw
fig = plt.figure(figsize=(7.0, 2.45))
def ax_at(l, w): return fig.add_axes([l, 0.185, w, 0.775])
ax = ax_at(0.075, 0.375)
for i in range(n_a):
    ax.loglog(E, Rr[:, i], color=MC[i], label=r'$s_%d/s_%d$' % (i+1, i))
ax.set_xlim(E.min(), E.max()); ax.set_ylim(0.1, 10.0)
ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(['0.1', '1', '10'])
ax.set_xlabel(r'$e$'); ax.set_ylabel(r'$s_i/s_{i-1}$')
ax.legend(loc='lower center', ncol=2, handlelength=1.3, borderpad=0.3, labelspacing=0.25, columnspacing=1.0)
ax = ax_at(0.585, 0.375)
for i in range(n_b):
    ax.loglog(xb, yb[:, i], color=MC[i], label=r'$s_%d/s_%d$' % (i+1, i))
ax.set_xlim(1e-6, 1e5); ax.set_ylim(0.25, 5.0)
ax.set_xticks([1e-6, 1e-3, 1e0, 1e3]); ax.set_yticks([0.3, 1, 3]); ax.set_yticklabels(['0.3', '1', '3'])
ax.set_xlabel(r'$K_Z$'); ax.set_ylabel(r'$s_i/s_{i-1}$')
ax.legend(loc='center right', handlelength=1.3, borderpad=0.3, labelspacing=0.25)
no_minor(fig)
for x, s in ((0.008, 'a'), (0.518, 'b')):
    fig.text(x, 0.945, s, fontsize=11, fontweight='bold', va='baseline')
fig.savefig(os.path.join(FIGURES, "Fig1.pdf"))
save_csv("Fig1a_enzyme_sweep", dict(free_enzyme=E, **{"s%d_over_s%d" % (i+1, i): Rr[:, i] for i in range(n_a)}), "Fig 1a: n=4, non-identical ladder, fixed free effector l=0.9")
save_csv("Fig1b_deadend_sweep", dict(K_Z=xb, **{"s%d_over_s%d" % (i+1, i): yb[:, i] for i in range(n_b)}), "Fig 1b: n=3, T_E/T_L=1.1e-9")
save_json("fig1_invariance", dict(enzyme_fold_range=fold_a, ratio_drift_a=drift_a, max_abs_dev_from_Ki_over_l=dev_a,
                                  deadend_ratio_drift=drift_b, TE_over_TL_b=TE_b/TL_b))
