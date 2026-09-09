"""Fig 6: the residual R(l) whose zeros are the positive steady states.  Run counterexamples.py first.

(a) both conditions of eq (cond) satisfied (n = 3 identical independent sites, Bhat = 1.5,
    Ghat = 0.8, T_E = 1, T_S = 4, T_L = 2): one zero, R increasing throughout;
(b) condition (i) broken by a spread of 2.89 in the catalytic ratios (S1_Data/S1_counterexamples.json):
    three zeros.
"""
import _path  # noqa: F401
import os, json, numpy as np
from math import comb
import matplotlib.pyplot as plt
from esbm import FIGURES, S1DATA, save_json, save_csv
from esbm.style import MC, use, tidy, panel
use()
def resid(n, A, c, Bh, Gh, TE, TS, TL, ls):
    u = 1.0/ls; ar = np.arange(n+1)[:, None]; ar2 = np.arange(n)[:, None]
    Mp = (c[:, None]*u**ar).sum(0); Vt = ((Gh*c[1:])[:, None]*u**ar2).sum(0)
    Q = ((Bh*c[:n])[:, None]*u**ar2).sum(0)+A*Vt
    a = (1+A*ls)*Q; b = (1+A*ls)*Mp+Q*TS-TE*Q; cc = -TE*Mp
    e = (-b+np.sqrt(np.maximum(b*b-4*a*cc, 0)))/(2*a); s0 = TS/(Mp+e*Q)
    return ls+A*e*ls+A*e*s0*Vt-TL, e, s0
ex = json.load(open(os.path.join(S1DATA, "S1_counterexamples.json")))["condition_i"]
n, A, c, Bh, Gh, TE, TS, TL = ex['n'], ex['A'], np.array(ex['c']), np.array(ex['Bhat']), np.array(ex['Ghat']), ex['TE'], ex['TS'], ex['TL']
rts = [float(r['l']) for r in ex['roots_60digit']]

fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
nA = 3; Kh = 1.0; AA = 1.0
cc0 = np.array([comb(nA, k)*Kh**k for k in range(nA+1)], float)
Bh0 = np.full(nA, 1.5); Gh0 = np.full(nA, 0.8); TE0, TS0, TL0 = 1.0, 4.0, 2.0
lg = np.logspace(-3, 3, 1500)
R0, e0, s00 = resid(nA, AA, cc0, Bh0, Gh0, TE0, TS0, TL0, lg)
m0 = np.isfinite(R0) & (e0 > 0) & (s00 > 0)
ax[0].plot(lg[m0], R0[m0], '-', color=MC[0]); ax[0].axhline(0, color='k', ls='--', lw=1.0)
i0 = np.where(np.diff(np.sign(R0[m0])) != 0)[0]
ax[0].plot(lg[m0][i0], np.zeros(len(i0)), 'o', color=MC[1], ms=7, zorder=5)
ax[0].set_xscale('log'); ax[0].set_xlabel(r'$l$'); ax[0].set_ylabel(r'$\mathcal{R}(l)$')
ax[0].set_xlim(1e-2, 1e2); ax[0].set_ylim(-1.5, 6); tidy(ax[0], 1, 0, nx=3); panel(ax[0], 'a', -0.20)
print("(a) zeros: %d ; R increasing: %s" % (len(i0), bool(np.all(np.diff(R0[m0]) > 0))))
lo, hi = min(rts)/40, max(rts)*40
lg2 = np.logspace(np.log10(lo), np.log10(hi), 3000)
R2, e2, s2 = resid(n, A, c, Bh, Gh, TE, TS, TL, lg2); m2 = np.isfinite(R2) & (e2 > 0) & (s2 > 0)
ax[1].plot(lg2[m2], R2[m2], '-', color=MC[0]); ax[1].axhline(0, color='k', ls='--', lw=1.0)
ax[1].plot(rts, np.zeros(len(rts)), 'o', color=MC[1], ms=7, zorder=5)
ax[1].set_xscale('log'); ax[1].set_xlabel(r'$l$'); ax[1].set_ylabel(r'$\mathcal{R}(l)$')
ax[1].set_xlim(lo, hi); tidy(ax[1], 1, 0, nx=3); panel(ax[1], 'b', -0.20)
plt.tight_layout(w_pad=2.0); plt.savefig(os.path.join(FIGURES, 'Fig6.pdf'), metadata={"CreationDate": None})
save_csv("Fig6a_residual", dict(l=lg[m0], R=R0[m0]), "Fig 6a: n=3, Bhat=1.5, Ghat=0.8, T_E=1, T_S=4, T_L=2")
save_csv("Fig6b_residual", dict(l=lg2[m2], R=R2[m2]), "Fig 6b: parameters in S1_counterexamples.json (condition_i)")
print("(b) n=%d spread %.2f zeros %d" % (n, ex['spread'], len(rts)))
save_json("fig6_uniqueness", dict(panel_a_zeros=len(i0), panel_a_increasing=bool(np.all(np.diff(R0[m0]) > 0)), panel_b_n=n, panel_b_spread=ex['spread'], panel_b_zeros=len(rts)))
