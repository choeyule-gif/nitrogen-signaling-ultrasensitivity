"""Companion manuscript (Choi 2026, stability): the scalar condition (P), Re F(i omega) >= 0, across n.
Not a figure of the PLOS Computational Biology manuscript; kept because the Results cite the result.

F is the Laplace transform of the parity autocorrelation of the auxiliary Markov chain
of the Results.  Following the reduction, Re F = ||r||^2_pi + Im <r,h>_pi with
r = delta*phi - 1, delta = <phi>_pi, and h the solution of the projected linear system.
The chain is built from the steady state of the identical-independent family.
"""
import os, sys; sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from esbm import FIGURES, save_json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

def chain(n, A=1.0, Bh=0.4, Gh=0.4, Kh=1.0, l=1.0):
    """Auxiliary chain at the identical-independent steady state.  Order: S_0..S_n, B_*, G_*."""
    m = 3*n+1
    b2 = b3 = g2 = g3 = 1.0
    Bhi = [(n-i)*Bh for i in range(n)]
    Ghi = [(i+1)*Gh for i in range(n)]
    b1 = [Bhi[i]*(b2+b3) for i in range(n)]
    g1 = [Ghi[i]*(g2+g3) for i in range(n)]
    e = 1.0; eps = A*e*l
    S = lambda k: k; B = lambda k: (n+1)+k; G = lambda k: (n+1)+n+k
    L = np.zeros((m, m))                       # row generator: L[from, to]
    for i in range(n):
        L[S(i),   B(i)] += b1[i]*e
        L[B(i),   S(i)] += b2
        L[B(i), S(i+1)] += b3
        L[S(i+1), G(i)] += g1[i]*eps
        L[G(i), S(i+1)] += g2
        L[G(i),   S(i)] += g3
    np.fill_diagonal(L, -L.sum(1))
    # stationary distribution
    w, V = np.linalg.eig(L.T)
    pi = np.real(V[:, np.argmin(np.abs(w))]); pi = np.abs(pi); pi /= pi.sum()
    q = -np.diag(L).copy()
    c = float(pi @ q); L, q = L/c, q/c            # mean exit rate 1
    phi = np.array([1.0]*(n+1) + [-1.0]*(2*n))
    return L, q, pi, phi, m

def ReF(L, q, pi, phi, m, om):
    P = np.diag(1.0/q) @ (L + np.diag(q))
    nu = pi*q; Sm = (np.eye(m)+P)/2
    T = np.diag(om/q); d = float(pi@phi); r = d*phi - np.ones(m)
    R2 = 1.0 - d*d
    if om == 0: return R2
    Pp = np.eye(m) - np.outer(phi, nu*phi)
    Am = 2*Sm + 1j*(Pp @ T @ (np.eye(m) - np.outer(phi, pi*phi)))
    g = Pp @ (T @ r)
    Dn = np.diag(np.sqrt(nu)); Di = np.diag(1.0/np.sqrt(nu))
    pt = Dn@phi; pt = pt/np.linalg.norm(pt)
    Zt = np.linalg.qr(np.column_stack([pt] + list(np.eye(m).T)))[0][:, 1:m]
    h = Di @ (Zt @ np.linalg.solve(Zt.conj().T@(Dn@Am@Di)@Zt, Zt.conj().T@(Dn@g)))
    return R2 + float(np.imag(np.sum(pi*np.conj(r)*h)))

def circ_threshold(L, q, pi, m):
    Jm = np.outer(pi, np.ones(m))*L - (np.outer(pi, np.ones(m))*L).T
    Pi = np.diag(1.0/np.sqrt(pi))
    return 0.5*np.linalg.norm(Pi@Jm@Pi, 2)

NS = [1, 2, 3, 4, 6]
OM = np.logspace(-3, 3, 300)
Z = np.zeros((len(NS), len(OM))); QMAX = []; THR = []
for j, n in enumerate(NS):
    L, q, pi, phi, m = chain(n)
    QMAX.append(q.max()); THR.append(circ_threshold(L, q, pi, m))
    for i, om in enumerate(OM): Z[j, i] = ReF(L, q, pi, phi, m, om)
    print("n=%2d  F(0)=%.5f  q_max=%.3f  circulation threshold=%.3f  min Re F=%.3e at om=%.3g"
          % (n, Z[j,0], q.max(), THR[-1], Z[j].min(), OM[np.argmin(Z[j])]))
print("global min over the grid: %.3e" % Z.min())
print("threshold / q_max:", " ".join("%.3f" % (t/qm) for t, qm in zip(THR, QMAX)))

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 10, 'pdf.fonttype': 42})
fig = plt.figure(figsize=(4.6, 3.35))
ax = fig.add_axes([0.0, 0.03, 0.99, 0.95], projection='3d')
lo = np.log10(OM)
for j, n in enumerate(NS):
    k = np.searchsorted(OM, THR[j])
    # proved: omega >= circulation threshold, and all omega at n = 1
    ax.plot(lo[k:], [n]*(len(OM)-k), Z[j][k:], color='#0072BD', lw=1.7, zorder=6)
    if n == 1:
        ax.plot(lo[:k+1], [n]*(k+1), Z[j][:k+1], color='#0072BD', lw=1.7, zorder=6)
    else:
        ax.plot(lo[:k+1], [n]*(k+1), Z[j][:k+1], color='#D95319', lw=2.6, zorder=8)
    ax.plot([lo[0]], [n], [Z[j][0]], 'o', color='#0072BD', ms=3.6, zorder=9)
ax.plot(lo, [0.3]*len(OM), [0]*len(OM), color='0.55', lw=0.8)
ax.set_xlabel(r'$\log_{10}\omega$', labelpad=-3)
ax.set_ylabel(r'$n$', labelpad=-4)
ax.set_zlabel(r'$Re\,F$', labelpad=-4)
ax.set_yticks(NS); ax.set_xticks([-3, -1, 1, 3]); ax.set_zlim(0, 1.05)
ax.set_zticks([0, 0.5, 1.0])
ax.tick_params(pad=-1); ax.view_init(elev=20, azim=-121)
for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane): pane.set_alpha(0.0)
fig.savefig(os.path.join(FIGURES, "FigS_ReF_companion.pdf"))
save_json("extras_ReF", dict(global_min=float(Z.min()), F0={int(n): float(Z[j,0]) for j, n in enumerate(NS)}))
