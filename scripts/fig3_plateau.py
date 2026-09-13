"""Fig 3: the plateau of eq (zeroorder) as a surface.

Local Hill coefficient of the TOTAL modified fraction over both totals at once:
T_E/T_S on one axis, T_S/K_m on the other, with K_m = 1/(n*Bhat), at n = 3, A*Khat = 1,
Bhat = Ghat = 0.4/Khat.  Reports the range over the grid and over the interior
T_E/T_S <= 1e-2 (the manuscript quotes 0.99745 to 1, width 2.5e-3).
"""
import _path  # noqa: F401  (makes esbm importable when run as a script)
import os, numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from esbm import FIGURES, RESULTS, save_json, save_csv
from esbm.model import nH_local_tot

n, Bh = 3, 0.4
Km = 1.0/(n*Bh)
nx, ny = 46, 46
X = np.logspace(-4, 1.3, nx)                 # T_E / T_S
Y = np.logspace(-2, 3, ny)                   # T_S / K_m
Z = np.empty((ny, nx))
for j, ts in enumerate(Y):
    TS = ts*Km
    for i, r in enumerate(X):
        Z[j, i] = nH_local_tot(r*TS, TS)
interior = Z[:, X <= 1e-2]
print("plateau: min %.5f max %.5f over the grid" % (Z.min(), Z.max()))
print("for T_E/T_S <= 1e-2 : %.6f to %.6f  (range %.1e)" % (interior.min(), interior.max(), interior.max()-interior.min()))
np.savez(os.path.join(RESULTS, "plateau_grid.npz"), X=X, Y=Y, Z=Z)
LXf, LYf = np.meshgrid(X, Y)
save_csv("Fig3_plateau_surface", dict(TE_over_TS=LXf.ravel(), TS_over_Km=LYf.ravel(), nH_local=Z.ravel()), "S3 Fig: n=3, A*Khat=1, Bhat=Ghat=0.4/Khat, K_m=1/(n Bhat)")

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 10, 'pdf.fonttype': 42})
fig = plt.figure(figsize=(4.4, 3.25))
ax = fig.add_axes([0.03, 0.05, 0.94, 0.92], projection='3d')
LX, LY = np.meshgrid(np.log10(X), np.log10(Y))
ax.plot_surface(LX, LY, Z, cmap='viridis', rstride=1, cstride=1, lw=0.25, edgecolor='0.35', antialiased=True, vmin=0.6, vmax=1.0)
ax.contour(LX, LY, Z, levels=[0.80, 0.90, 0.95, 0.99], colors='k', linewidths=0.7, offset=0.6)
ax.set_xlabel(r'$\log_{10}(T_E/T_S)$', labelpad=0); ax.set_ylabel(r'$\log_{10}(T_S/K_m)$', labelpad=0)
ax.set_zlabel(r'$n_H^{loc}$', labelpad=-2)
ax.set_zlim(0.6, 1.02); ax.set_zticks([0.6, 0.8, 1.0]); ax.set_xticks([-4, -2, 0]); ax.set_yticks([-2, 0, 2])
ax.tick_params(pad=-1)
ax.view_init(elev=18, azim=-60)
ax.xaxis.pane.set_alpha(0.0); ax.yaxis.pane.set_alpha(0.0); ax.zaxis.pane.set_alpha(0.0)
fig.savefig(os.path.join(FIGURES, "Fig3.pdf"), metadata={"CreationDate": None})
save_json("fig3_plateau", dict(grid_min=Z.min(), grid_max=Z.max(), interior_min=interior.min(), interior_max=interior.max(),
                               interior_width=interior.max()-interior.min()))
