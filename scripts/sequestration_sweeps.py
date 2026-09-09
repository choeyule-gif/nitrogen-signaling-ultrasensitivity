"""eq (zeroorder): what survives of a dependence on the totals is the sequestration correction.

1. Worked case (parameters of Fig 2c): apparent coefficient of the total modified fraction at
   T_E/T_S = 0.1 and after doubling the enzyme (manuscript: 0.988 -> 0.976).
2. Sweep of Bhat = Ghat over two decades and A over three decades at T_E/T_S = 0.1 and 0.2,
   n = 1 and 3 (manuscript: coefficient between 0.92 and 0.999, change on doubling 0.001-0.08).
3. Sweep of the target from K_m/60 to 100 K_m at the enzyme level of the titration
   (manuscript: n_H between 0.981 and 0.9996).
"""
import _path  # noqa: F401
import numpy as np
from esbm import save_json
from esbm.model import nH_local_tot
out = {}
print("1. worked case (n=3, A=1, Bhat=Ghat=0.4, Khat=1, T_S=1)")
a = nH_local_tot(0.1, 1.0); b = nH_local_tot(0.2, 1.0)
print("   T_E/T_S=0.1: %.5f ; 0.2: %.5f ; change %.5f" % (a, b, a-b))
out['worked'] = dict(nH_0p1=a, nH_0p2=b, change=a-b)
def sweep(As, TS, TE, label):
    rows = []; vals = []; chg = []
    for n in (1, 3):
        for B in (0.2, 2.0, 20.0):
            for A_ in As:
                kw = dict(n=n, A=A_, Bh=B, Gh=B, Kh=1.0)
                x = nH_local_tot(TE, TS, **kw); y = nH_local_tot(2*TE, TS, **kw)
                rows.append(dict(n=n, Bhat=B, A=A_, nH=x, nH_doubled=y, change=abs(x-y)))
                vals += [x, y]; chg.append(abs(x-y))
    print("   %s: coefficient range [%.3f, %.4f]; change on doubling [%.4f, %.3f]" % (label, min(vals), max(vals), min(chg), max(chg)))
    return dict(rows=rows, nH_min=min(vals), nH_max=max(vals), change_min=min(chg), change_max=max(chg))
print("2. sweeps of Bhat = Ghat (two decades) and A, at T_E/T_S = 0.1 and 0.2")
out['sweep_manuscript_design'] = sweep((1.0, 29.0), 0.5, 0.05, "A in {1, 29}, T_S = 0.5, T_E = 0.05 -> 0.10")
out['sweep_three_decades_A'] = sweep((0.1, 1.0, 10.0, 100.0), 1.0, 0.1, "A over three decades, T_S = 1, T_E = 0.1 -> 0.2")
print("3. target sweep K_m/60 .. 100 K_m at fixed T_E = 0.1 (K_m = 1/(n Bhat))")
n, Bh = 3, 0.4; Km = 1.0/(n*Bh); TE = 0.1
sw = [(TS, nH_local_tot(TE, TS)) for TS in np.logspace(np.log10(Km/60), np.log10(100*Km), 60)]
lo = min(v for _, v in sw); hi = max(v for _, v in sw)
print("   K_m = %.4f ; n_H in [%.4f, %.4f] (width %.4f)" % (Km, lo, hi, hi-lo))
out['target_sweep'] = dict(Km=Km, nH_min=lo, nH_max=hi, width=hi-lo)
save_json("sequestration_sweeps", out)
