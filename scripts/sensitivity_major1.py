"""Major 1: propagate experimental error on the three measured cascade coefficients.

The residual's non-monotonicity in PII is inherited entirely from the measured
n_H^90 = 4.48, 6.46, 5.23, which carry no reported error.  Assign each an independent
relative error sigma and ask how often the middle residual still exceeds both ends.
The digitisation bootstrap of the first layer is included jointly.
"""
import _path  # noqa: F401
import os
from esbm import RESULTS, save_json

import numpy as np
rng = np.random.default_rng(1)
MEAS = np.array([4.48, 6.46, 5.23])
COMP = np.array([4.9145, 3.7794, 3.5518])          # point composites
res = np.load(os.path.join(RESULTS, 'boot_residuals.npy'))   # bootstrap residuals from cascade_budget.py
comp_boot = MEAS/res                                # implied composites per replicate
print("point residual: %s" % np.round(MEAS/COMP, 3))
print("margin: middle residual minus larger end = %.3f" % (MEAS[1]/COMP[1] - max(MEAS[0]/COMP[0], MEAS[2]/COMP[2])))
print("the middle measurement 6.46 would have to fall to %.2f (-%.1f%%) for the peak to vanish"
      % (COMP[1]*MEAS[2]/COMP[2], 100*(1-COMP[1]*MEAS[2]/COMP[2]/MEAS[1])))
table = {}
print("\n%6s  %14s  %14s" % ("sigma", "P(peak middle)", "P(any non-mono)"))
for sig in (0.0, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20):
    N = 200000
    k = rng.integers(0, len(comp_boot), N)
    m = MEAS*np.exp(sig*rng.standard_normal((N, 3)))   # lognormal error on each measured value
    r = m/comp_boot[k]
    peak = (r[:, 1] > r[:, 0]) & (r[:, 1] > r[:, 2])
    mono = ((r[:, 0] <= r[:, 1]) & (r[:, 1] <= r[:, 2])) | ((r[:, 0] >= r[:, 1]) & (r[:, 1] >= r[:, 2]))
    print("%6.2f  %14.3f  %14.3f" % (sig, peak.mean(), 1-mono.mean())); table[str(sig)] = [float(peak.mean()), float(1-mono.mean())]

save_json("sensitivity_major1", dict(margin=float(MEAS[1]/COMP[1] - max(MEAS[0]/COMP[0], MEAS[2]/COMP[2])),
                                     middle_to_vanish=float(COMP[1]*MEAS[2]/COMP[2]), drop_percent=float(100*(1-COMP[1]*MEAS[2]/COMP[2]/MEAS[1])),
                                     n_boot=len(res), P_peak_vs_sigma=table))
