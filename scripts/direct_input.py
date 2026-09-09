"""'What is left': the first candidate, a direct glutamine input on the adenylyltransferase.

Re-implementation.  The original checks were run interactively and their scripts were not kept;
the manuscript quotes 3.0009 (PII path) + 0.4913 (direct path) = 3.4922 (computed directly) for the
additivity of the two log-slopes, and a ten-to-ninety gain 'close to one and usually below it',
0.95-1.2 at the measured adenylyltransferase constants.  This script states its construction:

Chain variable of the second layer: f = l_2/(l_1 phi(G)), with l_1 = (1-theta(G))^3 the fully
unmodified PII, l_2 = theta(G) the linear reading of the reverse effector, and phi(G) the direct
modulation of the forward/reverse ratio by glutamine acting on a single site with affinity K_G
and fold-activation R of the forward activity (phi = (1 + R G/K_G)/(1 + G/K_G)).
1. Additivity: dlog f/dlog G = [dlog(l_2/l_1)/dlog G] + [-dlog phi/dlog G], checked at the
   half-point of the composite.
2. Gain: the n_H^90 of the composite with and without the direct term, with Khat_2 pinned where
   the PII path alone puts it, maximised over K_G and R over four decades each, and evaluated at
   R = 3.7 (Jiang et al. 1998b) across K_G.
"""
import _path  # noqa: F401
import numpy as np
from scipy.optimize import brentq
from esbm import save_json, load_json
from esbm.model import n1090
cb = load_json("cascade_budget")['fit']
th0, thi, S, h = cb['theta0'], cb['theta_inf'], cb['S_mM'], cb['h']
theta = lambda G: thi + (th0-thi)/(1+(G/S)**h)
fPII = lambda G: theta(G)/(1-theta(G))**3
def phi(G, KG, R): return (1+R*G/KG)/(1+G/KG)
def composite(KG, R, K2):
    f = lambda G: fPII(G)/phi(G, KG, R)
    f0, fi = f(1e-12), f(1e12)
    return lambda G: ((f0-f(G))/(f0-fi))*((K2+fi)/(K2+f(G)))
# Khat_2 pinned where the PII path alone puts the half-point at the measured 0.530 mM (middle PII)
hp = 0.530
K2 = np.exp(brentq(lambda lk: np.log(brentq(lambda G: composite(1.0, 1.0, np.exp(lk))(G)-0.5, 1e-12, 1e12)/hp), -45, 45))
base = n1090(composite(1.0, 1.0, K2), 1e-12, 1e12)
print("PII path alone: Khat_2 = %.4g, n_H^90 = %.4f" % (K2, base))
# 1. additivity of the log-slopes at the composite's half-point, at the measured R = 3.7 and K_G = S
KG, R = S, 3.7
Gh_ = brentq(lambda G: composite(KG, R, K2)(G)-0.5, 1e-12, 1e12); d = 1e-5
ls = lambda F, G: (np.log(F(G*(1+d)))-np.log(F(G*(1-d))))/(2*d)
a1 = -ls(fPII, Gh_); a2 = ls(lambda G: phi(G, KG, R), Gh_); tot = -ls(lambda G: fPII(G)/phi(G, KG, R), Gh_)
print("additivity at the half-point: PII path %.4f + direct %.4f = %.4f ; computed directly %.4f" % (a1, a2, a1+a2, tot))
# 2. gain in the ten-to-ninety measure
best = (0, None)
grid = {}
for KG in np.logspace(np.log10(S)-2, np.log10(S)+2, 17):
    for RR in np.logspace(0, 4, 17):
        try: gval = n1090(composite(KG, RR, K2), 1e-12, 1e12)/base
        except Exception: continue
        grid[(KG, RR)] = gval
        if gval > best[0]: best = (gval, (KG, RR))
vals = np.array(list(grid.values()))
print("gain over four decades of K_G and R: min %.3f, max %.3f (at K_G=%.3g, R=%.3g); fraction below 1: %.2f" % (vals.min(), best[0], best[1][0], best[1][1], (vals < 1).mean()))
atR = {}
for KG in np.logspace(np.log10(S)-2, np.log10(S)+2, 9):
    atR[float(KG)] = n1090(composite(KG, 3.7, K2), 1e-12, 1e12)/base
print("at R = 3.7: gain from %.3f to %.3f across K_G" % (min(atR.values()), max(atR.values())))
save_json("direct_input", dict(K2=K2, base=base, additivity=dict(PII=a1, direct=a2, sum=a1+a2, direct_total=tot),
                               gain_min=vals.min(), gain_max=best[0], gain_argmax=best[1], fraction_below_1=float((vals < 1).mean()),
                               gain_at_R3p7=atR))
