"""Minor: the reassembly-mixture argument for R1, without assuming the answer.

The monovalent construct is 15% di/trivalent by competent trimer (27% of modifiable
subunits).  The published argument that this does not bias the ratio assumes the ladder
is binomial -- which is what is being measured.  Here the ladder is given a cooperativity
w, the mixture's apparent coefficient is computed, and the bias on nu_PII read off.
"""
import _path  # noqa: F401
from esbm import save_json

import numpy as np
from scipy.optimize import brentq
from math import comb

p = 6.0/42.0                                    # wild-type subunit fraction
wts = np.array([comb(3, k)*p**k*(1-p)**(3-k) for k in range(4)])
wts = wts[1:]/wts[1:].sum()                     # competent trimers, k = 1,2,3
subun = np.array([k*wts[k-1] for k in range(1, 4)]); subun /= subun.sum()
print("competent trimers  k=1,2,3 : %s   (di+tri %.1f%%)" % (np.round(wts,4), 100*wts[1:].sum()))
print("modifiable subunits k=1,2,3: %s   (di+tri %.1f%%)" % (np.round(subun,4), 100*subun[1:].sum()))

def theta_k(u, k, w):
    """mean occupancy of a k-site ladder with K_i = (k-i+1)/i * w^(i-1)."""
    K = np.array([(k-i+1)/i*w**(i-1) for i in range(1, k+1)])
    c = np.concatenate([[1.0], np.cumprod(K)])
    q = c*u**np.arange(k+1); q = q/q.sum()
    return (q*np.arange(k+1)).sum()/k

def n1090(F):
    a = brentq(lambda t: F(np.exp(t))-0.1, -60, 60)
    b = brentq(lambda t: F(np.exp(t))-0.9, -60, 60)
    return np.log(81)/abs(b-a)*np.log(81)/np.log(81)*1.0 if False else np.log(81)/abs(b-a)

rows = []
print("\n%6s %10s %12s %12s %12s" % ("w", "nu(n=3)", "pure n=1", "mixture", "bias on nu_PII"))
for w in (1.0, 1.1, 1.25, 1.5, 2.0):
    # nu of the wild-type trimer at its half-point
    f3 = lambda u: theta_k(u, 3, w)
    u50 = brentq(lambda t: f3(np.exp(t))-0.5, -60, 60)
    K = np.array([(3-i+1)/i*w**(i-1) for i in range(1, 4)])
    c = np.concatenate([[1.0], np.cumprod(K)]); q = c*np.exp(u50)**np.arange(4); q/=q.sum()
    i = np.arange(4); m=(q*i).sum(); v=(q*(i-m)**2).sum(); nu3 = 3*v/(m*(3-m))
    pure = n1090(lambda u: theta_k(u, 1, w))
    mix  = n1090(lambda u: sum(subun[k-1]*theta_k(u, k, w) for k in range(1, 4)))
    print("%6.2f %10.4f %12.4f %12.4f %12.4f" % (w, nu3, pure, mix, mix/pure)); rows.append(dict(w=w, nu3=float(nu3), pure=float(pure), mixture=float(mix), bias=float(mix/pure)))

save_json("mixture_bias", dict(competent_trimers=wts.tolist(), modifiable_subunits=subun.tolist(), di_tri_trimer_percent=float(100*wts[1:].sum()),
                               di_tri_subunit_percent=float(100*subun[1:].sum()), rows=rows))
