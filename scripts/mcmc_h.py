"""Posterior of the exponent h: adaptive Metropolis, four chains of 1e5, first fifth discarded
(Methods).  Flat priors on U0, Umin, log K, h; Jeffreys prior on the error scale.
Reports acceptance, Gelman-Rubin, effective sample size, median, 95% credible interval,
P(h < 1.5) and P(h < 1).  Writes results/mcmc_h.npy (thinned) for fig5_integer.py."""
import _path  # noqa: F401
import os, numpy as np
from esbm import RESULTS, save_json
from esbm.data import digitised_titration
rng = np.random.default_rng(20260901)
g, U = digitised_titration(); n = len(U)
def logpost(t):
    U0, Um, lK, h, ls = t
    if not (0 < h < 8 and -12 < lK < 12 and -8 < ls < 4 and 0 < U0 < 6 and -1 < Um < 3 and Um < U0): return -np.inf
    s = np.exp(ls); r = U-(Um+(U0-Um)*np.exp(lK)/(np.exp(lK)+g**h))
    return -n*ls-0.5*np.sum(r*r)/s**2
NS = 100_000; BURN = 20_000; NCH = 4
start = np.array([3.0, 0.45, np.log(0.313), 2.018, np.log(0.041)])
chains = []; accs = []
for c in range(NCH):
    t = start+rng.normal(0, [0.05, 0.05, 0.15, 0.08, 0.15])
    while not np.isfinite(logpost(t)): t = start+rng.normal(0, [0.02, 0.02, 0.05, 0.03, 0.05])
    lp = logpost(t); Cn = np.diag([0.02, 0.02, 0.06, 0.03, 0.06])**2; sd = 2.38**2/5
    S = np.empty((NS, 5)); acc = 0
    for i in range(NS):
        p = t+rng.multivariate_normal(np.zeros(5), sd*Cn+1e-12*np.eye(5))
        lq = logpost(p)
        if np.log(rng.random()) < lq-lp: t, lp = p, lq; acc += 1
        S[i] = t
        if i >= 500 and i % 50 == 0:           # Haario adaptation on a trailing window
            W = S[max(0, i-5000):i+1]; Cn = np.cov(W.T)+1e-10*np.eye(5)
    chains.append(S[BURN:]); accs.append(acc/NS)
    print("  chain %d: acceptance %.3f" % (c+1, acc/NS))
X = np.concatenate(chains); h = X[:, 3]
def R_hat(ch, k):
    L = len(ch[0]); mn = np.array([c[:, k].mean() for c in ch]); vr = np.array([c[:, k].var(ddof=1) for c in ch])
    B = L*mn.var(ddof=1); W = vr.mean(); V = (L-1)/L*W+B/L
    return np.sqrt(V/W)
def ESS(x):
    x = x-x.mean(); N = len(x); f = np.fft.rfft(x, 2*N); ac = np.fft.irfft(f*np.conj(f))[:N].real; ac /= ac[0]
    s = 0.0
    for k in range(1, min(N, 20000)):
        if ac[k] < 0.05: break
        s += ac[k]
    return N/(1+2*s)
lo, hi = np.percentile(h, [2.5, 97.5])
res = dict(acceptance_mean=float(np.mean(accs)), R_hat=R_hat(chains, 3), ESS=ESS(h), median=float(np.median(h)),
           mean=h.mean(), sd=h.std(), ci95=[lo, hi], P_h_lt_1p5=float(np.mean(h < 1.5)), P_h_lt_1=float(np.mean(h < 1)),
           n_samples=len(h))
print("  acceptance %.3f  R-hat %.5f  ESS %.0f" % (res['acceptance_mean'], res['R_hat'], res['ESS']))
print("  median %.4f  95%% CI [%.3f, %.3f]  P(h<1.5) = %.1e  P(h<1) = %.1e" % (res['median'], lo, hi, res['P_h_lt_1p5'], res['P_h_lt_1']))
np.save(os.path.join(RESULTS, "mcmc_h.npy"), h[::20])
save_json("mcmc_h", res)
