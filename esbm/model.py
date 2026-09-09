"""The effector-switched bifunctional network, eqs (1)-(3) of the manuscript.

    E + L  <->  EL                       (effector binding, a1/a2, A = a1/a2)
    E + S_i  <->  B_i  ->  E + S_{i+1}   (forward, b1_i, b2_i, b3_i)
    EL + S_{i+1}  <->  G_i  ->  EL + S_i (reverse, g1_i, g2_i, g3_i)

Two independent routes to the steady state are provided and used against each other:
(a) the kernel of the compartmental matrix at fixed free enzyme e and free effector l
    (build_Q / kernel / steady_state), which imposes no closed form; and
(b) the closed forms proved in the paper (chain, Adair polynomial, quartic).
Constants: Bhat_i = b1_i/(b2_i+b3_i), Ghat_i = g1_i/(g2_i+g3_i), rho_i = b3_i/g3_i,
K_i = Bhat_{i-1} b3_{i-1} / (A Ghat_{i-1} g3_{i-1}) so that s_i/s_{i-1} = K_i/l.
"""
import numpy as np
from math import comb
from scipy.optimize import brentq

# ----------------------------------------------------------------------------- (a) raw network
def build_Q(n, rates, e, eps):
    """Compartmental matrix Q (columns sum to zero) on the 3n+1 substrate-bearing species
    S_0..S_n, B_0..B_{n-1}, G_0..G_{n-1}, at fixed free enzyme e and fixed [EL] = eps."""
    b1, b2, b3, g1, g2, g3 = rates
    m = 3*n + 1
    S = lambda i: i; B = lambda i: n + 1 + i; G = lambda i: 2*n + 1 + i
    Q = np.zeros((m, m))
    for i in range(n):
        Q[B(i), S(i)] += b1[i]*e;      Q[S(i), S(i)] -= b1[i]*e
        Q[S(i), B(i)] += b2[i];        Q[B(i), B(i)] -= b2[i]
        Q[S(i+1), B(i)] += b3[i];      Q[B(i), B(i)] -= b3[i]
        Q[G(i), S(i+1)] += g1[i]*eps;  Q[S(i+1), S(i+1)] -= g1[i]*eps
        Q[S(i+1), G(i)] += g2[i];      Q[G(i), G(i)] -= g2[i]
        Q[S(i), G(i)] += g3[i];        Q[G(i), G(i)] -= g3[i]
    return Q

def kernel(Q):
    """The unique positive vector z with Qz = 0 and sum z = 1 (exact linear solve)."""
    m = Q.shape[0]
    A = np.vstack([Q[:-1, :], np.ones(m)]); b = np.zeros(m); b[-1] = 1.0
    return np.linalg.solve(A, b)

def steady_state(n, rates, e, l, A):
    """Normalised steady state at fixed (e, l): returns (s, beta, gamma)."""
    z = kernel(build_Q(n, rates, e, A*e*l))
    return z[:n+1], z[n+1:2*n+1], z[2*n+1:]

def chain_constants(rates, A):
    b1, b2, b3, g1, g2, g3 = [np.asarray(x, float) for x in rates]
    Bh = b1/(b2 + b3); Gh = g1/(g2 + g3); rho = b3/g3
    K = Bh*rho/(A*Gh)                      # K_{i+1}, i = 0..n-1
    return Bh, Gh, rho, K

def random_rates(n, rng, decades=1.5):
    return [10**rng.uniform(-decades, decades, n) for _ in range(6)]

def identical_independent_rates(n, Kh, A, Bh=None, Gh=None, kf=1.0, kr=1.0, kon=1.0):
    """Identical, independent sites: K_{i+1} = (n-i)/(i+1) Khat.
    Returns rate constants with Bhat_i = (n-i) Bhat, Ghat_i = (i+1) Ghat and rho = kf/kr."""
    Bh = 1.0 if Bh is None else Bh
    Gh = (Bh*kf/(A*Kh*kr)) if Gh is None else Gh
    b2 = np.full(n, kon); b3 = np.full(n, kf); g2 = np.full(n, kon); g3 = np.full(n, kr)
    b1 = np.array([(n-i)*Bh*(b2[i]+b3[i]) for i in range(n)])
    g1 = np.array([(i+1)*Gh*(g2[i]+g3[i]) for i in range(n)])
    return b1, b2, b3, g1, g2, g3

# ----------------------------------------------------------------------------- (b) closed forms
def adair(c, u):
    """Fractional occupancy theta = (1/n) dlogP/dlogu of the ladder P(u) = sum c_i u^i."""
    c = np.asarray(c, float); n = len(c) - 1; i = np.arange(n+1)
    w = c*u**i; return (i*w).sum()/(n*w.sum())

def nu(c, u):
    """Ladder factor nu = n Var(i)/(<i>(n-<i>)) of the distribution p_i ∝ c_i u^i."""
    c = np.asarray(c, float); n = len(c) - 1; i = np.arange(n+1)
    w = c*u**i; p = w/w.sum(); m = (i*p).sum(); v = ((i-m)**2*p).sum()
    return n*v/(m*(n-m))

def binomial_ladder(n, Kh):
    return np.array([comb(n, k)*Kh**k for k in range(n+1)], float)

def theta_tot(l, TE, TS, n=3, A=1.0, Bh=0.4, Gh=0.4, Kh=1.0):
    """Total modified fraction (free plus enzyme-bound target) at identical independent
    sites, from the full steady state; the parameters of Fig 2c, Fig 3 and Fig 7."""
    C = Bh + A*Gh*Kh
    u = Kh/l; P = (1+u)**n; Pm = (1+u)**(n-1); Pm2 = (1+u)**(n-2)
    f = lambda e: e*(1+A*l) + n*e*(TS/(P+n*e*C*Pm))*C*Pm - TE
    hi = 1.0
    while f(hi) < 0: hi *= 10
    e = brentq(f, 1e-300, hi, xtol=1e-18, rtol=8.9e-16)
    s0 = TS/(P + n*e*C*Pm)
    M = s0*(n*u*Pm + Bh*e*n*(n-1)*u*Pm2 + A*Gh*Kh*e*n*Pm2*(1+n*u))
    return M/(n*TS)

def nH_local_tot(TE, TS, d=1e-5, **kw):
    """Local Hill coefficient -dlog[theta/(1-theta)]/dlog l of theta_tot at its half-point."""
    g = lambda L: theta_tot(np.exp(L), TE, TS, **kw) - 0.5
    L = brentq(g, np.log(1e-8), np.log(1e8), xtol=1e-14)
    lg = lambda t: np.log(t/(1-t))
    return -(lg(theta_tot(np.exp(L+d), TE, TS, **kw)) - lg(theta_tot(np.exp(L-d), TE, TS, **kw)))/(2*d)

def n1090(F, lo=1e-14, hi=1e14):
    """The range measure log81/log(S10/S90) of an increasing response F on [0,1]."""
    a = brentq(lambda x: F(x)-0.1, lo, hi); b = brentq(lambda x: F(x)-0.9, lo, hi)
    return np.log(81)/np.log(b/a)

def hill4(x, U0, Um, S, h):
    """Four-parameter Hill function, decreasing in x (uridylylation state against glutamine)."""
    return Um + (U0-Um)/(1+(x/S)**h)
