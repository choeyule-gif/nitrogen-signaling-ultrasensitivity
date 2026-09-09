"""Numerical checks against the raw network (no closed form imposed).

1. eq (adair): theta from the kernel steady state equals (1/n) dlogP/dlogu with
   P(u) = sum c_i u^i, c_i = prod K_j, over 1,500 random parameter sets (n = 1..6, six rate
   constants per site over three decades, e and l over four decades).  Reports the maximum
   relative deviation (manuscript: 1e-15).
2. eq (two)/(ratio): the two-form network (forward on E.L1, reverse on E.L2), solved as the
   kernel of its compartmental matrix for n <= 4, recovers s_i/s_{i-1} = K_i l_1/l_2 and the
   ratios do not move over a thousandfold range of free (hence total) enzyme
   (manuscript: 8e-16).
3. eq (integer): with eps = A e l^h the local Hill coefficient of theta_free is h exactly.
4. At n = 2 with K_1 = 1e-3 and K_2 = 1e3, sup_u nu = 1.999.
"""
import _path  # noqa: F401
import numpy as np
from scipy.optimize import minimize_scalar
from esbm import save_json
from esbm.model import build_Q, kernel, chain_constants, adair, nu, random_rates, identical_independent_rates
rng = np.random.default_rng(1500)

print("1. Adair form against the kernel steady state")
worst = 0.0; count = 0; rejected = 0
while count < 1500:
    n = int(rng.integers(1, 7)); rates = random_rates(n, rng, 1.5); A = 10**rng.uniform(-1, 1)
    e = 10**rng.uniform(-2, 2); l = 10**rng.uniform(-2, 2)
    z = kernel(build_Q(n, rates, e, A*e*l)); s = z[:n+1]
    if z.min() <= 0 or z.max()/z.min() > 1e8:   # outside this window double precision loses the cancellations (Methods)
        rejected += 1; continue
    theta_raw = (np.arange(n+1)*s).sum()/(n*s.sum())
    Bh, Gh, rho, K = chain_constants(rates, A)
    c = np.concatenate([[1.0], np.cumprod(K)])
    theta_adair = adair(c, 1.0/l)
    ratios = s[1:]/s[:-1]
    dev = max(abs(theta_raw/theta_adair-1), np.max(np.abs(ratios*l/K-1)))
    worst = max(worst, dev); count += 1
print("   %d parameter sets (dynamic range < 1e8; %d rejected); worst relative deviation %.2e" % (count, rejected, worst))

# the same check on a well-conditioned family: rate constants over one decade, e and l over two,
# which is where double precision resolves the identity to rounding (the wide family above is
# limited by the condition number of the linear solve, of order 1e6)
worst_n = 0.0; count_n = 0; rej_n = 0
while count_n < 1500:
    n = int(rng.integers(1, 7)); rates = random_rates(n, rng, 0.5); A = 10**rng.uniform(-0.5, 0.5)
    e = 10**rng.uniform(-1, 1); l = 10**rng.uniform(-1, 1)
    z = kernel(build_Q(n, rates, e, A*e*l)); s = z[:n+1]
    if z.min() <= 0 or z.max()/z.min() > 1e6: rej_n += 1; continue
    theta_raw = (np.arange(n+1)*s).sum()/(n*s.sum())
    Bh, Gh, rho, K = chain_constants(rates, A); c = np.concatenate([[1.0], np.cumprod(K)])
    dev = max(abs(theta_raw/adair(c, 1.0/l)-1), np.max(np.abs((s[1:]/s[:-1])*l/K-1)))
    worst_n = max(worst_n, dev); count_n += 1
print("   well-conditioned family (rates over one decade, e and l over two, dynamic range < 1e6): %d sets, worst %.2e" % (count_n, worst_n))
# the same 1,500-set check with the kernel solved in 40-digit arithmetic (mpmath): the identity is
# then resolved far below double precision, which separates the algebra from the conditioning
from mpmath import mp
mp.dps = 40
def kernel_mp(n, rates, e, A, l, A2=None, l2=None):
    """Kernel of the compartmental matrix in mp arithmetic; eps = A e l formed in mp.
    With A2, l2 given, the two-form network (forward on E.L1, reverse on E.L2)."""
    b1, b2, b3, g1, g2, g3 = rates; m = 3*n+1
    S = lambda i: i; B = lambda i: n+1+i; G = lambda i: 2*n+1+i
    Q = mp.zeros(m, m); e = mp.mpf(e)
    if A2 is None:
        e_f, eps = e, mp.mpf(A)*e*mp.mpf(l)
    else:
        e_f, eps = mp.mpf(A)*e*mp.mpf(l), mp.mpf(A2)*e*mp.mpf(l2)
    for i in range(n):
        Q[B(i), S(i)] += mp.mpf(b1[i])*e_f; Q[S(i), S(i)] -= mp.mpf(b1[i])*e_f
        Q[S(i), B(i)] += mp.mpf(b2[i]);     Q[B(i), B(i)] -= mp.mpf(b2[i])
        Q[S(i+1), B(i)] += mp.mpf(b3[i]);   Q[B(i), B(i)] -= mp.mpf(b3[i])
        Q[G(i), S(i+1)] += mp.mpf(g1[i])*eps; Q[S(i+1), S(i+1)] -= mp.mpf(g1[i])*eps
        Q[S(i+1), G(i)] += mp.mpf(g2[i]);   Q[G(i), G(i)] -= mp.mpf(g2[i])
        Q[S(i), G(i)] += mp.mpf(g3[i]);     Q[G(i), G(i)] -= mp.mpf(g3[i])
    Aug = mp.zeros(m, m); rhs = mp.zeros(m, 1)
    for r in range(m-1):
        for c_ in range(m): Aug[r, c_] = Q[r, c_]
    for c_ in range(m): Aug[m-1, c_] = 1
    rhs[m-1] = 1
    z = mp.lu_solve(Aug, rhs)
    return [z[i] for i in range(n+1)]
rng2 = np.random.default_rng(1500); worst_mp = mp.mpf(0)
for t in range(1500):
    n = int(rng2.integers(1, 7)); rates = random_rates(n, rng2, 1.5); A = 10**rng2.uniform(-1, 1)
    e = 10**rng2.uniform(-2, 2); l = 10**rng2.uniform(-2, 2)
    s = kernel_mp(n, rates, e, A, l)
    b1, b2, b3, g1, g2, g3 = rates
    K = [mp.mpf(b1[i])/(mp.mpf(b2[i])+mp.mpf(b3[i]))*mp.mpf(b3[i])/(mp.mpf(A)*mp.mpf(g1[i])/(mp.mpf(g2[i])+mp.mpf(g3[i]))*mp.mpf(g3[i])) for i in range(n)]
    c = [mp.mpf(1)]
    for i in range(n): c.append(c[-1]*K[i])
    u = 1/mp.mpf(l); w = [c[i]*u**i for i in range(n+1)]
    theta_adair = sum(i*w[i] for i in range(n+1))/(n*sum(w)); theta_raw = sum(i*s[i] for i in range(n+1))/(n*sum(s))
    dev = max(abs(theta_raw/theta_adair-1), max(abs(s[i+1]/s[i]*mp.mpf(l)/K[i]-1) for i in range(n)))
    worst_mp = max(worst_mp, dev)
print("   the same wide family with the kernel solved in 40-digit arithmetic: 1500 sets, worst %s" % mp.nstr(worst_mp, 3))
worst2_mp = mp.mpf(0); drift2_mp = mp.mpf(0)
for n in (1, 2, 3, 4):
    for t in range(25):
        rates = random_rates(n, rng2, 1.0); A1, A2 = 10**rng2.uniform(-1, 1, 2); l1, l2 = 10**rng2.uniform(-1, 1, 2)
        b1, b2, b3, g1, g2, g3 = rates
        Kp = [mp.mpf(b1[i])/(mp.mpf(b2[i])+mp.mpf(b3[i]))*mp.mpf(b3[i])*mp.mpf(A1)/(mp.mpf(g1[i])/(mp.mpf(g2[i])+mp.mpf(g3[i]))*mp.mpf(g3[i])*mp.mpf(A2))*(mp.mpf(l1)/mp.mpf(l2)) for i in range(n)]
        Rm = []
        for e in np.logspace(-2, 1, 7):
            s = kernel_mp(n, rates, e, A1, l1, A2, l2); Rm.append([s[i+1]/s[i] for i in range(n)])
        worst2_mp = max(worst2_mp, max(abs(Rm[j][i]/Kp[i]-1) for j in range(len(Rm)) for i in range(n)))
        drift2_mp = max(drift2_mp, max(abs(Rm[j][i]/Rm[0][i]-1) for j in range(len(Rm)) for i in range(n)))
print("   two-form network in 40-digit arithmetic (n <= 4, 100 sets, thousandfold enzyme): worst |ratio/prediction - 1| = %s, drift %s" % (mp.nstr(worst2_mp, 3), mp.nstr(drift2_mp, 3)))
print("2. two-form network, n <= 4, over a thousandfold range of enzyme")
def build_two(n, rates, e, A1, A2, l1, l2):
    b1, b2, b3, g1, g2, g3 = rates; eps1, eps2 = A1*e*l1, A2*e*l2
    m = 3*n+1; S = lambda i: i; B = lambda i: n+1+i; G = lambda i: 2*n+1+i
    Q = np.zeros((m, m))
    for i in range(n):
        Q[B(i), S(i)] += b1[i]*eps1;   Q[S(i), S(i)] -= b1[i]*eps1
        Q[S(i), B(i)] += b2[i];        Q[B(i), B(i)] -= b2[i]
        Q[S(i+1), B(i)] += b3[i];      Q[B(i), B(i)] -= b3[i]
        Q[G(i), S(i+1)] += g1[i]*eps2; Q[S(i+1), S(i+1)] -= g1[i]*eps2
        Q[S(i+1), G(i)] += g2[i];      Q[G(i), G(i)] -= g2[i]
        Q[S(i), G(i)] += g3[i];        Q[G(i), G(i)] -= g3[i]
    return Q
worst2 = 0.0; drift2 = 0.0
for n in (1, 2, 3, 4):
    for t in range(50):
        rates = random_rates(n, rng, 0.5); A1, A2 = 10**rng.uniform(-0.5, 0.5, 2); l1, l2 = 10**rng.uniform(-0.5, 0.5, 2)
        b1, b2, b3, g1, g2, g3 = rates
        Kpred = (b1/(b2+b3))*b3*A1/((g1/(g2+g3))*g3*A2)*(l1/l2)
        R = []; ok = True
        for e in np.logspace(-2, 1, 13):                     # thousandfold
            z = kernel(build_two(n, rates, e, A1, A2, l1, l2)); s = z[:n+1]
            if z.min() <= 0 or z.max()/z.min() > 1e8: ok = False; break
            R.append(s[1:]/s[:-1])
        if not ok: continue
        R = np.array(R)
        worst2 = max(worst2, np.max(np.abs(R/Kpred-1))); drift2 = max(drift2, np.max(np.abs(R/R[0]-1)))
print("   worst |ratio/prediction - 1| = %.2e ; worst drift over 1000-fold enzyme = %.2e" % (worst2, drift2))

print("3. eps = A e l^h gives local Hill coefficient h (n = 4, identical independent sites)")
n = 4; rates = identical_independent_rates(n, Kh=1.0, A=1.0, Bh=1.0, kf=137., kr=4., kon=1e4); A = 1/29.
res3 = {}
for h in (1, 2, 3):
    ls = np.geomspace(50, 5000, 25); th = []
    for l in ls:
        z = kernel(build_Q(n, rates, 1.0, A*l**h)); s = z[:n+1]
        th.append((np.arange(n+1)*s).sum()/(n*s.sum()))
    th = np.array(th); slope = np.gradient(np.log(th/(1-th)), np.log(ls))
    res3[h] = [float(-slope.max()), float(-slope.min())]
    print("   h=%d: local n_H in [%.10f, %.10f]" % (h, -slope.max(), -slope.min()))

print("4. nu at n = 2 with K_1 = 1e-3, K_2 = 1e3")
c = [1.0, 1e-3, 1.0]
r = minimize_scalar(lambda lu: -nu(c, np.exp(lu)), bounds=(-20, 20), method='bounded', options=dict(xatol=1e-12))
print("   sup nu = %.6f" % (-r.fun))
save_json("verify_adair", dict(adair_sets=count, adair_rejected=rejected, adair_worst_rel_dev=worst, adair_wellconditioned_sets=count_n, adair_wellconditioned_worst=worst_n, adair_40digit_worst=float(worst_mp), twoform_40digit_worst=float(worst2_mp), twoform_40digit_drift=float(drift2_mp), twoform_worst_rel_dev=worst2, twoform_drift_1000fold=drift2,
                               integer_hill=res3, nu_sup_n2_K1_1e_minus3_K2_1e3=float(-r.fun)))
