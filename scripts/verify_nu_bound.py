"""eq (nubound) in high precision (manuscript: 'over random c_i spanning 28 decades no evaluation
of nu exceeded n, the bound was approached to 1e-20 as the intermediate c_i were driven to zero,
and for a binomial P the deviation of nu from 1 was below 6e-43 while any perturbation of a
single coefficient made nu non-constant')."""
import _path  # noqa: F401
import numpy as np
from math import comb
from mpmath import mp
from esbm import save_json
rng = np.random.default_rng(11)
def nu_mp(c, u):
    n = len(c)-1
    w = [c[k]*u**k for k in range(n+1)]; Z = sum(w); p = [x/Z for x in w]
    m = sum(k*p[k] for k in range(n+1)); v = sum((k-m)**2*p[k] for k in range(n+1))
    if m <= 0 or m >= n: return None
    return n*v/(m*(n-m))
out = {}
print("1) nu <= n (50 digits, c_i over 28 decades)")
mp.dps = 50; exceed = 0; worst_gap = None
for n in (2, 3, 4, 6, 12):
    mx = mp.mpf(0); bad = 0
    for _ in range(3000):
        c = [mp.mpf(1)]+[mp.mpf(10)**mp.mpf(float(x)) for x in rng.uniform(-14, 14, n)]
        for lu in rng.uniform(-16, 16, 25):
            r = nu_mp(c, mp.mpf(10)**mp.mpf(float(lu)))
            if r is None: continue
            if r > mx: mx = r
            if r > n+mp.mpf('1e-30'): bad += 1
    exceed += bad
    print("   n=%2d: max nu = %s  exceedances %d" % (n, mp.nstr(mx, 12), bad))
out['exceedances'] = exceed
print("2) equality approached as c_1..c_{n-1} -> 0 (80 digits)")
mp.dps = 80; gaps = {}
for n in (2, 3, 12):
    for ex in (6, 20, 60):
        eps = mp.mpf(10)**(-ex); c = [mp.mpf(1)]+[eps]*(n-1)+[mp.mpf(1)]
        best = mp.mpf(0)
        for lu in np.linspace(-3, 3, 400):
            r = nu_mp(c, mp.mpf(10)**mp.mpf(float(lu)))
            if r and r > best: best = r
        gaps["n%d_eps1e-%d" % (n, ex)] = float(n-best)
        print("   n=%2d eps=1e-%-3d: n - max nu = %s" % (n, ex, mp.nstr(n-best, 6)))
out['gap_to_n'] = gaps
print("3) nu == 1 iff binomial (50 digits)")
mp.dps = 50; dev = {}
for n in (2, 3, 5, 12):
    c = [mp.mpf(comb(n, k))*mp.mpf('2.7')**k for k in range(n+1)]
    d = max(abs(nu_mp(c, mp.mpf(10)**mp.mpf(float(lu)))-1) for lu in np.linspace(-8, 8, 300))
    dev[n] = float(d); print("   n=%2d binomial: max|nu-1| = %s" % (n, mp.nstr(d, 6)))
pert = {}
for p in ('1.0001', '1.01', '1.5'):
    n = 3; c = [mp.mpf(comb(n, k)) for k in range(n+1)]; c[1] *= mp.mpf(p)
    vals = [nu_mp(c, mp.mpf(10)**mp.mpf(float(lu))) for lu in np.linspace(-6, 6, 600)]; vals = [v for v in vals if v]
    pert[p] = [float(min(vals)), float(max(vals))]
    print("   n=3, c1 x %s: nu in [%s, %s] -> not constant" % (p, mp.nstr(min(vals), 10), mp.nstr(max(vals), 10)))
out['binomial_max_dev'] = dev; out['perturbed_range'] = pert; out['binomial_max_dev_overall'] = max(dev.values())
save_json("verify_nu_bound", out)
