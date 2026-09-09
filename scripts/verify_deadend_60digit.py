"""eq (nogo) in the raw species variables with a dead-end complex Z = E.S_0 whose association
constant sweeps eleven decades (1e-6..1e5), n = 3, solved by 60-digit Newton iteration on the
six raw unknowns (e, l, s_0..s_3) with the chain relation NOT imposed.
Manuscript: the chain holds to 1.3e-29; the ratios themselves drift by 5.4e-5 at
T_E/T_L = 1.1e-4 and by 5.4e-10 at 1.1e-9, in proportion to T_E/T_L (eq. sandwich)."""
import _path  # noqa: F401
import mpmath as mp
from esbm import save_json
mp.mp.dps = 60
n = 3; A = mp.mpf(1); Bh = mp.mpf('0.4'); Gh = mp.mpf('0.4'); Kh = mp.mpf(1)
rho = A*Gh*Kh/Bh
g3 = [mp.mpf(1)/rho]*n; b3 = [mp.mpf(1)]*n
Bhi = [(n-i)*Bh for i in range(n)]; Ghi = [(i+1)*Gh for i in range(n)]
TS = mp.mpf(1)
def eqs_factory(KZ, TE, TL):
    def eqs(*v):
        e, l = mp.exp(v[0]), mp.exp(v[1]); s = [mp.exp(x) for x in v[2:]]
        eps = A*e*l
        b = [Bhi[i]*e*s[i] for i in range(n)]; g = [Ghi[i]*eps*s[i+1] for i in range(n)]
        z = KZ*e*s[0]
        out = [b3[i]*b[i]-g3[i]*g[i] for i in range(n)]
        out += [e+eps+sum(b)+sum(g)+z-TE, l+eps+sum(g)-TL, sum(s)+sum(b)+sum(g)+z-TS]
        return out
    return eqs
Ki = [(n-i+1)/mp.mpf(i)*Kh for i in range(1, n+1)]
res = {}
for TEs in ('1e-4', '1e-9'):
    TE = mp.mpf(TEs); TL = mp.mpf('0.9')+TE
    v = [mp.log(TE/2), mp.log(mp.mpf('0.9'))]+[mp.log(mp.mpf('0.5')/(2**i)) for i in range(n+1)]
    out = []
    for lk in [x/4 for x in range(-24, 21)]:
        v = mp.findroot(eqs_factory(mp.mpf(10)**lk, TE, TL), tuple(v), tol=mp.mpf(10)**-45)
        v = [v[i] for i in range(n+3)]
        e, l = mp.exp(v[0]), mp.exp(v[1]); s = [mp.exp(v[i+2]) for i in range(n+1)]
        rat = [s[i+1]/s[i] for i in range(n)]
        dev = max(abs(rat[i]-Ki[i]/l)/(Ki[i]/l) for i in range(n))
        out.append((lk, l, rat, dev))
    l0 = out[0][1]; r0 = out[0][2]
    dl = max(abs(o[1]/l0-1) for o in out); dr = max(max(abs(o[2][i]/r0[i]-1) for i in range(n)) for o in out); dv = max(o[3] for o in out)
    print("T_E=%s  T_E/T_L=%s" % (TEs, mp.nstr(TE/TL, 3)))
    print("   chain s_i/s_(i-1) = K_i/l holds to      %s" % mp.nstr(dv, 3))
    print("   free effector l drifts over 11 decades: %s" % mp.nstr(dl, 3))
    print("   ratios themselves drift:                %s" % mp.nstr(dr, 3))
    res[TEs] = dict(TE_over_TL=float(TE/TL), chain_rel_dev=float(dv), l_drift=float(dl), ratio_drift=float(dr))
save_json("verify_deadend_60digit", res)
