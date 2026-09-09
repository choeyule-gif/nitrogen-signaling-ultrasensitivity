"""Methods, 'Deficiency and the ACR criterion', and the buffering example of 'Two kinds of robustness'.

1. n = 1: complexes, linkage classes, rank of the stoichiometric matrix, deficiency (= 1);
   the deficiency of each linkage class (0, 0, 0, so the Deficiency One Algorithm is not
   available); the non-terminal complexes and the count of pairs differing in exactly one
   species (0, so the Shinar-Feinberg criterion does not apply).  General n: |C| = 4n+4,
   l = 3, rank = 3n+1, delta = n, checked for n = 1..6.
2. Complex-mediated bifunctionality (reverse activity on the enzyme-target complex): s_1 is
   the same at every steady state, s_1 = b3(g2+g3)/(g1 g3), for six pairs of totals (ACR).
3. Effector-switched network at n = 1 (parameters A = 1.857, Bhat = 1, Ghat = 0.85, Khat = 0.9,
   T_E = 0.5, l = 1.5): raising T_S sixteenfold moves s_0 from 0.267 to 4.765 and e from
   0.113 to 0.033 -- the ratio, not a concentration, is what is held (manuscript numbers).
"""
import _path  # noqa: F401
import itertools, numpy as np, mpmath as mp
from scipy.optimize import brentq
from esbm import save_json
out = {}

def network(n):
    species = ['E', 'L', 'EL'] + ['S%d' % i for i in range(n+1)] + ['B%d' % i for i in range(n)] + ['G%d' % i for i in range(n)]
    cplx = [('E', 'L'), ('EL',)]; rxn = [(0, 1), (1, 0)]
    for i in range(n):
        k = len(cplx)
        cplx += [('E', 'S%d' % i), ('B%d' % i,), ('E', 'S%d' % (i+1))]
        rxn += [(k, k+1), (k+1, k), (k+1, k+2)]
        k = len(cplx)
        cplx += [('EL', 'S%d' % (i+1)), ('G%d' % i,), ('EL', 'S%d' % i)]
        rxn += [(k, k+1), (k+1, k), (k+1, k+2)]
    # merge identical complexes
    uniq = []; idx = {}
    for c in cplx:
        key = tuple(sorted(c))
        if key not in idx: idx[key] = len(uniq); uniq.append(key)
    rxn = [(idx[tuple(sorted(cplx[a]))], idx[tuple(sorted(cplx[b]))]) for a, b in rxn]
    return species, uniq, rxn
def analyse(n):
    species, cplx, rxn = network(n)
    def vec(c):
        v = np.zeros(len(species))
        for s in c: v[species.index(s)] += 1
        return v
    S = np.array([vec(cplx[b])-vec(cplx[a]) for a, b in rxn]).T
    rank = np.linalg.matrix_rank(S)
    par = list(range(len(cplx)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in rxn:
        ra, rb = find(a), find(b)
        if ra != rb: par[ra] = rb
    classes = {}
    for i in range(len(cplx)): classes.setdefault(find(i), []).append(i)
    lc = len(classes)
    # deficiency of each linkage class
    dcls = []
    for members in classes.values():
        Sm = np.array([vec(cplx[b])-vec(cplx[a]) for a, b in rxn if a in members]).T
        dcls.append(len(members)-1-np.linalg.matrix_rank(Sm))
    # terminal strong linkage classes: complexes from which every reachable complex reaches back
    reach = {i: {i} for i in range(len(cplx))}
    changed = True
    while changed:
        changed = False
        for a, b in rxn:
            new = reach[a] | reach[b]
            if new != reach[a]: reach[a] = new; changed = True
    terminal = [i for i in range(len(cplx)) if all(i in reach[j] for j in reach[i])]
    nonterm = [i for i in range(len(cplx)) if i not in terminal]
    pairs = [(i, j) for i, j in itertools.combinations(nonterm, 2) if np.count_nonzero(vec(cplx[i])-vec(cplx[j])) == 1]
    return dict(complexes=len(cplx), linkage_classes=lc, rank=int(rank), deficiency=len(cplx)-lc-int(rank),
                class_deficiencies=dcls, nonterminal=[cplx[i] for i in nonterm], sf_pairs=len(pairs))
print("1. deficiency")
for n in range(1, 7):
    r = analyse(n)
    print("   n=%d: |C|=%d  l=%d  rank=%d  delta=%d  class deficiencies %s  non-terminal %d  SF pairs %d"
          % (n, r['complexes'], r['linkage_classes'], r['rank'], r['deficiency'], r['class_deficiencies'], len(r['nonterminal']), r['sf_pairs']))
    if n == 1: out['n1'] = r; print("      non-terminal complexes:", ['+'.join(c) for c in r['nonterminal']])
    out.setdefault('general', {})[n] = dict(complexes=r['complexes'], linkage_classes=r['linkage_classes'], rank=r['rank'], deficiency=r['deficiency'])
    assert r['complexes'] == 4*n+4 and r['linkage_classes'] == 3 and r['rank'] == 3*n+1 and r['deficiency'] == n

print("2. complex-mediated bifunctionality: ACR in s_1")
b1, b2, b3, g1, g2, g3 = 1.7, 0.9, 2.3, 3.1, 0.7, 1.9
s1 = b3*(g2+g3)/(g1*g3); be = b1/(b2+b3); ga = g1*be*s1/(g2+g3)
vals = []
for TE, TS in [(0.5, 10.), (0.5, 50.), (3.0, 10.), (3.0, 80.), (0.1, 200.), (20., 5.)]:
    f = lambda s0: s0+s1+(be+2*ga)*s0*TE/(1+(be+ga)*s0)-TS
    s0 = brentq(f, 1e-14, TS, xtol=1e-15, rtol=8.9e-16); e = TE/(1+(be+ga)*s0)
    # verify the full balances of the complex-mediated network at this state
    B = be*e*s0; G = ga*e*s0
    resid = max(abs(b1*e*s0-(b2+b3)*B), abs(g1*B*s1-(g2+g3)*G), abs(b3*B-g3*G), abs(e+B+G-TE), abs(s0+s1+B+2*G-TS))
    vals.append(dict(TE=TE, TS=TS, s0=s0, e=e, s1=s1, residual=resid))
    print("   T_E=%5.1f T_S=%6.1f : s1* = %.10f  (s0* %.4g, e* %.4g; balance residual %.1e)" % (TE, TS, s1, s0, e, resid))
out['acr_s1'] = s1; out['acr_cases'] = vals

print("3. effector-switched network at n = 1: ratio held, concentrations free")
mp.mp.dps = 40
A = mp.mpf('1.857'); Bh = mp.mpf(1); Gh = mp.mpf('0.85'); Kh = mp.mpf('0.9'); TE = mp.mpf('0.5'); l = mp.mpf('1.5'); C = Bh+A*Gh*Kh
def solve(TS):
    F = lambda e: e*(1+A*l)+e*(TS/(1+Kh/l+e*C))*C-TE
    hi = mp.mpf(1)
    while F(hi) < 0: hi *= 10
    e = mp.findroot(F, (mp.mpf(0), hi), solver='bisect', tol=mp.mpf(10)**-35)
    s0 = TS/(1+Kh/l+e*C); return e, s0, s0*Kh/l
TS0 = mp.exp(mp.findroot(lambda t: solve(mp.exp(t))[1]-mp.mpf('0.2669365'), mp.mpf(0)))
e1, s01, s11 = solve(TS0); e16, s016, s116 = solve(16*TS0)
print("   T_S x1 : e=%s s0=%s s1/s0=%s" % (mp.nstr(e1, 6), mp.nstr(s01, 6), mp.nstr(s11/s01, 8)))
print("   T_S x16: e=%s s0=%s s1/s0=%s" % (mp.nstr(e16, 6), mp.nstr(s016, 6), mp.nstr(s116/s016, 8)))
out['buffering'] = dict(TS=float(TS0), e_x1=float(e1), s0_x1=float(s01), e_x16=float(e16), s0_x16=float(s016), ratio_x1=float(s11/s01), ratio_x16=float(s116/s016))
save_json("deficiency_acr", out)
