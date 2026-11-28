#!/usr/bin/env python3
"""Referee check of Claims P1, P2, P3 and of the Bihan-formula part of S.
Usage: python3 -I check_P.py <V14_polytope_deg16.json> <system_N14_deg16.json> [n_random_u]
Everything is exact integer arithmetic.  The stored 'orientation' is read ONLY in the
final comparison step (clearly marked); no result depends on it."""
import sys, random
from itertools import combinations
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import load, det, sgn, det_cols, lifted

cert = load(sys.argv[1])
syst = load(sys.argv[2])
NU = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
P = cert['polytope_points']
Cv = cert['C_int']
claimed = sorted(tuple(f) for f in cert['facets'])
V = len(P)
ok = True


def report(name, cond, extra=''):
    global ok
    ok &= bool(cond)
    print(('PASS ' if cond else 'FAIL ') + name + (' : ' + extra if extra else ''))


print('== consistency of the two certificate files ==')
report('system.A == polytope_points (both files)', syst['A'] == P and cert['system']['A'] == P)
Ccols = [[syst['C'][r][j] for r in range(4)] for j in range(V)]
report('columns of system.C == C_int (both files)', Ccols == Cv and cert['system']['C'] == syst['C'])
report('14 distinct points', len(set(map(tuple, P))) == 14)

print('== P1: general position, convex position, facets ==')
# all 2002 maximal minors of A-hat
five = {}
zero5 = []
for S in combinations(range(V), 5):
    d = lifted(P, S)
    five[S] = d
    if d == 0:
        zero5.append(S)
report('(U_A) all C(14,5)=2002 5x5 minors of A-hat nonzero', not zero5, 'zero: %s' % zero5[:5] if zero5 else 'min |minor| = %d' % min(abs(d) for d in five.values()))


def side(T, i):
    """sign of det[(1,a_T0),...,(1,a_T3),(1,a_i)] for sorted T, i not in T."""
    S = tuple(sorted(T + (i,)))
    pos = S.index(i)
    return sgn(five[S]) * (-1) ** (4 - pos)

# brute-force hull: a 4-set is a facet iff all other points lie strictly on one side
mine = []
outsign = {}
for T in combinations(range(V), 4):
    ss = {side(T, i) for i in range(V) if i not in T}
    if len(ss) == 1 and 0 not in ss:
        mine.append(T)
        outsign[T] = ss.pop()      # sign det A-hat_{T u {i}} , i last, same for every i off T
report('facets of conv(A) recomputed by brute force == claimed list', mine == claimed,
       'mine %d, claimed %d, mine-claimed %s, claimed-mine %s' % (len(mine), len(claimed), sorted(set(mine) - set(claimed)), sorted(set(claimed) - set(mine))))
report('claimed facets are 56 distinct sorted 4-sets', len(set(claimed)) == 56 and all(list(f) == sorted(f) for f in cert['facets']))
verts = sorted({v for T in mine for v in T})
report('convex position: every point is a vertex of a facet', verts == list(range(V)))
# independent convex-position test: each point strictly separated => (implied by the above in general position).
edges = {e for T in mine for e in combinations(T, 2)}
tris = {}
for T in mine:
    for i in range(4):
        tris.setdefault(T[:i] + T[i + 1:], []).append((T, i))
fvec = (len(verts), len(edges), len(tris), len(mine))
report('f-vector == (14,70,112,56)', fvec == (14, 70, 112, 56), str(fvec))
report('every triangle lies in exactly 2 facets', all(len(v) == 2 for v in tris.values()))
report('Euler: V-E+T-F == 0', fvec[0] - fvec[1] + fvec[2] - fvec[3] == 0)
# links of vertices are 2-spheres (closed surface with chi = 2) -> combinatorial 3-manifold check
linkok = True
for v in range(V):
    L = [tuple(x for x in T if x != v) for T in mine if v in T]
    le = {}
    for t in L:
        for e in combinations(t, 2):
            le[e] = le.get(e, 0) + 1
    lv = {x for t in L for x in t}
    # connectedness
    adj = {x: set() for x in lv}
    for a, b in le:
        adj[a].add(b); adj[b].add(a)
    seen = set(); st = [next(iter(lv))]
    while st:
        x = st.pop()
        if x in seen: continue
        seen.add(x); st.extend(adj[x] - seen)
    linkok &= all(c == 2 for c in le.values()) and len(lv) - len(le) + len(L) == 2 and seen == lv
report('every vertex link is a connected closed surface with chi=2 (a 2-sphere)', linkok)
# dual graph connected
fa = {T: set() for T in mine}
for v_ in tris.values():
    if len(v_) == 2:
        fa[v_[0][0]].add(v_[1][0]); fa[v_[1][0]].add(v_[0][0])
seen = set(); st = [mine[0]]
while st:
    x = st.pop()
    if x in seen: continue
    seen.add(x); st.extend(fa[x] - seen)
report('dual graph connected (orientation unique up to global sign)', len(seen) == len(mine))

print('== orientation, two ways ==')
# (a) combinatorial: propagate o_t (-1)^i = - o_t' (-1)^j across triangles
o_comb = {mine[0]: 1}
st = [mine[0]]
consistent = True
while st:
    T = st.pop()
    for i in range(4):
        tri = T[:i] + T[i + 1:]
        for (T2, j) in tris[tri]:
            if T2 == T: continue
            want = -o_comb[T] * (-1) ** i * (-1) ** j
            if T2 in o_comb:
                consistent &= (o_comb[T2] == want)
            else:
                o_comb[T2] = want; st.append(T2)
report('combinatorial coherent orientation exists (complex orientable)', consistent and len(o_comb) == len(mine))
# (b) geometric: o_t = sign det A-hat_{t u {i}} (i last) -- well defined since all other points on one side
o_geom = dict(outsign)
r = {o_comb[T] * o_geom[T] for T in mine}
report('combinatorial orientation == +-(geometric orientation)', len(r) == 1, 'ratio %s' % r)
# geometric orientation is itself coherent (direct check of the definition)
coh = all(o_geom[a[0]] * (-1) ** a[1] == -o_geom[b[0]] * (-1) ** b[1] for a, b in (tuple(v) for v in tris.values()))
report('geometric orientation satisfies the coherence condition directly', coh)

print('== P2: vertex map, degree, q ==')
dC = {T: det_cols([Cv[v] for v in T]) for T in mine}
report('all 56 facet determinants of C nonzero', all(d != 0 for d in dC.values()), 'min |det| = %d' % min(abs(d) for d in dC.values()))
allmin = {T: det_cols([Cv[v] for v in T]) for T in combinations(range(V), 4)}
zeros4 = [T for T, d in allmin.items() if d == 0]
report('(U_C) all C(14,4)=1001 4x4 minors of C nonzero', not zeros4, ('zero: %s' % zeros4[:5]) if zeros4 else 'min |minor| = %d' % min(abs(d) for d in allmin.values()))


def generic(u):
    return all(det_cols([Cv[a], Cv[b], Cv[c], u]) != 0 for a, b, c in combinations(range(V), 3))


def in_open_cone(T, u):
    D = dC[T]
    cols = [Cv[v] for v in T]
    for k in range(4):
        c2 = cols[:k] + [u] + cols[k + 1:]
        if sgn(det_cols(c2)) != sgn(D):
            return False
    return True


def degrees(u):
    cov = [T for T in mine if in_open_cone(T, u)]
    d_comb = sum(o_comb[T] * sgn(dC[T]) for T in cov)
    d_geom = sum(o_geom[T] * sgn(dC[T]) for T in cov)
    bihan = (-1) ** 4 * sum(sgn(dC[T]) * o_geom[T] for T in cov)
    return d_comb, d_geom, bihan, len(cov)

for nm, o in (('combinatorial', o_comb), ('geometric', o_geom)):
    s = [o[T] * sgn(dC[T]) for T in mine]
    print('   facet signs with %s orientation: +%d / -%d' % (nm, s.count(1), s.count(-1)))
sg = [o_geom[T] * sgn(dC[T]) for T in mine]
report('all 56 facets have the same sign s_t', len(set(sg)) == 1)

rng = random.Random(20261127)
us = []
# structured u: +-e_i perturbed, +- each C_v perturbed, barycentres of facets perturbed, then random at 3 scales
for i in range(4):
    for s in (1, -1):
        u = [rng.randint(-3, 3) for _ in range(4)]; u[i] += s * 1000; us.append(u)
for v in range(V):
    for s in (1, -1):
        us.append([s * 1000 * Cv[v][k] + rng.randint(-7, 7) for k in range(4)])
for T in mine:
    us.append([1000 * sum(Cv[v][k] for v in T) + rng.randint(-7, 7) for k in range(4)])
for T in mine:   # near a wall: 3 vectors of a facet, tiny push to either side
    for s in (1, -1):
        us.append([10 ** 6 * sum(Cv[v][k] for v in T[:3]) + s * Cv[T[3]][k] + rng.randint(-1, 1) for k in range(4)])
for scale in (5, 100, 10 ** 6):
    for _ in range(NU):
        us.append([rng.randint(-scale, scale) for _ in range(4)])
res = {}
ngen = 0
covers = {}
for u in us:
    if not generic(u):
        continue
    ngen += 1
    r_ = degrees(u)
    res[r_[:3]] = res.get(r_[:3], 0) + 1
    covers[r_[3]] = covers.get(r_[3], 0) + 1
print('   generic u tested: %d of %d ; (deg_comb, deg_geom, bihan) -> count: %s' % (ngen, len(us), res))
print('   number of facet cones containing u -> count: %s' % covers)
report('degree independent of u, same for both orientations up to global sign', len(res) == 1)
(dc, dg, bh), = res.keys() if len(res) == 1 else [(None, None, None)]
report('|deg| == 16', dc is not None and abs(dc) == 16 and abs(dg) == 16, 'deg_comb=%s deg_geom=%s' % (dc, dg))
report("Bihan's formula value == -16 at every generic u tested", bh == -16, str(bh))
q = sum(1 for T in mine if o_geom[T] * sgn(dC[T]) != sgn(dg)) if dg else None
report('q == 0', q == 0, 'q=%s' % q)

print('== P3: Theorem M arithmetic ==')
from fractions import Fraction as Fr
F_, E_ = fvec[3], fvec[1]
rhs = 2 * F_ - Fr(E_, 2) - Fr(q or 0, 2)
report('4|deg| < 2F - E/2 - q/2', dg is not None and 4 * abs(dg) < rhs, '%d < %s' % (4 * abs(dg or 0), rhs))
print('   largest |deg| allowed by Theorem M with q=0 for this f-vector: %d' % max(d for d in range(100) if 4 * d < 2 * F_ - Fr(E_, 2)))

print('== comparison with the authors\' stored orientation (not used above) ==')
o_theirs = {tuple(f): s for f, s in zip(cert['facets'], cert['orientation'])}
r = {o_theirs[T] * o_comb[T] for T in mine} if set(o_theirs) == set(mine) else None
report("authors' stored orientation == +-(mine)", r is not None and len(r) == 1, 'ratio %s' % r)
r2 = {o_theirs[T] * o_geom[T] for T in mine} if r else None
print('   stored orientation relative to geometric (sign det A-hat, extra index last): %s' % r2)
report('stored degree/q/f fields', cert['degree'] == 16 and cert['q'] == 0 and cert['f'] == [14, 70, 112, 56])
print('ALL PASS' if ok else 'SOME CHECKS FAILED')
