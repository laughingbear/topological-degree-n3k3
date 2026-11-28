#!/usr/bin/env python3
"""Referee checks for E1 (symmetry) and E4 (smaller numbers).
Usage: python3 -I E_checks.py <V14_polytope_deg16.json> [out_small.json]
The integer parts are exact.  The part about the exactly symmetric (unrounded) configuration uses
mpmath at 60 digits with the base points taken from 'source' as printed; it is NOT an exact check
(the true symmetric configuration has algebraic coordinates) and is reported as numerical evidence."""
import sys, os, json, random
from itertools import combinations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load, det, sgn, det_cols, lifted
import mpmath as mp

cert = load(sys.argv[1])
P = cert['polytope_points']; Cv = cert['C_int']; V = 14
claimed = sorted(tuple(f) for f in cert['facets'])
src = cert['source']


def analyse(P, Cv, nu=40, seed=7):
    """exact: returns dict(UA, UC, facets, degs) for integer points P and integer vectors Cv."""
    V = len(P)
    five = {S: lifted(P, S) for S in combinations(range(V), 5)}
    UA = all(d != 0 for d in five.values())
    facets = {}
    for T in combinations(range(V), 4):
        ss = set()
        for i in range(V):
            if i in T: continue
            S = tuple(sorted(T + (i,))); pos = S.index(i)
            ss.add(sgn(five[S]) * (-1) ** (4 - pos))
        if len(ss) == 1 and 0 not in ss:
            facets[T] = ss.pop()
    minors = {T: det_cols([Cv[v] for v in T]) for T in combinations(range(V), 4)}
    UC = all(d != 0 for d in minors.values())
    res = dict(UA=UA, UC=UC, facets=facets, nvert=len({v for T in facets for v in T}), degs=None,
               facet_dets_ok=all(minors[T] != 0 for T in facets))
    if not res['facet_dets_ok']:
        return res
    rng = random.Random(seed); degs = set(); n = 0
    while n < nu:
        u = [rng.randint(-10 ** 4, 10 ** 4) for _ in range(4)]
        if any(det_cols([Cv[a], Cv[b], Cv[c], u]) == 0 for a, b, c in combinations(range(V), 3)):
            continue
        n += 1; tot = 0
        for T, eps in facets.items():
            D = minors[T]; cols = [Cv[v] for v in T]
            if all(sgn(det_cols(cols[:k] + [u] + cols[k + 1:])) == sgn(D) for k in range(4)):
                tot += sgn(D) * eps
        degs.add(tot)
    res['degs'] = degs
    res['q'] = min(sum(1 for T, e in facets.items() if sgn(minors[T]) * e > 0), sum(1 for T, e in facets.items() if sgn(minors[T]) * e < 0))
    return res


print('== E1: symmetry of the rounded (claimed) object ==')
rot = lambda v: (v + 1) % 7 if v < 7 else 7 + (v - 7 + 1) % 7
img = sorted(tuple(sorted(rot(v) for v in T)) for T in claimed)
print('facet list invariant under the rotation i -> i+1 in each block: %s' % (img == claimed))
orbits = set()
for T in claimed:
    o = []; t = T
    for _ in range(7):
        o.append(t); t = tuple(sorted(rot(v) for v in t))
    orbits.add(min(o))
print('number of Z_7-orbits of facets: %d (each of size 7: %s)' % (len(orbits), 7 * len(orbits) == len(claimed)))
refl = lambda v: (-v) % 7 if v < 7 else 7 + (-(v - 7)) % 7
img2 = sorted(tuple(sorted(refl(v) for v in T)) for T in claimed)
print('facet list invariant under the reflection i -> -i (dihedral symmetry): %s' % (img2 == claimed))

print('== E1: exactly symmetric configuration (60-digit floating point, NOT exact) ==')
mp.mp.dps = 60


def orbit(base, p, q):
    out = []
    for b_ in base:
        b_ = [mp.mpf(repr(x)) for x in b_]
        for i in range(7):
            t1 = 2 * mp.pi * p * i / 7; t2 = 2 * mp.pi * q * i / 7
            out.append([mp.cos(t1) * b_[0] - mp.sin(t1) * b_[1], mp.sin(t1) * b_[0] + mp.cos(t1) * b_[1],
                        mp.cos(t2) * b_[2] - mp.sin(t2) * b_[3], mp.sin(t2) * b_[2] + mp.cos(t2) * b_[3]])
    return out


Xs = orbit(src['X'], src['p1'], src['q1']); Ys = orbit(src['Y'], src['p2'], src['q2'])
mdet = lambda cols: mp.det(mp.matrix([[cols[j][i] for j in range(len(cols))] for i in range(len(cols))]))
five = {S: mdet([[mp.mpf(1)] + Xs[i] for i in S]) for S in combinations(range(V), 5)}
print('symmetric points: smallest |5x5 minor| = %s (largest %s)' % (mp.nstr(min(abs(d) for d in five.values()), 5), mp.nstr(max(abs(d) for d in five.values()), 5)))
fs = {}
for T in combinations(range(V), 4):
    ss = set()
    for i in range(V):
        if i in T: continue
        S = tuple(sorted(T + (i,))); pos = S.index(i)
        ss.add(int(mp.sign(five[S])) * (-1) ** (4 - pos))
    if len(ss) == 1:
        fs[T] = ss.pop()
print('symmetric configuration: %d facets; same facet list as the rounded one: %s' % (len(fs), sorted(fs) == claimed))
m4 = {T: mdet([Ys[v] for v in T]) for T in combinations(range(V), 4)}
small = sorted((abs(d), T) for T, d in m4.items())
scale4 = max(abs(d) for d in m4.values())
nz = [T for a, T in small if a < scale4 * mp.mpf(10) ** -40]
print('symmetric vectors: %d of 1001 4x4 minors vanish to 40 digits (relative); smallest nonvanishing/largest = %s' % (len(nz), mp.nstr(min(a for a, T in small if T not in nz) / scale4, 5)))
print('   vanishing minors that are facets: %d' % sum(1 for T in nz if T in fs))
if not any(T in fs for T in nz):
    rng = random.Random(3); degs = set()
    for _ in range(40):
        u = [mp.mpf(rng.randint(-1000, 1000)) for _ in range(4)]
        tot = 0
        for T, eps in fs.items():
            cols = [Ys[v] for v in T]; D = m4[T]
            if all(mp.sign(mdet(cols[:k] + [u] + cols[k + 1:])) == mp.sign(D) for k in range(4)):
                tot += int(mp.sign(D)) * eps
        degs.add(tot)
    print('symmetric configuration: degree (geometric orientation, 40 random u) = %s; facet signs all equal: %s' % (sorted(degs), len({int(mp.sign(m4[T])) * e for T, e in fs.items()}) == 1))

print('== E4: smaller numbers (exact verification of every candidate) ==')
Xf = [[float(v) for v in r] for r in Xs]; Yf = [[float(v) for v in r] for r in Ys]
mx = max(abs(v) for r in Xf for v in r); my = max(abs(v) for r in Yf for v in r)


def rounded(M, mxv, s):
    return [[round(v / mxv * s) for v in r] for r in M]


def good(Ai, Ci):
    if len(set(map(tuple, Ai))) < 14: return None
    r = analyse(Ai, Ci)
    if r['UA'] and r['UC'] and r['nvert'] == 14 and r['degs'] is not None and len(r['degs']) == 1 and abs(next(iter(r['degs']))) == 16:
        return r
    return None


Cfull = cert['C_int']
bestA = None
for sa in range(4, 59):
    Ai = rounded(Xf, mx, sa)
    r = good(Ai, Cfull)
    if r:
        bestA = (sa, Ai, r); break
print('smallest exponent scale (with the original C) giving (U_A),(U_C), 14 vertices, |deg|=16: %s' % (bestA[0] if bestA else None))
bestC = None
if bestA:
    for sc in range(2, 21):
        Ci = rounded(Yf, my, sc)
        r = good(bestA[1], Ci)
        if r:
            bestC = (sc, Ci, r); break
    print('smallest coefficient scale with that A: %s' % (bestC[0] if bestC else None))
if bestA and bestC:
    Ai = bestA[1]; mn = [min(a[k] for a in Ai) for k in range(4)]
    Ai = [[a[k] - mn[k] for k in range(4)] for a in Ai]; Ci = bestC[1]
    r = analyse(Ai, Ci, nu=300, seed=11)
    print('SMALL SYSTEM: exponents in [0,%d], coefficients in [%d,%d]' % (max(max(a) for a in Ai), min(min(c) for c in Ci), max(max(c) for c in Ci)))
    print('  A =', Ai); print('  C (one vector per monomial) =', Ci)
    print('  exact: (U_A)=%s (U_C)=%s vertices=%d facets=%d Bihan-sign degree over 300 generic u = %s (value of formula = geometric-orientation sum, n=4) q=%d same facet list as claimed: %s'
          % (r['UA'], r['UC'], r['nvert'], len(r['facets']), sorted(r['degs']), r['q'], sorted(r['facets']) == claimed))
    if len(sys.argv) > 2:
        json.dump(dict(N=14, n=4, k=9, A=Ai, C=[[Ci[j][r_] for j in range(14)] for r_ in range(4)],
                       note='referee-found smaller version; verified exactly by E_checks.py (formula value %s)' % sorted(r['degs'])), open(sys.argv[2], 'w'))
        print('  written to', sys.argv[2])
