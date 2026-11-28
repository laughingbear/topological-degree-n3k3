#!/usr/bin/env python3
"""Claim S, part 1 (exact rational arithmetic): an explicit bound on the region that can contain
positive roots.  Writes nothing; prints a verdict.
Usage: python3 -I S_outer.py <system.json> <T0 as fraction, e.g. 1/2>

Log coordinates: x = exp(y), G(y) = sum_i C_i exp(<a_i - b, y>), b an integer vector (dividing F by
the monomial x^b does not change the positive roots).  We prove G(y) != 0 whenever |y|_inf >= T0.

Write y = t d with |d|_inf = 1, t >= T0.  The set {|d|_inf = 1} is 8 cubes of dimension 3; each is
subdivided into dyadic boxes D.  For a box D let [lo_i, hi_i] be the exact range of <a_i-b, d> on D,
L = max_i lo_i, and S a set of at most 4 indices containing every i with hi_i >= L.  Put
g = L - max_{i not in S} hi_i > 0.  Then for every d in D the largest monomial m_max is attained in S,
and every monomial outside S is <= exp(-g t) m_max.  Choose a row vector lam with lam.C_i = 1 for i in S
(possible: the columns C_S are independent by (U_C)).  Then
   lam.G(y) >= m_max * (1 - exp(-g T0) * sum_{i not in S} |lam.C_i|)
which is > 0 as soon as exp(g T0) > W := sum_{i not in S} |lam.C_i|.  exp(x) is bounded below by a
partial sum of its series (x > 0), so the test is exact rational arithmetic."""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load

syst = load(sys.argv[1])
T0 = Fr(sys.argv[2])
A = syst['A']; N = len(A)
C = [[syst['C'][r][j] for r in range(4)] for j in range(N)]       # columns
b = [round(sum(a[k] for a in A) / N) for k in range(4)]
A0 = [[a[k] - b[k] for k in range(4)] for a in A]


def solve(M, rhs):
    """exact solution of square system M x = rhs (Fractions), None if singular"""
    n = len(M)
    M = [[Fr(v) for v in row] + [Fr(r)] for row, r in zip(M, rhs)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


lam_cache = {}


def lam_for(S):
    """row vector lam with lam.C_i = 1 for i in S (least-norm one if |S| < 4)"""
    if S not in lam_cache:
        gram = [[sum(C[i][k] * C[j][k] for k in range(4)) for j in S] for i in S]
        w = solve(gram, [1] * len(S))
        lam = [sum(w[t] * C[S[t]][k] for t in range(len(S))) for k in range(4)]
        assert all(sum(lam[k] * C[i][k] for k in range(4)) == 1 for i in S)
        W = sum(abs(sum(lam[k] * C[i][k] for k in range(4))) for i in range(N) if i not in S)
        lam_cache[S] = W
    return lam_cache[S]


def exp_lower(x):
    s = Fr(1); term = Fr(1)
    for n in range(1, 40):
        term = term * x / n
        s += term
    return s


def test(face, sgn_, lo3, hi3):
    """box D: coordinate `face` fixed to sgn_, other three in [lo3, hi3] (Fractions)."""
    others = [k for k in range(4) if k != face]
    rng = []
    for a in A0:
        base = a[face] * sgn_
        l = h = Fr(base)
        for t, k in enumerate(others):
            v1 = a[k] * lo3[t]; v2 = a[k] * hi3[t]
            l += min(v1, v2); h += max(v1, v2)
        rng.append((l, h))
    L = max(r[0] for r in rng)
    order = sorted(range(N), key=lambda i: -rng[i][1])
    for k in range(1, 5):
        if rng[order[k]][1] >= L:
            continue                                   # S must contain every possible maximiser
        S = tuple(sorted(order[:k]))
        g = L - rng[order[k]][1]
        W = lam_for(S)
        if W == 0 or exp_lower(g * T0) > W:
            return True
    return False


count = 0; maxdepth = 0
stack = []
for face in range(4):
    for s in (1, -1):
        stack.append((face, s, (Fr(-1),) * 3, (Fr(1),) * 3, 0))
while stack:
    face, s, lo3, hi3, depth = stack.pop()
    if test(face, s, lo3, hi3):
        count += 1; maxdepth = max(maxdepth, depth)
        continue
    if depth > 40:
        print('FAIL: could not exclude direction box', face, s, lo3, hi3); sys.exit(1)
    # bisect all three coordinates
    mid = tuple((l + h) / 2 for l, h in zip(lo3, hi3))
    for m in range(8):
        l2 = tuple(lo3[t] if not (m >> t) & 1 else mid[t] for t in range(3))
        h2 = tuple(mid[t] if not (m >> t) & 1 else hi3[t] for t in range(3))
        stack.append((face, s, l2, h2, depth + 1))
print('monomial shift b =', b)
print('PASS: no positive root with |log x|_inf >= %s  (%d direction boxes, max depth %d, exact arithmetic)' % (T0, count, maxdepth))
