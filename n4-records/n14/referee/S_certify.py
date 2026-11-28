#!/usr/bin/env python3
"""Claim S, part 3: certification of each positive root in EXACT rational interval arithmetic.
Usage: python3 -I S_certify.py <system.json> <output of S_inner>

For every root box reported by S_inner we polish the centre with floating-point Newton (a guess only),
round it to a rational point xh, and run the Krawczyk test for the ORIGINAL polynomial map
     F(x) = sum_i C_i x^{a_i}        (integer exponents, x in R^4_{>0})
on the box X = prod_k [xh_k (1-delta), xh_k (1+delta)] using Fractions only:
     K = xh - Y F(xh) + (I - Y J(X)) (X - xh),   K inside int(X)  =>  exactly one root of F in X.
J(X) is an interval enclosure of the Jacobian (monomials are monotone on the positive orthant, so the
enclosure of x^e on X is [prod lo^e, prod hi^e], exact).  The sign of det J on X is obtained from an
interval Laplace expansion; if it excludes 0 the Jacobian sign at the root is certified."""
import sys, os, math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load

syst = load(sys.argv[1])
A = syst['A']; N = len(A)
C = [[syst['C'][r][j] for r in range(4)] for j in range(N)]
b = [round(sum(a[k] for a in A) / N) for k in range(4)]
A0 = [[a[k] - b[k] for k in range(4)] for a in A]
boxes = []
for line in open(sys.argv[2]):
    if line.startswith('ROOT'):
        t = line.split()
        v = list(map(float, t[2:]))
        boxes.append((int(t[1]), [(v[2 * k], v[2 * k + 1]) for k in range(4)]))


def fnewton(y):
    import numpy as np
    A0n = np.array(A0, float); Cn = np.array(C, float).T
    y = np.array(y, float)
    for _ in range(60):
        e = A0n @ y; M = e.max(); m = np.exp(e - M)
        g = Cn @ m; J = (Cn * m) @ A0n
        d = np.linalg.solve(J, -g)
        y = y + d
        if abs(d).max() < 1e-15:
            break
    return [float(v) for v in y]


class IV:
    __slots__ = ('lo', 'hi')

    def __init__(s, lo, hi=None):
        s.lo = Fr(lo); s.hi = Fr(lo if hi is None else hi)

    def __add__(s, o): return IV(s.lo + o.lo, s.hi + o.hi)
    def __sub__(s, o): return IV(s.lo - o.hi, s.hi - o.lo)

    def __mul__(s, o):
        p = (s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi)
        return IV(min(p), max(p))


def idet(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    s = IV(0)
    for k in range(n):
        minor = [row[1:] for i, row in enumerate(M) if i != k]
        t = M[k][0] * idet(minor)
        s = s + t if k % 2 == 0 else s - t
    return s


def certify(yc, delta=Fr(1, 2 ** 34), bits=44):
    xh = [Fr(round(math.exp(v) * 2 ** bits), 2 ** bits) for v in yc]
    lo = [x * (1 - delta) for x in xh]; hi = [x * (1 + delta) for x in xh]

    def mono(pt, e):
        r = Fr(1)
        for k in range(4):
            r *= pt[k] ** e[k]
        return r
    Fx = [sum(C[i][r] * mono(xh, A[i]) for i in range(N)) for r in range(4)]
    JX = [[IV(0) for _ in range(4)] for _ in range(4)]
    Jh = [[Fr(0)] * 4 for _ in range(4)]
    for i in range(N):
        for j in range(4):
            if A[i][j] == 0:
                continue
            e = list(A[i]); e[j] -= 1
            iv = IV(mono(lo, e), mono(hi, e)); mh = mono(xh, e)
            for r in range(4):
                JX[r][j] = JX[r][j] + IV(C[i][r] * A[i][j]) * iv
                Jh[r][j] += C[i][r] * A[i][j] * mh
    # approximate inverse of Jh: scale rows/cols in floats, then convert back to Fractions (a choice only)
    import numpy as np
    sc = max(abs(v) for row in Jh for v in row)
    Y = np.linalg.inv(np.array([[float(v / sc) for v in row] for row in Jh]))
    Y = [[Fr(float(Y[a][c])) / sc for c in range(4)] for a in range(4)]
    inside = True
    for a in range(4):
        K = IV(xh[a] - sum(Y[a][c] * Fx[c] for c in range(4)))
        for c in range(4):
            yj = IV(0)
            for t in range(4):
                yj = yj + IV(Y[a][t]) * JX[t][c]
            K = K + (IV(1 if a == c else 0) - yj) * IV(lo[c] - xh[c], hi[c] - xh[c])
        inside &= (K.lo > lo[a] and K.hi < hi[a])
    d = idet(JX)
    sign = 1 if d.lo > 0 else (-1 if d.hi < 0 else 0)
    return inside, sign, lo, hi


ok = True; tot = 0; cert = []
print('%d root boxes read' % len(boxes))
for s_cpp, bx in boxes:
    yc = fnewton([(l + h) / 2 for l, h in bx])
    inbox = all(l < v < h for v, (l, h) in zip(yc, bx))
    inside, sign, lo, hi = certify(yc)
    # the certified x-box must lie inside the S_inner y-box (compare via exact bounds on exp):
    # exp(l) < lo  and hi < exp(h) are checked with floats widened by a relative 1e-12 margin.
    contained = all(math.exp(l) * (1 + 1e-12) < float(lo[k]) and float(hi[k]) < math.exp(h) * (1 - 1e-12) for k, (l, h) in enumerate(bx))
    good = inside and sign != 0 and sign == s_cpp and inbox and contained
    ok &= good; tot += sign
    cert.append((lo, hi))
    print(('PASS' if good else 'FAIL') + '  x = (%s)  sign det Jac F = %+d  [Krawczyk %s, in S_inner box %s]' % (', '.join('%.9f' % float((l + h) / 2) for l, h in zip(lo, hi)), sign, inside, contained))
# pairwise disjointness of the certified boxes (exact)
dis = all(any(cert[i][1][k] < cert[j][0][k] or cert[j][1][k] < cert[i][0][k] for k in range(4)) for i in range(len(cert)) for j in range(i))
ok &= dis
npos = sum(1 for s, _ in boxes if s > 0)
print('certified boxes pairwise disjoint (exact): %s' % dis)
print('certified nondegenerate positive roots: %d ; Jacobian signs +%d / -%d ; signed sum %+d' % (len(boxes), npos, len(boxes) - npos, tot))
print('ALL PASS' if ok else 'SOME FAILED')
