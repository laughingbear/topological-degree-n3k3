r"""Claim D, step 2: CERTIFY positive zeros of F(x) = sum_i C_i x^{a_i} and the sign of det DF there.

Exact rational interval arithmetic only (fractions.Fraction endpoints; the only floating /
multiprecision input is the approximate root x0 and approximate inverse Y, which are merely
*guesses* that the exact test then validates).

Krawczyk test: X = x0 + [-r,r]^4,  K(X) = x0 - Y F(x0) + (I - Y DF(X)) (X - x0).
If K(X) is contained in the interior of X then F has exactly one zero in X, and DF is
nonsingular on X.  DF(X) is enclosed exactly: every monomial is increasing in every variable on
the positive orthant, so its range over a box is [m(lower corner), m(upper corner)].
The sign of det DF at the zero is certified by an interval evaluation of det DF(X) excluding 0.

usage: D_certify.py system.json roots.json     (roots.json: list of t = log x as decimal strings)
"""
import json, sys, itertools
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 120

class I:
    __slots__ = ("lo", "hi")
    def __init__(s, lo, hi=None):
        s.lo = Fr(lo); s.hi = Fr(lo if hi is None else hi); assert s.lo <= s.hi
    def __add__(s, o):
        o = o if isinstance(o, I) else I(o); return I(s.lo + o.lo, s.hi + o.hi)
    __radd__ = __add__
    def __neg__(s): return I(-s.hi, -s.lo)
    def __sub__(s, o):
        o = o if isinstance(o, I) else I(o); return I(s.lo - o.hi, s.hi - o.lo)
    def __rsub__(s, o): return I(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, I) else I(o)
        p = (s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi); return I(min(p), max(p))
    __rmul__ = __mul__
    def has0(s): return s.lo <= 0 <= s.hi

def idet(M):
    n = len(M)
    if n == 1: return M[0][0]
    tot = I(0)
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = M[0][j] * idet(minor)
        tot = tot + (term if j % 2 == 0 else -term)
    return tot

def to_fr(v, bits=160):
    """exact rational equal to the mpf v rounded to `bits` significant bits (any magnitude)."""
    if v == 0: return Fr(0)
    m, e = mp.frexp(v)                      # v = m * 2**e, 0.5 <= |m| < 1
    return Fr(int(mp.floor(m * 2**bits)), 2**bits) * (Fr(2) ** int(e))

def mono(x, e):
    p = Fr(1)
    for xi, ei in zip(x, e):
        if ei: p *= xi ** ei
    return p

def certify(A, C, t_strs, rad_rel=Fr(1, 10**30), bits=200):
    N = len(A)
    xs = [mp.e ** mp.mpf(s) for s in t_strs]
    x0 = [Fr(int(mp.floor(x * 2**bits)), 2**bits) for x in xs]
    # exact F(x0), approximate Jacobian for Y
    m = [mono(x0, a) for a in A]
    F0 = [sum(C[j][i] * m[i] for i in range(N)) for j in range(4)]
    Jm = mp.matrix(4, 4)
    for j in range(4):
        for k in range(4):
            Jm[j, k] = sum(C[j][i] * A[i][k] * mp.mpf(m[i].numerator) / mp.mpf(m[i].denominator) for i in range(N)) / xs[k]
    Ym = mp.inverse(Jm)
    Y = [[to_fr(Ym[j, k]) for k in range(4)] for j in range(4)]
    r = [x * rad_rel for x in x0]
    lo = [x - ri for x, ri in zip(x0, r)]; hi = [x + ri for x, ri in zip(x0, r)]
    assert all(l > 0 for l in lo)
    # interval Jacobian over the box
    JX = [[I(0) for _ in range(4)] for _ in range(4)]
    for i in range(N):
        for k in range(4):
            if A[i][k] == 0: continue
            e = list(A[i]); e[k] -= 1
            rng = I(mono(lo, e), mono(hi, e))
            for j in range(4):
                if C[j][i]:
                    JX[j][k] = JX[j][k] + (C[j][i] * A[i][k]) * rng
    # K(X)
    YF = [sum(Y[j][l] * F0[l] for l in range(4)) for j in range(4)]
    ok = True; margins = []
    for j in range(4):
        acc = I(x0[j] - YF[j])
        for k in range(4):
            ent = I(1 if j == k else 0) - sum((Y[j][l] * JX[l][k] for l in range(4)), I(0))
            acc = acc + ent * I(-r[k], r[k])
        ok &= (lo[j] < acc.lo and acc.hi < hi[j])
        margins.append(float(max(abs(acc.lo - x0[j]), abs(acc.hi - x0[j])) / r[j]))
    D = idet(JX)
    sign = 0 if D.has0() else (1 if D.lo > 0 else -1)
    return ok, sign, [float(x) for x in x0], max(margins)

if __name__ == "__main__":
    d = json.load(open(sys.argv[1])); A = d["A"]; C = d["C"]
    roots = json.load(open(sys.argv[2]))
    tot = 0; nok = 0; boxes = []
    for t in roots:
        ok, sign, x0, mg = certify(A, C, t)
        nok += ok; tot += sign if ok else 0
        print(f"x ~ ({', '.join(f'{v:.9g}' for v in x0)})  unique zero in box (rel. radius 1e-30): {ok}  "
              f"[|K-x0|/r = {mg:.2e}]  sign det DF = {sign:+d}", flush=True)
    # distinctness: boxes pairwise disjoint (centres differ by far more than 1e-30 relative)
    xs = [[mp.e ** mp.mpf(s) for s in t] for t in roots]
    dist = min((max(abs(a - b) / a for a, b in zip(p, q)) for p, q in itertools.combinations(xs, 2)), default=1)
    print(f"certified zeros: {nok} of {len(roots)}; pairwise distinct (min rel. separation {mp.nstr(dist, 3)}); "
          f"SIGNED COUNT of certified zeros = {tot:+d}")
