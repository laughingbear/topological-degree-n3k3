"""Independent referee code for the topdeg n=4 claims.  Written from the definitions in CLAIMS.md sec. 0 only.

Conventions (my own):
  * A facet is a sorted 4-tuple.  A coherent orientation is a map facet -> +-1 such that the two facets
    containing each triangle induce opposite orientations on it (boundary rule: deleting position k gives
    the face with sign (-1)^k).
  * Exact arithmetic: vertex vectors may be integers (python int) -> determinants exact; signs of Cramer
    coefficients exact.
"""
import itertools, math
from fractions import Fraction
import numpy as np


# ---------------------------------------------------------------- combinatorics
def coherent_orientation(facets):
    facets = [tuple(sorted(f)) for f in facets]
    tri = {}
    for i, f in enumerate(facets):
        for k in range(4):
            tri.setdefault(f[:k] + f[k + 1:], []).append((i, k))
    for t, lst in tri.items():
        if len(lst) != 2:
            raise ValueError(f"triangle {t} in {len(lst)} facets: not a pseudomanifold")
    o = [0] * len(facets)
    o[0] = 1
    queue = [0]
    while queue:
        i = queue.pop()
        f = facets[i]
        for k in range(4):
            (a, ka), (b, kb) = tri[f[:k] + f[k + 1:]]
            if a == i:
                j, kj, ki = b, kb, ka
            else:
                j, kj, ki = a, ka, kb
            # induced: o_i (-1)^ki  must equal  - o_j (-1)^kj
            want = -o[i] * (-1) ** ki * (-1) ** kj
            if o[j] == 0:
                o[j] = want
                queue.append(j)
            elif o[j] != want:
                raise ValueError("non-orientable")
    if 0 in o:
        raise ValueError("disconnected")
    return facets, o


def edges_of(facets):
    return sorted({e for f in facets for e in itertools.combinations(f, 2)})


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


# ---------------------------------------------------------------- exact linear algebra
def det(M):
    """exact determinant (Bareiss) of a square list-of-lists of ints/Fractions; floats also fine."""
    M = [list(r) for r in M]
    n = len(M)
    if any(isinstance(x, float) for r in M for x in r):
        return float(np.linalg.det(np.array(M, float)))
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j])
                M[i][j] = M[i][j] // prev if not isinstance(M[i][j], Fraction) else M[i][j] / prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def sgn(x):
    return (x > 0) - (x < 0)


def cols(C, idx):
    """matrix whose columns are C[i] for i in idx (C: list of 4-vectors) -> list of rows"""
    return [[C[i][r] for i in idx] for r in range(4)]


def colmat(vecs):
    return [[v[r] for v in vecs] for r in range(len(vecs[0]))]


# ---------------------------------------------------------------- the invariants
class VertexMap:
    def __init__(self, facets, C):
        """C: list (length V) of 4-vectors (ints preferred)."""
        self.facets, self.o = coherent_orientation(facets)
        self.C = [list(c) for c in C]
        self.V = len(C)
        self.E_list = edges_of(self.facets)
        self.E = len(self.E_list)
        self.F = len(self.facets)
        self.d = [det(cols(self.C, f)) for f in self.facets]
        if any(x == 0 for x in self.d):
            raise ValueError("not a vertex map (singular facet)")
        self.s = [self.o[i] * sgn(self.d[i]) for i in range(self.F)]

    def flip(self):
        self.o = [-x for x in self.o]
        self.s = [-x for x in self.s]

    # --- degree by cone counting at u (exact if u integer)
    def cone_count(self, u):
        tot = 0
        for i, f in enumerate(self.facets):
            D = self.d[i]
            ok = True
            for k in range(4):
                vecs = [self.C[j] for j in f]
                vecs[k] = u
                Dk = det(colmat(vecs))
                if Dk == 0:
                    raise ValueError("u on a triangle cone")
                if sgn(Dk) != sgn(D):
                    ok = False
                    break
            if ok:
                tot += self.s[i]
        return tot

    def degree(self, rng, trials=6, R=10**6):
        vals = set()
        for _ in range(trials):
            u = [int(x) for x in rng.integers(-R, R, 4)]
            try:
                vals.add(self.cone_count(u))
            except ValueError:
                continue
        if len(vals) != 1:
            raise RuntimeError(f"degree not constant: {vals}")
        return vals.pop()

    # --- link degree w_v: genuinely a 3-dim computation in C_v^perp
    def link_degree(self, v, rng, trials=4, R=10**6):
        Cv = self.C[v]
        nv = sum(x * x for x in Cv)

        def P(x):  # nv * orthogonal projection to Cv^perp (integer if inputs integer)
            dot = sum(a * b for a, b in zip(x, Cv))
            return [nv * a - dot * b for a, b in zip(x, Cv)]

        tris = []
        for i, f in enumerate(self.facets):
            if v in f:
                rest = [x for x in f if x != v]
                # orientation of link triangle (rest sorted) induced by o_i via t = [v, rest]
                ol = self.o[i] * perm_sign([v] + rest)
                tris.append((rest, ol))
        vals = set()
        for _ in range(trials):
            w = P([int(x) for x in rng.integers(-R, R, 4)])  # random direction in Cv^perp
            tot = 0
            bad = False
            for rest, ol in tris:
                Pv = [P(self.C[j]) for j in rest]
                D = det(colmat([Cv] + Pv))  # orientation of Cv^perp: (Cv, basis) positive
                inside = True
                for k in range(3):
                    vecs = list(Pv)
                    vecs[k] = w
                    Dk = det(colmat([Cv] + vecs))
                    if Dk == 0:
                        bad = True
                    if sgn(Dk) != sgn(D):
                        inside = False
                if inside:
                    tot += ol * sgn(D)
            if not bad:
                vals.add(tot)
        if len(vals) != 1:
            raise RuntimeError(f"link degree not constant {vals}")
        return vals.pop()

    # --- edge winding w_e: angle sum of the projected link polygon, ordered by the orientation
    def edge_winding(self, e):
        a, b = e
        Ca = np.array(self.C[a], float)
        Cb = np.array(self.C[b], float)
        Q, _ = np.linalg.qr(np.column_stack([Ca, Cb, np.eye(4)]))
        e1, e2 = Q[:, 2], Q[:, 3]
        if np.linalg.det(np.column_stack([Ca, Cb, e1, e2])) < 0:
            e2 = -e2
        # oriented link cycle: facet t = [a,b,x,y] with orientation -> directed segment x->y
        succ = {}
        segs = []
        for i, f in enumerate(self.facets):
            if a in f and b in f:
                x, y = [z for z in f if z not in e]
                ol = self.o[i] * perm_sign([a, b, x, y])
                if ol < 0:
                    x, y = y, x
                succ[x] = y
                segs.append((x, y, self.s[i] * ol))
        # walk the cycle
        start = next(iter(succ))
        cyc = [start]
        while succ[cyc[-1]] != start:
            cyc.append(succ[cyc[-1]])
        assert len(cyc) == len(succ)
        pts = {z: (np.dot(self.C[z], e1), np.dot(self.C[z], e2)) for z in cyc}
        tot = 0.0
        for k in range(len(cyc)):
            p, q = pts[cyc[k]], pts[cyc[(k + 1) % len(cyc)]]
            tot += math.atan2(p[0] * q[1] - p[1] * q[0], p[0] * q[0] + p[1] * q[1])
        w = tot / (2 * math.pi)
        wr = round(w)
        assert abs(w - wr) < 1e-7, w
        pe = sum(1 for i, f in enumerate(self.facets) if a in f and b in f and self.s[i] > 0)
        me = len(cyc)
        return wr, pe, me


# ---------------------------------------------------------------- measures
def solid_angle_3(v1, v2, v3):
    """normalised solid angle (fraction of S^2) of the cone spanned by three vectors in R^3."""
    a, b, c = (np.asarray(x, float) / np.linalg.norm(x) for x in (v1, v2, v3))
    num = abs(np.dot(a, np.cross(b, c)))
    den = 1 + np.dot(a, b) + np.dot(b, c) + np.dot(c, a)
    om = 2 * math.atan2(num, den)
    if om < 0:
        om += 4 * math.pi
    return om / (4 * math.pi)


def orth_basis_perp(Cv):
    Q, _ = np.linalg.qr(np.column_stack([np.asarray(Cv, float), np.eye(4)]))
    return Q[:, 1:4]


def orthant_probs_4(M, tol=1e-9):
    """M: 4x4 basis (columns).  u ~ N(0,I), c = M^{-1}u ~ N(0,S).  Returns dict sign-pattern -> probability
    for all 16 orthants, via scipy's Genz algorithm."""
    from scipy.stats import multivariate_normal
    Mi = np.linalg.inv(np.asarray(M, float))
    S = Mi @ Mi.T
    out = {}
    for pat in itertools.product([1, -1], repeat=4):
        Dg = np.diag(pat)
        S2 = Dg @ S @ Dg  # P(pat_k c_k > 0 all k) = P(-pat c < 0) = cdf at 0
        mv = multivariate_normal(mean=np.zeros(4), cov=S2, allow_singular=False)
        out[pat] = float(mv.cdf(np.zeros(4), abseps=tol, releps=tol, maxpts=4_000_000))
    return out
