"""4D gluing check (vertex-map level, exact rationals): connected sum of two copies of the degree-6 vertex map on 9 vertices
along a positive facet of copy 1 and a negative facet of copy 2  ->  14 vertices, degree should be 12 in absolute value."""
import json, itertools, random
from fractions import Fraction as Fr
from verify_exact import structure, det, cols, sgn
S = json.load(open('data/system_deg6.json')); A, C = S['A'], S['C']
UA, fac = structure(A)                       # (sorted G, eps)
col = lambda C, j: [C[r][j] for r in range(4)]
sig = {G: e * sgn(det(cols(C, G))) for G, e in fac}
def par(t): return (-1) ** sum(1 for i in range(4) for j in range(i + 1, 4) if t[i] > t[j])
def solve(M, b):                             # M 4x4 (list of rows) , exact
    d = det(M); return [det([[M[r][j] if j != c else b[r] for j in range(4)] for r in range(4)]) / d for c in range(4)]
def degree(ofac, Cc, u):
    t = 0
    for G, e in ofac:
        M = [[Cc[j][r] for j in G] for r in range(4)]
        lam = solve(M, u)
        assert all(l != 0 for l in lam)
        if all(l > 0 for l in lam): t += e * sgn(det(M))
    return t
def glue(fac1, C1, fac2, C2, t1, t2):
    """fac: list (sorted tuple, eps); C: dict vertex -> vector. t1 positive in 1, t2 negative in 2. returns glued (fac, C)."""
    e1 = dict(fac1)[t1]; e2 = dict(fac2)[t2]
    phi = dict(zip(t2, t1)) if e1 * e2 == -1 else dict(zip(t2, (t1[1], t1[0], t1[2], t1[3])))
    # M with M*C2[x] = C1[phi(x)]
    B2 = [[C2[x][r] for x in t2] for r in range(4)]              # columns C2[x]
    Mrows = []
    for r in range(4):                                           # row r of M: M = B1 * B2^{-1}  -> solve B2^T m = B1row
        b = [C1[phi[x]][r] for x in t2]
        Mrows.append(solve([[B2[c][k] for c in range(4)] for k in range(4)], b))
    assert det(Mrows) > 0
    n1 = max(C1) + 1; new = {}
    for v in C2:
        if v in phi: new[v] = phi[v]
        else: new[v] = n1; n1 += 1
    Cg = dict(C1)
    for v in C2:
        if v not in phi: Cg[new[v]] = [sum(Mrows[r][k] * C2[v][k] for k in range(4)) for r in range(4)]
    out = [(G, e) for G, e in fac1 if G != t1]
    for G, e in fac2:
        if G == t2: continue
        H = tuple(new[v] for v in G); out.append((tuple(sorted(H)), e * par(H)))
    return out, Cg
C1 = {v: [Fr(x) for x in col(C, v)] for v in range(9)}
s = sum(sig.values()); print('copy: facets', len(fac), 'sum of signs', s)
pos = [G for G in sig if sig[G] > 0]; neg = [G for G in sig if sig[G] < 0]
random.seed(1); us = [[random.randint(-40, 40) for _ in range(4)] for _ in range(6)]
print('deg copy', {degree(fac, C1, u) for u in us})
cur_f, cur_C = fac, C1
for step in range(1, 4):
    sg = {G: e * sgn(det([[cur_C[j][r] for j in G] for r in range(4)])) for G, e in cur_f}
    t1 = next(G for G in sg if sg[G] > 0)
    cur_f, cur_C = glue(cur_f, cur_C, fac, C1, t1, neg[0])
    V = len(cur_C)
    uc = all(det([[cur_C[j][r] for j in T] for r in range(4)]) != 0 for T in itertools.combinations(range(V), 4))
    print('glued', step + 1, 'copies: V =', V, 'facets', len(cur_f), 'all 4-subsets independent:', uc, 'deg', {degree(cur_f, cur_C, u) for u in us})
