"""Exact verification (rational arithmetic), any n: Bihan arXiv:2611.34336 Thm 1.3.
  Assume (U_A),(U_C). For u not in any cone spanned by n-1 columns of C,
     deg F = (-1)^n * sum_{facets G of P, u in cone(C_G)} sign det C^G * sign det Ahat^{G u {i_G}}
  (columns of G increasing, a column i_G off the facet placed last)."""
import itertools, json, sys
from fractions import Fraction as Fr
def det(M):
    M = [[Fr(x) for x in r] for r in M]; n = len(M); d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return Fr(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n): M[r][k] -= f * M[c][k]
    return d
sgn = lambda x: (x > 0) - (x < 0)
cols = lambda M, idx: [[M[r][j] for j in idx] for r in range(len(M))]
def structure(A):
    N, n = len(A), len(A[0])
    Ah = [[1] * N] + [[A[i][r] for i in range(N)] for r in range(n)]
    UA = all(det(cols(Ah, s)) != 0 for s in itertools.combinations(range(N), n + 1))
    facets = []
    for G in itertools.combinations(range(N), n):
        sides = {sgn(det(cols(Ah, list(G) + [i]))) for i in range(N) if i not in G}
        if len(sides) == 1 and 0 not in sides: facets.append((G, sides.pop()))
    return UA, facets
def degree(facets, C, u):
    n = len(C)
    N = len(C[0])
    for s in itertools.combinations(range(N), n - 1):
        if det([[C[r][j] for j in s] + [u[r]] for r in range(n)]) == 0: raise ValueError('u not generic')
    total = 0; used = []
    for G, eps in facets:
        CG = cols(C, G); dCG = det(CG)
        lam = [det([[CG[r][j] if j != c else u[r] for j in range(n)] for r in range(n)]) / dCG for c in range(n)]
        if all(l > 0 for l in lam): total += sgn(dCG) * eps; used.append(G)
    return (-1) ** n * total, used
def UC(C):
    n = len(C); N = len(C[0])
    return all(det(cols(C, s)) != 0 for s in itertools.combinations(range(N), n))
if __name__ == '__main__':
    S = json.load(open(sys.argv[1])); A, C = S['A'], S['C']
    UA, facets = structure(A)
    print('N =', len(A), 'n =', len(A[0]), '(U_A) =', UA, '(U_C) =', UC(C), 'facets =', len(facets))
    import random; random.seed(1); degs = []
    us = [[1, 2, 3, 5][:len(C)], [-3, 1, 7, 2][:len(C)], [5, -2, -1, 9][:len(C)]] + [[random.randint(-50, 50) for _ in C] for _ in range(40)]
    for u in us:
        try: d, used = degree(facets, C, u)
        except ValueError: continue
        degs.append(d)
    print('deg F over', len(degs), 'generic u:', sorted(set(degs)))
