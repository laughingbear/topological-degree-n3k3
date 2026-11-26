"""Exact verification (rational arithmetic) that the fewnomial system below has topological degree +-3
on the positive orthant, via Bihan's Theorem A (arXiv 2611.34336, Thm 2.3 = [Degree1, Thm A]):

  Assume (U_A) and (U_C). For u not in any cone spanned by n-1 columns of C,
     deg F = (-1)^n * sum_{facets G of P, u in cone(C_s, s in G)} sign det C^G * sign det A^{G u {i_G}}
  where A has columns (1, a_i), column i_G (a point off the facet) placed last, columns of G in increasing order.
"""
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

def check(A, C, u):
    N, n = len(A), len(A[0])
    cols = lambda M, idx: [[M[r][j] for j in idx] for r in range(len(M))]
    Ah = [[1] * N] + [[A[i][r] for i in range(N)] for r in range(n)]          # (n+1) x N
    UA = all(det(cols(Ah, s)) != 0 for s in itertools.combinations(range(N), n + 1))
    UC = all(det(cols(C, s)) != 0 for s in itertools.combinations(range(N), n))
    # facets of P: n-subsets with all other points strictly on one side
    facets = []
    for G in itertools.combinations(range(N), n):
        sides = {sgn(det(cols(Ah, list(G) + [i]))) for i in range(N) if i not in G}
        if len(sides) == 1 and 0 not in sides: facets.append(G)
    # u generic: not in any cone spanned by n-1 columns  (n=3: u not in a plane through two columns)
    for s in itertools.combinations(range(N), n - 1):
        if det([[C[r][j] for j in s] + [u[r]] for r in range(n)]) == 0: raise ValueError('u not generic')
    total, used = 0, []
    for G in facets:
        CG = cols(C, G); dCG = det(CG)
        # u in open cone(C_G):  solve CG * lam = u by Cramer, need all lam > 0
        lam = [det([[CG[r][j] if j != c else u[r] for j in range(n)] for r in range(n)]) / dCG for c in range(n)]
        if all(l > 0 for l in lam):
            iG = next(i for i in range(N) if i not in G)
            eps = sgn(det(cols(Ah, list(G) + [iG])))
            total += sgn(dCG) * eps; used.append(G)
    return UA, UC, facets, (-1) ** n * total, used

if __name__ == '__main__':
    S = json.load(open(sys.argv[1]))
    A, C = S['A'], S['C']
    for u in ([1, 2, 3], [-3, 1, 7], [5, -2, -1], [0, 1, -4], [11, 13, -17]):
        try:
            UA, UC, facets, d, used = check(A, C, u)
        except ValueError:
            print('u', u, 'not generic, skipped'); continue
        print(f'u={u}: (U_A)={UA} (U_C)={UC} facets={len(facets)} deg F = {d}  (facets used: {used})')
