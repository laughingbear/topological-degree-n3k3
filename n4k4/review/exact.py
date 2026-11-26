"""Referee's own exact verifier (integer arithmetic only, Laplace determinants). Written from Bihan Thm 1.3:
   deg F = (-1)^n sum_{facets G, u in cone(C_G)} sign det C_G * sign det Ahat_{G u {i}}   (G increasing, i off-facet last)."""
import itertools, json, sys, random
from functools import lru_cache

def det(M):
    n = len(M)
    if n == 1: return M[0][0]
    if n == 2: return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    s = 0
    for j in range(n):
        if M[0][j] == 0: continue
        s += (-1) ** j * M[0][j] * det([r[:j] + r[j + 1:] for r in M[1:]])
    return s
sg = lambda x: (x > 0) - (x < 0)

def analyse(A, C):
    """A: list of N points in Z^n (rows). C: n x N integer matrix (rows)."""
    N, n = len(A), len(A[0])
    hat = lambda idx: [[1] * len(idx)] + [[A[i][r] for i in idx] for r in range(n)]   # columns (1,a_i)
    Ccols = lambda idx: [[C[r][i] for i in idx] for r in range(n)]
    UA = all(det(hat(s)) != 0 for s in itertools.combinations(range(N), n + 1))
    UC = all(det(Ccols(s)) != 0 for s in itertools.combinations(range(N), n))
    facets = {}
    for G in itertools.combinations(range(N), n):
        side = [sg(det(hat(list(G) + [i]))) for i in range(N) if i not in G]
        if all(x == 1 for x in side) or all(x == -1 for x in side): facets[G] = side[0]
    verts = sorted(set(v for G in facets for v in G))
    return UA, UC, facets, verts

def is_sphere_like(facets, n):
    """every (n-1)-subset of a facet lies in exactly two facets; coherent orientation check for eps."""
    ridge = {}
    for G, e in facets.items():
        for i in range(n):
            ridge.setdefault(G[:i] + G[i + 1:], []).append((e * (-1) ** i))
    return all(len(v) == 2 for v in ridge.values()), all(len(v) == 2 and v[0] == -v[1] for v in ridge.values())

def degree(facets, C, u):
    n = len(C)
    N = len(C[0])
    col = lambda i: [C[r][i] for r in range(n)]
    # genericity: u not in the span of any n-1 columns (stronger than needed)
    for s in itertools.combinations(range(N), n - 1):
        cols = [col(i) for i in s] + [u]
        if det([[cols[c][r] for c in range(n)] for r in range(n)]) == 0: return None
    tot = 0; contrib = {}
    for G, eps in facets.items():
        cols = [col(i) for i in G]
        d = det([[cols[c][r] for c in range(n)] for r in range(n)])
        inside = True
        for k in range(n):
            cc = cols[:k] + [u] + cols[k + 1:]
            dk = det([[cc[c][r] for c in range(n)] for r in range(n)])
            if sg(dk) != sg(d): inside = False; break          # lambda_k = dk/d must be > 0
        if inside: tot += sg(d) * eps; contrib[G] = sg(d) * eps
    return (-1) ** n * tot, contrib

if __name__ == '__main__':
    fam = json.load(open(sys.argv[1]))
    random.seed(20261126)
    us = [[random.randint(-97, 97) for _ in range(4)] for _ in range(60)]
    A = fam['A']
    for d, C in sorted(fam['C_by_degree'].items()):
        UA, UC, facets, verts = analyse(A, C)
        two, coh = is_sphere_like(facets, 4)
        degs = [degree(facets, C, u) for u in us]; gen = [x[0] for x in degs if x is not None]
        # also at u = -C_v for each vertex
        dv = []
        for v in verts:
            r = degree(facets, C, [-C[k][v] for k in range(4)]); dv.append(None if r is None else r[0])
        fs = [e * sg(det([[C[r][i] for i in G] for r in range(4)])) for G, e in facets.items()]
        print(f'label {d}: U_A={UA} U_C={UC} hull vertices={len(verts)} facets={len(facets)} ridges-in-2={two} eps coherent={coh} '
              f'deg over {len(gen)} generic u: {sorted(set(gen))}; at u=-C_v: {sorted(set(dv), key=str)}; facet signs +{fs.count(1)}/-{fs.count(-1)}')
