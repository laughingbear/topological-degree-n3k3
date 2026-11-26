"""Cold-review verifier, written from scratch (does not import the authors' code).

Inputs: system.json with A (N x n exponent vectors) and C (n x N coefficient matrix).
Checks, in exact rational arithmetic (sympy):
  * (U_A): every (n+1)x(n+1) minor of Ahat (columns (1,a_i)) is nonzero
  * (U_C): every n x n minor of C is nonzero
  * facets of P = conv(A) (via sign-uniformity of Ahat minors), vertex set
  * deg F by Bihan's Theorem A (arXiv:2611.34336, Thm 2.3):
        deg F = (-1)^n sum_{facets G, u in cone(C_G)} sign det C^G * sign det Ahat^{G (increasing) + [i_G] last}
    for many u (random integer vectors), checking u genericity exactly.
  * an independent float computation: signed spherical area of the PL map on the
    outward-oriented boundary complex (orientation computed by outward normal, not by A-minor signs).
"""
import sys, json, itertools, random
import numpy as np
from sympy import Matrix, sign as ssign

def main(path, nu=25, seed=7):
    S = json.load(open(path))
    A = [list(map(int, a)) for a in S['A']]
    C = Matrix(S['C'])
    N, n = len(A), len(A[0])
    assert C.shape == (n, N)
    Ah = Matrix([[1]*N] + [[A[i][r] for i in range(N)] for r in range(n)])
    UA = all(Ah[:, list(s)].det() != 0 for s in itertools.combinations(range(N), n+1))
    UC = all(C[:, list(s)].det() != 0 for s in itertools.combinations(range(N), n))
    facets = []
    for G in itertools.combinations(range(N), n):
        sg = {ssign(Ah[:, list(G)+[i]].det()) for i in range(N) if i not in G}
        if 0 not in sg and len(sg) == 1:
            facets.append(G)
    verts = sorted(set(itertools.chain(*facets)))
    print(f'N={N} n={n}  (U_A)={UA}  (U_C)={UC}  #facets={len(facets)}  vertices={verts}')
    # Euler / simplicial check for n=3
    if n == 3:
        V = len(verts); E = len({tuple(sorted(e)) for G in facets for e in itertools.combinations(G, 2)})
        print(f'  V={V} E={E} F={len(facets)}  V-E+F={V-E+len(facets)}  (F==2V-4: {len(facets)==2*V-4})')
    rnd = random.Random(seed)
    degs = []
    tries = 0
    while len(degs) < nu and tries < 10*nu:
        tries += 1
        u = Matrix([rnd.randint(-50, 50) for _ in range(n)])
        # genericity: u not in span of any n-1 columns (n=3: det(C_i, C_j, u) != 0)
        if any(Matrix.hstack(C[:, list(s)], u).det() == 0 for s in itertools.combinations(range(N), n-1)):
            continue
        tot = 0
        for G in facets:
            CG = C[:, list(G)]
            lam = CG.LUsolve(u)
            if all(l > 0 for l in lam):
                iG = next(i for i in range(N) if i not in G)
                tot += ssign(CG.det()) * ssign(Ah[:, list(G)+[iG]].det())
        degs.append(int((-1)**n * tot))
    print(f'  Theorem A degree over {len(degs)} generic u: {sorted(set(degs))}')
    # independent float method: outward orientation + signed spherical areas
    if n == 3:
        P = np.array(A, float); cen = P[verts].mean(0)
        Cf = np.array(C.tolist(), float).T
        tot = 0.0; signs = []
        for (i, j, k) in facets:
            nr = np.cross(P[j]-P[i], P[k]-P[i])
            if nr @ (P[i]-cen) < 0: j, k = k, j
            a, b, c = (Cf[x]/np.linalg.norm(Cf[x]) for x in (i, j, k))
            num = a @ np.cross(b, c); den = 1 + a@b + b@c + c@a
            tot += 2*np.arctan2(num, den); signs.append(int(np.sign(num)))
        print(f'  signed-area degree = {tot/(4*np.pi):.12f}; triangle orientation signs: {signs}')
    return UA, UC, facets, degs

if __name__ == '__main__':
    main(sys.argv[1])
