"""Referee's positive controls (own encoder).  (a) random real C on random Lutz spheres, V = 7,8,9; (b) the exact systems.
Fix chi to the true signs: must be SAT at d = |deg| (after a reflection if deg < 0) and UNSAT at d + 1.
Also a sensitivity test: with the WRONG sign on the middle Grassmann-Pluecker term, true chirotopes must be rejected."""
import sys, json, itertools, random, numpy as np
from pysat.solvers import Solver
from sph import parse_lutz, orient
from enc import Enc
import exact
L = parse_lutz('/home/claude/ext/3_manifolds_all.txt', '/home/claude/ext/3_manifolds_all_type.txt')
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 2026); ntr = int(sys.argv[1]) if __name__ == "__main__" else 0
def fdeg(eps, C, u):                        # C: V x 4 float
    t = 0
    for G, e in eps.items():
        M = C[list(G)].T; lam = np.linalg.solve(M, u)
        if (lam > 0).all(): t += e * np.sign(np.linalg.det(M))
    return int(round(t))
def check(E, eps, chi, d):
    out = []
    for dd in (d, d + 1):
        cl, nv = E.formula(eps, dd)
        if any(len(c) == 0 for c in cl): out.append(False); continue
        with Solver(name='cadical153', bootstrap_with=cl) as s:
            out.append(s.solve(assumptions=[E.x[S] if chi[S] > 0 else -E.x[S] for S in E.x]))
    return out
bad = 0; n = 0; hist = {}; wrong_rej = 0
for V in (7, 8, 9):
    E = Enc(V); Ewrong = Enc(V, gp_middle_sign=+1)
    names = [k for k, (v, F, t) in L.items() if v == V and t == 'S^3']
    for tr in range(ntr):
        nm = names[rng.integers(len(names))]; eps = orient(L[nm][1]); C = rng.normal(size=(V, 4))
        ds = {fdeg(eps, C, rng.normal(size=4)) for _ in range(5)}; assert len(ds) == 1, ds
        d = ds.pop()
        if d < 0: C[:, rng.integers(4)] *= -1; d = -d; assert fdeg(eps, C, rng.normal(size=4)) == d
        chi = {S: np.sign(np.linalg.det(C[list(S)].T)) for S in E.x}
        r = check(E, eps, chi, d); n += 1; hist[(V, d)] = hist.get((V, d), 0) + 1
        if r != [True, False]: bad += 1; print('MISMATCH', V, nm, d, r)
        if tr < 20:
            with Solver(name='cadical153', bootstrap_with=Ewrong.gp) as s:
                wrong_rej += not s.solve(assumptions=[Ewrong.x[S] if chi[S] > 0 else -Ewrong.x[S] for S in Ewrong.x])
print('random controls', n, 'mismatches', bad, 'hist', dict(sorted(hist.items())), '| wrong-GP-sign formula rejected true chirotope in', wrong_rej, 'of', 60)
# exact systems
E = Enc(9); sys6 = json.load(open('/home/claude/author/deg6_systems.json')); fam = json.load(open('/home/claude/author/D44_family.json'))
items = [('deg6:' + k, s['A'], s['C']) for k, s in sys6.items()] + [('family:' + k, fam['A'], C) for k, C in fam['C_by_degree'].items()]
random.seed(3); bad2 = 0
for tag, A, C in items:
    UA, UC, facets, verts = exact.analyse(A, C); assert UA and UC and len(verts) == 9
    u = [random.randint(-9999, 9999) for _ in range(4)]; d = exact.degree(facets, C, u)[0]
    if d < 0: C = [[-x for x in C[0]]] + C[1:]; d2 = exact.degree(facets, C, u)[0]; assert d2 == -d; d = d2
    chi = {S: exact.sg(exact.det([[C[r][i] for i in S] for r in range(4)])) for S in E.x}
    r = check(E, facets, chi, d); bad2 += r != [True, False]
    print(tag, 'facets', len(facets), '|deg| =', d, 'SAT at d, d+1:', r)
print('exact-system controls', len(items), 'mismatches', bad2)
