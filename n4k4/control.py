"""Positive control for the SAT encoding: for real C (random, and the degree-6 optima), fix the chirotope variables
to the true signs; the formula with bound d = |true degree| must be SAT, and with d+1 UNSAT."""
import itertools, json, sys, numpy as np, subprocess
from pysat.solvers import Cadical153
from omsat import build
sub = list(itertools.combinations(range(9), 4))
def load(V):
    out = []
    for l in open(f'pool{V}.txt').read().split('\n')[1:]:
        if l.strip():
            z = list(map(int, l.split())); out.append([(sub[z[1 + 2 * k]], z[2 + 2 * k]) for k in range(z[0])])
    return out
def truedeg(ofac, C, u):
    t = 0
    for G, eps in ofac:
        M = C[list(G)].T; lam = np.linalg.solve(M, u)
        if (lam > 0).all(): t += eps * np.sign(np.linalg.det(M))
    return int(round(t))
rng = np.random.default_rng(5); bad = 0; n = 0; hist = {}
for V in (7, 8, 9):
    P = load(V)
    for trial in range(int(sys.argv[1])):
        i = rng.integers(len(P)); ofac = P[i]; C = rng.normal(size=(V, 4))
        d = truedeg(ofac, C, rng.normal(size=4))
        if d < 0: C[:, 0] *= -1; d = -d            # reflect -> degree >= 0
        assert truedeg(ofac, C, rng.normal(size=4)) == d
        ass_sign = {S: np.linalg.det(C[list(S)].T) > 0 for S in itertools.combinations(range(V), 4)}
        ok = []
        for dd in (d, d + 1):
            cnf, X = build(V, ofac, dd)
            if any(len(c) == 0 for c in cnf.clauses): ok.append(False); continue
            with Cadical153(bootstrap_with=cnf.clauses) as s:
                ok.append(s.solve(assumptions=[X[S] if ass_sign[S] else -X[S] for S in X]))
        n += 1; hist[(V, d)] = hist.get((V, d), 0) + 1
        if ok != [True, False]: bad += 1; print('MISMATCH', V, i, d, ok)
print('controls', n, 'mismatches', bad, 'degree histogram', dict(sorted(hist.items())))
