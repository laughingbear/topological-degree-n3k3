"""Controls at higher degree: hill-climb real C on random spheres (own code), then the fixed-chirotope SAT test."""
import sys, numpy as np
from controls import *      # reuses L, fdeg, check (re-runs the base controls with 0 trials)
rng = np.random.default_rng(99); hist = {}; bad = 0; n = 0
for V, runs in ((7, 25), (8, 40), (9, 90)):
    E = Enc(V); names = [k for k, (v, F, t) in L.items() if v == V and t == 'S^3']
    for run in range(runs):
        nm = names[rng.integers(len(names))]; eps = orient(L[nm][1]); C = rng.normal(size=(V, 4)); u = rng.normal(size=4)
        cur = abs(fdeg(eps, C, u))
        for step in range(400):
            v = rng.integers(V); old = C[v].copy()
            C[v] = rng.normal(size=4) if rng.random() < .3 else C[v] + rng.choice([.5, .15, .04]) * rng.normal(size=4)
            e = abs(fdeg(eps, C, u))
            if e >= cur: cur = e
            else: C[v] = old
        ds = {fdeg(eps, C, rng.normal(size=4)) for _ in range(6)}
        if len(ds) != 1: continue
        d = ds.pop()
        if d < 0: C[:, 0] *= -1; d = -d
        chi = {S: np.sign(np.linalg.det(C[list(S)].T)) for S in E.x}
        r = check(E, eps, chi, d); n += 1; hist[(V, d)] = hist.get((V, d), 0) + 1
        if r != [True, False]: bad += 1; print('MISMATCH', V, nm, d, r)
print('hill-climbed controls', n, 'mismatches', bad, 'hist', dict(sorted(hist.items())))
