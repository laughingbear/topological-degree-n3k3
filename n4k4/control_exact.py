"""Positive controls near the bound (cold-review correction 2): for exact integer systems (A, C) -- the nine degree-6 systems
and the family of degrees 0..6 -- fix the chirotope variables to the true signs (exact integers).
Expected: UNSAT at d+1 for both reflections; at d = |deg| SAT for exactly one reflection (both if d = 0)."""
import json, itertools
from pysat.solvers import Cadical153
from omsat import build
from verify_exact import structure, degree, det, cols
def _deg(fac, C, u):
    try: return degree(fac, C, u)[0]
    except ValueError: return None
systems = [(f'deg6_systems[{k}]', s['A'], s['C']) for k, s in json.load(open('data/deg6_systems.json')).items()]
F = json.load(open('data/D44_family.json')); systems += [(f'family[{k}]', F['A'], C) for k, C in F['C_by_degree'].items()]
bad = 0; hist = {}
for name, A, C in systems:
    UA, fac = structure(A); d = abs(next(x for x in map(lambda u: _deg(fac, C, u), ([1, 2, 3, 5], [-3, 1, 7, 2], [5, -2, -1, 9], [7, -11, 13, -3], [2, 9, -4, -6])) if x is not None))
    sign = {T: det(cols(C, T)) > 0 for T in itertools.combinations(range(9), 4)}
    res = {}
    for flip in (False, True):
        for dd in (d, d + 1):
            cnf, X = build(9, fac, dd)
            with Cadical153(bootstrap_with=cnf.clauses) as sv:
                res[(flip, dd)] = sv.solve(assumptions=[X[T] if (sign[T] != flip) else -X[T] for T in X])
    at_d = (res[(False, d)] and res[(True, d)]) if d == 0 else (res[(False, d)] != res[(True, d)])
    ok = at_d and not res[(False, d + 1)] and not res[(True, d + 1)]
    hist[d] = hist.get(d, 0) + 1; bad += not ok
    print(name, '|deg| =', d, 'OK' if ok else 'MISMATCH')
print('controls', len(systems), 'mismatches', bad, 'degree histogram', dict(sorted(hist.items())))
