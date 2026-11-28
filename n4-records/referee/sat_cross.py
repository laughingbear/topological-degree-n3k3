r"""Cross-check of Claim B with a CDCL solver (CaDiCaL via python-sat) on the CNF produced by
OUR encoder (chiro.gp_terms / unit_conditions).  Each GP relation -> 16 clauses of width 6:
forbid every assignment in which the three signed terms coincide."""
import itertools, sys, os
from pysat.solvers import Cadical153
from common import load, check_closed_manifold, orientation
from chiro import gp_terms, unit_conditions
for path in sys.argv[1:]:
    d = load(path); F = [tuple(sorted(f)) for f in d["facets"]]; check_closed_manifold(F); o = orientation(F)
    V = 1 + max(max(f) for f in F)
    X = {S: i + 1 for i, S in enumerate(itertools.combinations(range(V), 5))}
    for eps in (1, -1):
        s = Cadical153(); n = 0
        for k, v, _ in unit_conditions(F, o, eps):
            s.add_clause([v * X[k]]); n += 1
        for sg in itertools.combinations(range(V), 3):
            rest = [x for x in range(V) if x not in sg]
            for q in itertools.combinations(rest, 4):
                T = gp_terms(sg, q)
                for target in (1, -1):
                    # all three terms == target: choose signs (v1,v2) of each term with c*v1*v2 == target
                    opts = [[(a, b) for a in (1, -1) for b in (1, -1) if c * a * b == target] for c, _, _ in T]
                    for combo in itertools.product(*opts):
                        cl = []
                        for (c, k1, k2), (a, b) in zip(T, combo):
                            cl += [-a * X[k1], -b * X[k2]]
                        s.add_clause(cl); n += 1
        r = s.solve(); s.delete()
        print(f"{os.path.basename(path)} eps={eps:+d}: {n} clauses -> {'SAT' if r else 'UNSAT'}", flush=True)
