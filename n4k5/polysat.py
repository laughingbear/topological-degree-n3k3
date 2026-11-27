"""Necessary condition for polytopality of an oriented simplicial 3-sphere K on V vertices (Bokowski-Sturmfels style):
a rank-5 uniform chirotope chi on the V points (signs of det[1,a_i]) with  chi(F, w) = eps_F  for every facet F and w not in F,
satisfying the 3-term Grassmann-Pluecker relations.  UNSAT => K is not the boundary of a 4-polytope in general position
(hence not polytopal: simplicial polytopes can be perturbed to general position)."""
import itertools, json, sys
from pysat.solvers import Cadical153
def parity(t): return sum(1 for i in range(len(t)) for j in range(i + 1, len(t)) if t[i] > t[j]) % 2
_gp = {}
def gp5(V):
    if V in _gp: return _gp[V]
    X = {S: i + 1 for i, S in enumerate(itertools.combinations(range(V), 5))}
    def chi(t):
        l = X[tuple(sorted(t))]; return -l if parity(t) else l
    cl = []
    for sig in itertools.combinations(range(V), 3):
        rest = [z for z in range(V) if z not in sig]
        for a, b, c, d in itertools.combinations(rest, 4):
            terms = [((chi(sig + (a, b)), chi(sig + (c, d))), 0), ((chi(sig + (a, c)), chi(sig + (b, d))), 1), ((chi(sig + (a, d)), chi(sig + (b, c))), 0)]
            for want in (True, False):
                opts = []
                for (l1, l2), neg in terms:
                    eq = (want == (neg == 0)); opts.append([(l1, l2), (-l1, -l2)] if eq else [(l1, -l2), (-l1, l2)])
                for combo in itertools.product(*opts): cl.append([-l for pair in combo for l in pair])
    _gp[V] = (X, cl, chi); return _gp[V]
def polytopal_sat(V, ofac):
    """ofac: list of (sorted facet, eps) with eps a coherent orientation. returns model dict or None"""
    X, cl, chi = gp5(V); units = []
    for G, eps in ofac:
        for w in range(V):
            if w in G: continue
            l = chi(tuple(G) + (w,)); units.append([l if eps > 0 else -l])
    with Cadical153(bootstrap_with=cl + units) as s:
        if not s.solve(): return None
        m = set(s.get_model()); return {S: (X[S] in m) for S in X}
