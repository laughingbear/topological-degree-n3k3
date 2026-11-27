"""Same formula as omsat.build, with the (sphere-independent) Grassmann-Pluecker clauses cached per V.
check_same() asserts that the clause lists are identical to omsat.build's."""
import itertools
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from omsat import parity
_cache = {}

def _gp(V):
    if V in _cache: return _cache[V]
    X = {S: i + 1 for i, S in enumerate(itertools.combinations(range(V), 4))}
    def chi(t):
        l = X[tuple(sorted(t))]
        return -l if parity(t) else l
    cl = []
    for sig in itertools.combinations(range(V), 2):
        rest = [z for z in range(V) if z not in sig]
        for a, b, c, dd in itertools.combinations(rest, 4):
            terms = [((chi(sig + (a, b)), chi(sig + (c, dd))), 0), ((chi(sig + (a, c)), chi(sig + (b, dd))), 1), ((chi(sig + (a, dd)), chi(sig + (b, c))), 0)]
            for want in (True, False):
                opts = []
                for (l1, l2), neg in terms:
                    need_equal = (want == (neg == 0))
                    opts.append([(l1, l2), (-l1, -l2)] if need_equal else [(l1, -l2), (-l1, l2)])
                for combo in itertools.product(*opts):
                    cl.append([-l for pair in combo for l in pair])
    _cache[V] = (X, cl)
    return _cache[V]

def build_fast(V, ofac, d):
    """returns (clauses, X); clauses is a list (the cached GP prefix is shared, do not mutate)"""
    X, gp = _gp(V)
    pool = IDPool(start_from=len(X) + 1); cl = []
    def chi(t):
        l = X[tuple(sorted(t))]
        return -l if parity(t) else l
    for v in range(V):
        lits = []; nneg = 0
        for G, eps in ofac:
            if v in G: continue
            y = pool.id(('y', v, G)); conds = []
            for i in range(4):
                t = list(G); t[i] = v
                a, b = chi(tuple(t)), chi(G)
                e = pool.id(('e', v, G, i))
                cl.extend([[-e, a, b], [-e, -a, -b], [e, a, -b], [e, -a, b]]); conds.append(e)
            cl.extend([[-y, c] for c in conds]); cl.append([y] + [-c for c in conds])
            s = chi(G) if eps > 0 else -chi(G)
            p = pool.id(('p', v, G)); n = pool.id(('n', v, G))
            cl.extend([[-p, y], [-p, s], [p, -y, -s], [-n, y], [-n, -s], [n, -y, s]])
            lits += [p, -n]; nneg += 1
        if d + nneg > len(lits): return None, X          # empty clause: trivially UNSAT
        cl.extend(CardEnc.atleast(lits=lits, bound=d + nneg, vpool=pool, encoding=EncType.seqcounter).clauses)
    return gp, cl, X

def check_same(V, ofac, d):
    from omsat import build
    cnf, X0 = build(V, ofac, d); r = build_fast(V, ofac, d)
    if r[0] is None: return any(len(c) == 0 for c in cnf.clauses)
    gp, cl, X = r
    return X == X0 and cnf.clauses == gp + cl
