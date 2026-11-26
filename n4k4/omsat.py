"""SAT relaxation: is there a uniform rank-4 chirotope chi on the V vertices of an oriented simplicial 3-sphere K
such that the vertex-map degree is >= d ?   Degree evaluated at u = -C_v for EVERY vertex v:
   deg = sum_{facets G not containing v, -C_v in cone(C_G)} eps_G * chi(G).
Chirotope axioms: 3-term Grassmann-Pluecker relations (sufficient for uniform oriented matroids).
UNSAT  =>  no vertex map (K,C) with C satisfying (U_C) has degree >= d   (realizable chirotopes are a subset)."""
import itertools, json, sys
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType

def parity(t):
    return sum(1 for i in range(len(t)) for j in range(i + 1, len(t)) if t[i] > t[j]) % 2

def build(V, ofac, d, verts=None):
    """ofac: list of (sorted 4-tuple, eps).  returns CNF"""
    pool = IDPool(); cnf = CNF()
    X = {S: pool.id(('x', S)) for S in itertools.combinations(range(V), 4)}
    def chi(t):                      # literal that is true iff det(C_t) > 0, for ordered tuple t of distinct elts
        l = X[tuple(sorted(t))]
        return -l if parity(t) else l
    # GP: for sigma={x,y}, a<b<c<d: [ab][cd] - [ac][bd] + [ad][bc] = 0  -> three signed terms not all equal
    for sig in itertools.combinations(range(V), 2):
        rest = [z for z in range(V) if z not in sig]
        for a, b, c, dd in itertools.combinations(rest, 4):
            terms = [((chi(sig + (a, b)), chi(sig + (c, dd))), 0), ((chi(sig + (a, c)), chi(sig + (b, dd))), 1), ((chi(sig + (a, dd)), chi(sig + (b, c))), 0)]
            # term positive iff (l1 == l2) xor neg
            for want in (True, False):            # forbid all three terms having sign 'want'
                # term has sign want  iff  (l1 xor l2) == (want == neg ? ...)
                opts = []
                for (l1, l2), neg in terms:
                    pos_when_equal = (neg == 0)
                    need_equal = (want == pos_when_equal)
                    # assignments making the term have sign 'want': need_equal -> (l1,l2) in {TT,FF}; else {TF,FT}
                    opts.append([(l1, l2), (-l1, -l2)] if need_equal else [(l1, -l2), (-l1, l2)])
                for combo in itertools.product(*opts):
                    cnf.append([-l for pair in combo for l in pair])
    for v in (range(V) if verts is None else verts):
        lits = []; nneg = 0
        for G, eps in ofac:
            if v in G: continue
            y = pool.id(('y', v, G))
            conds = []
            for i in range(4):
                t = list(G); t[i] = v
                # need chi(t) == -chi(G):  literal e true iff they differ
                a, b = chi(tuple(t)), chi(G)
                e = pool.id(('e', v, G, i))
                cnf.extend([[-e, a, b], [-e, -a, -b], [e, a, -b], [e, -a, b]])
                conds.append(e)
            cnf.extend([[-y, c] for c in conds]); cnf.append([y] + [-c for c in conds])
            s = chi(G) if eps > 0 else -chi(G)        # contribution sign positive
            p = pool.id(('p', v, G)); n = pool.id(('n', v, G))
            cnf.extend([[-p, y], [-p, s], [p, -y, -s], [-n, y], [-n, -s], [n, -y, s]])
            lits += [p, -n]; nneg += 1
        if d + nneg > len(lits): cnf.append([]); continue
        enc = CardEnc.atleast(lits=lits, bound=d + nneg, vpool=pool, encoding=EncType.seqcounter)
        cnf.extend(enc.clauses)
    return cnf, X

if __name__ == '__main__':
    from pysat.solvers import Cadical153
    V = int(sys.argv[1]); d = int(sys.argv[2]); sel = json.loads(sys.argv[3]) if len(sys.argv) > 3 else None
    lines = open(f'pool{V}.txt').read().split('\n')[1:]
    res = {}
    for i, l in enumerate(lines):
        if not l.strip() or (sel is not None and i not in sel): continue
        z = list(map(int, l.split())); sub = list(itertools.combinations(range(9), 4))
        ofac = [(sub[z[1 + 2 * k]], z[2 + 2 * k]) for k in range(z[0])]
        cnf, X = build(V, ofac, d)
        if any(len(c) == 0 for c in cnf.clauses): r = False
        else:
            with Cadical153(bootstrap_with=cnf.clauses) as s: r = s.solve()
        res[i] = r
        print(i, 'SAT' if r else 'UNSAT', flush=True)
    print('SAT:', [i for i, r in res.items() if r])
