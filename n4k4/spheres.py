"""Enumerate combinatorial simplicial 3-spheres with V vertices by bistellar 2-3 / 3-2 flip closure
(up to isomorphism, canonical form via nauty).  Known counts (Altshuler-Steinberg, Lutz): 1,2,5,39,1296 for V=5..9."""
import itertools, sys, json
import pynauty

def canon(F, V):
    F = [tuple(sorted(f)) for f in F]
    n = V + len(F)
    adj = {V + i: list(f) for i, f in enumerate(F)}
    g = pynauty.Graph(n, adjacency_dict=adj, vertex_coloring=[set(range(V)), set(range(V, n))])
    lab = pynauty.canon_label(g)
    pos = {lab[i]: i for i in range(n)}
    return tuple(sorted(tuple(sorted(pos[v] for v in f)) for f in F))

def flips(F):
    F = set(F)
    tri = {}
    edges = {}
    for f in F:
        for t in itertools.combinations(f, 3): tri.setdefault(t, []).append(f)
        for e in itertools.combinations(f, 2): edges.setdefault(e, []).append(f)
    out = []
    for t, fs in tri.items():                     # 2 -> 3
        f1, f2 = fs
        d = next(x for x in f1 if x not in t); e = next(x for x in f2 if x not in t)
        if tuple(sorted((d, e))) in edges: continue
        a, b, c = t
        G = set(F); G.discard(f1); G.discard(f2)
        for p, q in ((a, b), (b, c), (a, c)): G.add(tuple(sorted((p, q, d, e))))
        out.append(G)
    for e, fs in edges.items():                   # 3 -> 2
        if len(fs) != 3: continue
        others = set().union(*fs) - set(e)
        if len(others) != 3: continue
        t = tuple(sorted(others))
        if t in tri: continue
        G = set(F)
        for f in fs: G.discard(f)
        for x in e: G.add(tuple(sorted(t + (x,))))
        out.append(G)
    return out

def cyclic(V):
    """boundary of cyclic 4-polytope C(V,4) via Gale evenness"""
    F = []
    for S in itertools.combinations(range(V), 4):
        Sset = set(S); ok = True
        rest = [x for x in range(V) if x not in Sset]
        for i in rest:
            for j in rest:
                if i < j and sum(1 for s in S if i < s < j) % 2: ok = False
        if ok: F.append(S)
    return F

def enumerate_all(V):
    start = canon(cyclic(V), V)
    seen = {start}; frontier = [start]
    while frontier:
        nxt = []
        for T in frontier:
            for U in flips(T):
                c = canon(U, V)
                if c not in seen: seen.add(c); nxt.append(c)
        frontier = nxt
        print(V, 'so far', len(seen), file=sys.stderr, flush=True)
    return sorted(seen)

if __name__ == '__main__':
    res = {}
    for V in range(5, int(sys.argv[1]) + 1):
        L = enumerate_all(V); res[V] = L
        print(V, len(L), flush=True)
    json.dump({str(k): v for k, v in res.items()}, open('spheres.json', 'w'))
