"""For N = 8 (n = 3, k = 4): search for |deg| = 4 by enumerating ALL labelled triangulations of S^2 on 8 vertices
(via edge-flip closure) and, for random vertex images C, maximising the degree over all of them.
Degree of the vertex map for oriented triangulation T: sum of signed solid angles / 4pi."""
import numpy as np, itertools, sys, pickle, os

V = 8
def neighbors_ok(tris):
    return True

def canon(tris):
    return frozenset(frozenset(t) for t in tris)

def flips(T):
    """all triangulations reachable by one edge flip. T: frozenset of frozenset triangles."""
    edges = {}
    for t in T:
        for e in itertools.combinations(sorted(t), 2):
            edges.setdefault(e, []).append(t)
    E = set(edges)
    out = []
    for e, (t1, t2) in edges.items():
        a, b = e
        c = next(iter(t1 - {a, b})); d = next(iter(t2 - {a, b}))
        if c == d or tuple(sorted((c, d))) in E: continue
        newT = set(T); newT.discard(t1); newT.discard(t2)
        newT.add(frozenset((a, c, d))); newT.add(frozenset((b, c, d)))
        out.append(frozenset(newT))
    return out

def all_triangulations():
    if os.path.exists('tri8.pkl'): return pickle.load(open('tri8.pkl', 'rb'))
    # start: bipyramid over hexagon 0..5 with apexes 6,7
    start = frozenset(frozenset((i, (i + 1) % 6, a)) for i in range(6) for a in (6, 7))
    seen = {start}; frontier = [start]
    while frontier:
        nxt = []
        for T in frontier:
            for U in flips(T):
                if U not in seen: seen.add(U); nxt.append(U)
        frontier = nxt
        print('triangulations so far', len(seen), flush=True)
    L = list(seen); pickle.dump(L, open('tri8.pkl', 'wb')); return L

def orient(T):
    """consistently orient the triangles of a sphere triangulation (returns list of ordered triples)."""
    T = [tuple(sorted(t)) for t in T]
    oriented = {T[0]: T[0]}; stack = [T[0]]
    edge2tri = {}
    for t in T:
        for e in itertools.combinations(t, 2): edge2tri.setdefault(frozenset(e), []).append(t)
    while stack:
        t = stack.pop(); o = oriented[t]
        for k in range(3):
            a, b = o[k], o[(k + 1) % 3]
            for u in edge2tri[frozenset((a, b))]:
                if u == t or u in oriented: continue
                c = next(x for x in u if x not in (a, b))
                oriented[u] = (b, a, c); stack.append(u)     # neighbour traverses edge as b->a
    return [oriented[t] for t in T]

if __name__ == '__main__':
    Ts = all_triangulations()
    print('labelled triangulations:', len(Ts), flush=True)
    O = np.array([orient(T) for T in Ts])            # (M, 12, 3)
    # index triples -> position in table of signed solid angles for ordered (i<j<k) with sign of permutation
    trip = list(itertools.combinations(range(V), 3)); tid = {t: i for i, t in enumerate(trip)}
    idx = np.zeros(O.shape[:2], int); sgn = np.zeros(O.shape[:2])
    for m in range(O.shape[0]):
        for f in range(12):
            a, b, c = O[m, f]; s = tuple(sorted((a, b, c))); idx[m, f] = tid[s]
            # parity of permutation (a,b,c) -> sorted
            perm = [s.index(x) for x in (a, b, c)]
            inv = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
            sgn[m, f] = -1 if inv % 2 else 1
    rng = np.random.default_rng(int(sys.argv[1])); trials = int(sys.argv[2]); best = 0
    def solid(a, b, c):
        a, b, c = a / np.linalg.norm(a), b / np.linalg.norm(b), c / np.linalg.norm(c)
        return 2 * np.arctan2(a @ np.cross(b, c), 1 + a @ b + b @ c + c @ a)
    from collections import Counter
    hist = Counter()
    for t in range(trials):
        C = rng.normal(size=(V, 3))
        w = np.array([solid(C[i], C[j], C[k]) for (i, j, k) in trip])
        deg = (w[idx] * sgn).sum(1) / (4 * np.pi)          # degree for every labelled triangulation
        d = np.abs(np.round(deg)).astype(int); m = d.max(); hist[m] += 1
        if m > best or m >= 4:
            best = max(best, m); j = int(d.argmax())
            print(f'trial {t}: max |deg| = {m}', flush=True)
            if m >= 4:
                np.save('n8_deg4.npy', {'C': C, 'T': O[j]}); print('FOUND degree 4'); break
    print('best', best, 'histogram of per-trial max', dict(hist))
