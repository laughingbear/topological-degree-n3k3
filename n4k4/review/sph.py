"""Referee's sphere utilities: parse Lutz's list, nauty-free canonical form, orientation, sanity checks."""
import re, itertools, json
from collections import deque

def parse_lutz(path, typepath):
    txt = open(path).read(); types = dict(re.findall(r'(manifold_3_\d+_\d+):\s*(\S.*)', open(typepath).read()))
    out = {}
    for name, body in re.findall(r'(manifold_3_\d+_\d+)=\[((?:[\[\]\d,\s])*?)\]\s*(?=manifold|\Z)', txt):
        F = [tuple(int(x) - 1 for x in f.split(',')) for f in re.findall(r'\[([\d,\s]+?)\]', body)]
        V = int(name.split('_')[2]); assert sorted(set(v for f in F for v in f)) == list(range(V)), name
        out[name] = (V, sorted(tuple(sorted(f)) for f in F), types[name].strip())
    return out

def canon(F):
    """Canonical form of a strongly connected 3-pseudomanifold without nauty: every (facet, vertex order) seed determines a
    labelling by breadth-first gluing across triangles; take the lexicographically smallest relabelled facet list."""
    F = [tuple(sorted(f)) for f in F]; tri = {}
    for f in F:
        for t in itertools.combinations(f, 3): tri.setdefault(t, []).append(f)
    assert all(len(v) == 2 for v in tri.values())
    nV = len(set(v for f in F for v in f)); best = None
    for f0 in F:
        for order in itertools.permutations(f0):
            lab = {v: i for i, v in enumerate(order)}; q = deque([f0]); seen = {f0}
            while q:
                f = q.popleft()
                for t in sorted(itertools.combinations(sorted(f, key=lab.get), 3), key=lambda t: [lab[x] for x in t]):
                    a, b = tri[tuple(sorted(t))]; g = b if a == f else a
                    if g in seen: continue
                    seen.add(g); w = next(x for x in g if x not in t)
                    if w not in lab: lab[w] = len(lab)
                    q.append(g)
            assert len(seen) == len(F) and len(lab) == nV
            c = tuple(sorted(tuple(sorted(lab[x] for x in f)) for f in F))
            if best is None or c < best: best = c
    return best

def perm_sign(t):
    return (-1) ** sum(1 for i in range(len(t)) for j in range(i + 1, len(t)) if t[i] > t[j])

def orient(F):
    """Coherent orientation by propagating ORDERED tetrahedra: the neighbour across the face opposite a is obtained by replacing a
    with the new vertex b in the same slot and reversing the sign.  Returns {sorted facet: eps}, eps = sign of the oriented
    tetrahedron relative to the increasing order."""
    F = [tuple(sorted(f)) for f in F]; tri = {}
    for f in F:
        for t in itertools.combinations(f, 3): tri.setdefault(t, []).append(f)
    assert all(len(v) == 2 for v in tri.values()), 'not a pseudomanifold'
    ordered = {F[0]: (F[0], 1)}; q = deque([F[0]])
    while q:
        f = q.popleft(); tup, s = ordered[f]
        for p in range(4):
            a = tup[p]; t = tuple(sorted(x for x in tup if x != a)); g = next(h for h in tri[t] if h != f)
            b = next(x for x in g if x not in t); ntup = tup[:p] + (b,) + tup[p + 1:]; ns = -s
            eps_new = ns * perm_sign(ntup)
            if g in ordered:
                gt, gs = ordered[g]; assert gs * perm_sign(gt) == eps_new, 'not orientable / inconsistent'
            else: ordered[g] = (ntup, ns); q.append(g)
    assert len(ordered) == len(F), 'not strongly connected'
    return {f: s * perm_sign(tup) for f, (tup, s) in ordered.items()}

def manifold_checks(F):
    """f-vector, Euler characteristic, vertex links are 2-spheres (connected closed surface with chi = 2)."""
    F = [tuple(sorted(f)) for f in F]; V = sorted(set(v for f in F for v in f))
    E = set(e for f in F for e in itertools.combinations(f, 2)); T = set(t for f in F for t in itertools.combinations(f, 3))
    ok = len(V) - len(E) + len(T) - len(F) == 0
    for v in V:
        L = [tuple(x for x in f if x != v) for f in F if v in f]
        le = {}
        for t in L:
            for e in itertools.combinations(t, 2): le.setdefault(e, []).append(t)
        lv = set(x for t in L for x in t)
        if not all(len(x) == 2 for x in le.values()) or len(lv) - len(le) + len(L) != 2: ok = False
        seen = {L[0]}; st = [L[0]]
        while st:
            t = st.pop()
            for e in itertools.combinations(t, 2):
                for g in le[e]:
                    if g not in seen: seen.add(g); st.append(g)
        if len(seen) != len(L): ok = False
    return ok, (len(V), len(E), len(T), len(F))
