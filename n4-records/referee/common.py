"""Referee code, written from CLAIMS.md section 0. No code or stored signs from ours/ are used.
All arithmetic is on Python ints (exact)."""
import itertools, json
from collections import defaultdict


def det(M):
    """Exact determinant of a square integer matrix (list of rows) by Laplace expansion."""
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    tot = 0
    for j in range(n):
        if M[0][j] == 0:
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        tot += (-1) ** j * M[0][j] * det(minor)
    return tot


def sgn(x):
    return (x > 0) - (x < 0)


def load(path):
    return json.load(open(path))


def faces(facets, k):
    out = defaultdict(list)
    for idx, f in enumerate(facets):
        for s in itertools.combinations(f, k):
            out[s].append(idx)
    return out


def check_closed_manifold(facets):
    """Returns f-vector. Raises AssertionError if not a closed connected combinatorial
    3-manifold candidate: every triangle in exactly 2 facets, every vertex link a 2-sphere
    (connected closed surface with Euler characteristic 2), every edge link a single cycle."""
    facets = [tuple(sorted(f)) for f in facets]
    assert all(len(set(f)) == 4 for f in facets), "facet with repeated vertex"
    assert len(set(facets)) == len(facets), "duplicate facet"
    V = sorted({v for f in facets for v in f})
    assert V == list(range(len(V))), "vertex labels not 0..V-1"
    tri = faces(facets, 3)
    assert all(len(v) == 2 for v in tri.values()), "triangle not in exactly 2 facets"
    edges = faces(facets, 2)
    # vertex links
    for v in V:
        L = [tuple(x for x in f if x != v) for f in facets if v in f]
        lv = sorted({x for t in L for x in t})
        le = faces(L, 2)
        assert all(len(x) == 2 for x in le.values()), f"link of {v}: edge not in 2 triangles"
        chi = len(lv) - len(le) + len(L)
        assert chi == 2, f"link of {v} has Euler characteristic {chi}"
        assert connected(lv, list(le)), f"link of {v} disconnected"
        # each vertex of the link has a single-cycle link (link is a surface, not pinched)
        for w in lv:
            cyc = [tuple(x for x in t if x != w) for t in L if w in t]
            assert single_cycle(cyc), f"link of edge {v},{w} is not one cycle"
    # global connectivity (strongly connected through triangles)
    adj = defaultdict(list)
    for a, b in tri.values():
        adj[a].append(b); adj[b].append(a)
    seen = {0}; st = [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); st.append(y)
    assert len(seen) == len(facets), "not strongly connected"
    f = (len(V), len(edges), len(tri), len(facets))
    assert f[0] - f[1] + f[2] - f[3] == 0
    return f


def connected(verts, edges):
    adj = defaultdict(list)
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    seen = {verts[0]}; st = [verts[0]]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y); st.append(y)
    return len(seen) == len(verts)


def single_cycle(edges):
    vs = sorted({x for e in edges for x in e})
    deg = defaultdict(int)
    for a, b in edges:
        deg[a] += 1; deg[b] += 1
    return len(edges) >= 3 and all(deg[v] == 2 for v in vs) and connected(vs, edges)


def orientation(facets):
    """Coherent orientation per CLAIMS.md: o_t (-1)^i = - o_t' (-1)^j across each triangle.
    Normalised to o = +1 on facets[0]. Verifies consistency on every triangle."""
    facets = [tuple(sorted(f)) for f in facets]
    tri = defaultdict(list)
    for idx, f in enumerate(facets):
        for i in range(4):
            tri[f[:i] + f[i + 1:]].append((idx, i))
    nb = defaultdict(list)
    for (t, i), (t2, j) in tri.values():
        rel = -((-1) ** i) * ((-1) ** j)  # o_t2 = rel * o_t
        nb[t].append((t2, rel)); nb[t2].append((t, rel))
    o = {0: 1}; st = [0]
    while st:
        x = st.pop()
        for y, rel in nb[x]:
            if y not in o:
                o[y] = rel * o[x]; st.append(y)
    assert len(o) == len(facets)
    for (t, i), (t2, j) in tri.values():
        assert o[t] * (-1) ** i == -o[t2] * (-1) ** j, "not orientable"
    return [o[i] for i in range(len(facets))]


def facet_signs(facets, o, C):
    s = []
    for f, ot in zip(facets, o):
        d = det([C[v] for v in sorted(f)])
        assert d != 0, f"facet {f}: dependent vectors"
        s.append(ot * sgn(d))
    return s


def in_open_cone(vecs, u):
    """u in open cone(vecs) (4 independent integer vectors)? Cramer, exact.
    Returns True/False; raises ValueError if u is on the boundary of the cone's span walls
    (some coordinate zero), i.e. u not generic for this facet."""
    D = det(vecs)
    inside = True
    for k in range(4):
        M = list(vecs); M[k] = u
        Dk = det(M)
        if Dk == 0:
            raise ValueError("u not generic")
        if sgn(Dk) != sgn(D):
            inside = False
    return inside


def degree(facets, s, C, u):
    tot = 0
    for f, st in zip(facets, s):
        if in_open_cone([C[v] for v in sorted(f)], u):
            tot += st
    return tot
