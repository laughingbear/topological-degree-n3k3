"""Build fewnomial systems (n = 3) attaining |deg F| = floor((3N-8)/4) for any N >= 5 by gluing gadgets.

Gluing: host polytope P (degree d>0, a facet t with positive image triangle), gadget cap Q (7 vertices, degree +3,
exactly one negative image triangle t_g). Q is realized as a flat cap over its facet t_g; an affine map puts t_g onto t
(outward-flattened), a linear map M (det>0) sends the gadget's coefficient columns at t_g onto the host's at t.
Facets t and t_g disappear; vertices: V + 4; degree: d + 3 (degree is additive since the removed image triangles
coincide with opposite orientation). Exact verification of every output with Bihan's Theorem A (verify_exact.check)."""
import numpy as np, json, sys, os, itertools, networkx as nx
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
sys.path.insert(0, HERE)
from fractions import Fraction
from scipy.spatial import ConvexHull
from verify_exact import check
from deg import solid

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)

def outward_facets(A):
    h = ConvexHull(A); c = A.mean(0); F = []
    for s in h.simplices:
        i, j, k = (int(x) for x in s)
        if np.dot(np.cross(A[j] - A[i], A[k] - A[i]), A[i] - c) < 0: j, k = k, j
        F.append((i, j, k))
    return F, h

def degree(F, C):
    return sum(solid(C[a], C[b], C[c]) for a, b, c in F) / (4 * np.pi)

def sig(C, t):
    a, b, c = t; return np.sign(np.linalg.det(np.array([C[a], C[b], C[c]])))

# ---------- gadget: realize its triangulation as a flat cap over the negative triangle ----------
G = json.load(open(os.path.join(DATA, 'gadget.json')))
gC = np.array(G['C'], float)                   # integer gadget, degree +3, one negative triangle
gT = [tuple(int(x) for x in t) for t in G['T']]
assert round(degree(gT, gC)) == 3
neg = [t for t in gT if sig(gC, t) < 0]; assert len(neg) == 1
tg = neg[0]                                    # outward-oriented (in gadget sphere) negative triangle
others = [v for v in range(7) if v not in tg]
Gg = nx.Graph(); [Gg.add_edges_from([(a, b), (b, c), (c, a)]) for a, b, c in gT]

def realize_cap():
    base = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    for _ in range(200000):
        # 4 points above the base triangle, projecting inside it
        u = rng.dirichlet([1, 1, 1], size=4); z = rng.uniform(0.05, 1, size=4)
        P = np.vstack([base, np.column_stack([u[:, 1], u[:, 2], z])])
        F, h = outward_facets(P)
        if len(h.vertices) < 7: continue
        H = nx.Graph(); [H.add_edges_from([(a, b), (b, c), (c, a)]) for a, b, c in F]
        if frozenset((0, 1, 2)) not in {frozenset(f) for f in F}: continue
        gm = nx.isomorphism.GraphMatcher(Gg, H)
        for phi in gm.isomorphisms_iter():          # gadget label -> cap point index
            if {phi[tg[0]], phi[tg[1]], phi[tg[2]]} != {0, 1, 2}: continue
            # faces must match (graph iso of a triangulation preserves faces) and orientation must match
            mapped = {frozenset(phi[x] for x in t) for t in gT}
            if mapped != {frozenset(f) for f in F}: continue
            # orientation check: outward orientation of cap facet == image of gadget's oriented triangle
            Fo = {frozenset(f): f for f in F}
            ok = True
            for t in gT:
                f = Fo[frozenset(phi[x] for x in t)]; img = tuple(phi[x] for x in t)
                rot = [img, img[1:] + img[:1], img[2:] + img[:2]]
                if f not in rot: ok = False; break
            if ok:
                return P, phi
    raise RuntimeError('no cap realization')

capP, phi = realize_cap()
print('gadget cap realized', flush=True)

# ---------- host bases ----------
def load(fn):
    S = json.load(open(fn)); return np.array(S['A'], float), np.array(S['C'], float).T
def base(V):
    """base systems from the repository data (all exactly verified)."""
    if V == 7: return load(os.path.join(DATA, 'system_deg3.json'))
    if V == 8: return load(os.path.join(DATA, 'system_N8_deg4.json'))
    S = json.load(open(os.path.join(DATA, 'sharp_family.json')))[str(V)]
    return np.array(S['A'], float), np.array(S['C'], float).T

def glue(A, C):
    """glue one gadget onto a positively oriented facet of host (A, C) with positive degree."""
    F, h = outward_facets(A)
    t = next(f for f in F if sig(C, f) > 0)
    a, b, c = t
    # gadget tg=(x,y,z) outward in gadget; glue x->a, y->c, z->b (orientation reversal)
    x, y, z = tg
    target = {x: a, y: c, z: b}
    # linear map M: M gC[x] = C[a], M gC[y] = C[c], M gC[z] = C[b]
    Sg = np.column_stack([gC[x], gC[y], gC[z]]); Sh = np.column_stack([C[a], C[c], C[b]])
    M = Sh @ np.linalg.inv(Sg); assert np.linalg.det(M) > 0
    # affine placement of cap: base point phi[x],phi[y],phi[z] -> A[a],A[c],A[b]; height -> outward normal * eps
    nrm = np.cross(A[b] - A[a], A[c] - A[a]); nrm /= np.linalg.norm(nrm)
    if np.dot(nrm, A[a] - A.mean(0)) < 0: nrm = -nrm
    best = None
    for eps in [0.3, 0.1, 0.03, 0.01, 0.003]:
        newA = list(A); newC = list(C)
        for v in others:
            p = capP[phi[v]]                        # (u1,u2,h): barycentric-ish coords over base (0,0),(1,0),(0,1)
            bx = np.array([1 - p[0] - p[1], p[0], p[1]])   # weights of base points 0,1,2
            basepts = {phi[x]: A[a], phi[y]: A[c], phi[z]: A[b]}
            q = bx[0] * basepts[0] + bx[1] * basepts[1] + bx[2] * basepts[2] + eps * p[2] * np.linalg.norm(A[b] - A[a]) * nrm
            newA.append(q); newC.append(M @ gC[v])
        newA = np.array(newA); newC = np.array(newC)
        Fn, hn = outward_facets(newA)
        if len(hn.vertices) == len(newA) and round(degree(Fn, newC)) == round(degree(F, C)) + 3:
            return newA, newC
    raise RuntimeError('gluing failed')

def integerize(A, C, sa, sc):
    Ai = np.round(A * sa).astype(int); Ci = np.round(C / np.abs(C).max() * sc).astype(int)
    Ai = Ai - Ai.min(0)
    return Ai, Ci

if __name__ == '__main__':
    Nmax = int(sys.argv[1]); out = {}
    for N in range(5, Nmax + 1):
        V0 = {0: 8, 1: 5, 2: 6, 3: 7}[N % 4]
        if N < V0: continue
        A, C = base(V0)
        if round(degree(outward_facets(A)[0], C)) < 0: C = -C
        for _ in range((N - V0) // 4): A, C = glue(A, C)
        target = (3 * N - 8) // 4
        Ff, _ = outward_facets(A); dfl = degree(Ff, C)
        done = False
        for sa in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096]:
            for sc in [100, 1000, 10000, 100000]:
                Ai, Ci = integerize(A * (1.0 / np.abs(A).max()) * 20, C, sa, sc)
                if len(ConvexHull(Ai.astype(float)).vertices) < N: continue
                res = []
                for u in ([7, 11, -13], [1, 2, 3], [-5, 3, 8]):
                    try: res.append(check(Ai.tolist(), Ci.T.tolist(), u))
                    except ValueError: pass
                ds = {r[3] for r in res}
                UA = all(r[0] for r in res); UC = all(r[1] for r in res); d = next(iter(ds)) if len(ds) == 1 else None
                if res and UA and UC and d is not None and abs(d) == target:
                    out[N] = {'A': Ai.tolist(), 'C': Ci.T.tolist(), 'deg': d}
                    print(f'N={N}: target {target}, float deg {dfl:.4f}, exact deg {d}, scale A {sa}, C {sc}', flush=True)
                    done = True; break
            if done: break
        if not done: print(f'N={N}: float deg {dfl:.4f}, integer certification failed', flush=True)
    json.dump(out, open(sys.argv[3] if len(sys.argv) > 3 else 'sharp_family_rebuilt.json', 'w'))
