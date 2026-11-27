"""Check the mixed-sign identity and bound for vertex maps on simplicial 3-spheres (n = 4).
  4 deg = sum_v w_v - X + Y,   sum_v w_v = sum_e w_e - (p - q),
  X (Y) = sum over positive (negative) facets of D_t,  D_t = measure of the 4 orthants of t with exactly one negative coordinate,
  bound: 4 deg < 2F - E/2 - q/2   (deg > 0, q = # negative facets).
"""
import json, sys, itertools, numpy as np
from scipy.spatial import ConvexHull
rng = np.random.default_rng(1)

def orient(facets):
    facets = [tuple(sorted(f)) for f in facets]
    ridge = {}
    for i, f in enumerate(facets):
        for k in range(4):
            ridge.setdefault(f[:k] + f[k+1:], []).append((i, (-1) ** k))
    o = {0: 1}; stack = [0]
    while stack:
        i = stack.pop(); f = facets[i]
        for k in range(4):
            r = f[:k] + f[k+1:]; (a, sa), (b, sb) = ridge[r]
            j, sj, si = (b, sb, sa) if a == i else (a, sa, sb)
            want = -o[i] * si * sj
            if j in o: assert o[j] == want
            else: o[j] = want; stack.append(j)
    return facets, [o[i] for i in range(len(facets))]

def analyse(facets, C, name, nmc=400000):
    C = np.array(C, float); V = C.shape[1]
    facets, o = orient(facets)
    F = len(facets); M = [C[:, list(f)] for f in facets]
    s = [o[i] * int(np.sign(np.linalg.det(M[i]))) for i in range(F)]
    Minv = [np.linalg.inv(m) for m in M]
    def count(u, t_idx, free):   # signed count of facets whose coefficients at positions not in 'free' are > 0
        tot = 0
        for i in t_idx:
            c = Minv[i] @ u
            if all(c[k] > 0 for k in range(4) if facets[i][k] not in free): tot += s[i]
        return tot
    degs = {count(rng.standard_normal(4), range(F), ()) for _ in range(40)}
    assert len(degs) == 1, degs; deg = degs.pop()
    if deg < 0: s = [-x for x in s]; deg = -deg
    p = s.count(1); q = s.count(-1)
    wv = []
    for v in range(V):
        idx = [i for i in range(F) if v in facets[i]]
        vals = {count(rng.standard_normal(4), idx, (v,)) for _ in range(12)}; assert len(vals) == 1; wv.append(vals.pop())
    edges = sorted({e for f in facets for e in itertools.combinations(f, 2)}); E = len(edges)
    we = []
    for e in edges:
        idx = [i for i in range(F) if e[0] in facets[i] and e[1] in facets[i]]
        vals = {count(rng.standard_normal(4), idx, e) for _ in range(12)}; assert len(vals) == 1; we.append(vals.pop())
    U = rng.standard_normal((nmc, 4))
    X = Y = 0.0; Dmax = 0; vol = 0.0
    for i in range(F):
        neg = ((U @ Minv[i].T) < 0).sum(1)
        D = (neg == 1).mean(); vol += s[i] * (neg == 0).mean(); Dmax = max(Dmax, D)
        if s[i] > 0: X += D
        else: Y += D
    assert F == E - V
    bound = 2 * F - E / 2 - q / 2
    print(f"{name}: V={V} E={E} F={F} p={p} q={q} deg={deg} (MC vol {vol:.3f})  sum w_v={sum(wv)}  sum w_e={sum(we)}  "
          f"Gram: {sum(wv)} = {sum(we) - (p - q)}  X={X:.3f} Y={Y:.3f}  sum_w - X + Y = {sum(wv) - X + Y:.3f} vs 4deg={4*deg}  "
          f"Dmax={Dmax:.3f}  bound 4deg < {bound}  {'OK' if 4*deg < bound else 'VIOLATION'}")
    assert sum(wv) == sum(we) - (p - q)
    assert abs(sum(wv) - X + Y - 4 * deg) < 0.15, "identity fails"
    assert 4 * deg < bound and Dmax < 0.5

def hull_facets(A):
    return [tuple(sorted(int(x) for x in f)) for f in ConvexHull(np.array(A, float)).simplices]

if __name__ == "__main__":
    T = "data/examples/"
    d = json.load(open(T + "v10/vertexmap_deg9_nonpolytopal.json"))
    analyse([f for f, _ in d["oriented_facets"]], d["C"], "deg9 sphere101272")
    d = json.load(open(T + "v10/system_deg8_N10_small.json")); analyse(hull_facets(d["A"]), d["C"], "deg8 N=10")
    d = json.load(open(T + "v10/D45_family.json"))
    for k, C in d["C_by_degree"].items(): analyse(hull_facets(d["A"]), C, f"D45 family deg {k}", 100000)
    d = json.load(open(T + "deg6_systems.json"))
    for k, sy in d.items(): analyse(hull_facets(sy["A"]), sy["C"], f"deg6 sys {k}", 100000)
    # random maps on random polytopes
    for trial in range(30):
        V = int(rng.integers(6, 14)); A = rng.standard_normal((V, 4)); A /= np.linalg.norm(A, axis=1)[:, None]; C = rng.integers(-9, 10, (4, V))
        try: analyse(hull_facets(A), C, f"random {trial}", 60000)
        except np.linalg.LinAlgError: pass
    print("ALL OK")
