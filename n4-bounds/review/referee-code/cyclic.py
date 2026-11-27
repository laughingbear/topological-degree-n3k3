"""Item 5(a): Theorem C by direct cone counting, independent of the attached scripts.
 - facets of C(V,4) generated from Gale's evenness condition (not from the {i,i+1,j,j+1} recipe),
   and cross-checked against a Qhull convex hull of V points on the trigonometric moment curve
   (cos t, sin t, cos 2t, sin 2t), which is combinatorially cyclic;
 - orientation from my own coherent-orientation routine (no moment-curve sign formula);
 - C_i = (-1)^i gamma(pi i / V), gamma(t) = (cos t, sin t, cos 3t, sin 3t), in 60-digit arithmetic;
 - degree = signed count of facet cones containing u, for several random u AND for u = gamma(tau),
   tau in (0, pi/V); every determinant is checked to be far from 0 relative to working precision.
"""
import itertools, sys
import mpmath as mp
import numpy as np
from scipy.spatial import ConvexHull
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from core import coherent_orientation

mp.mp.dps = 60


def gale_facets(V):
    out = []
    for S in itertools.combinations(range(V), 4):
        Sset = set(S)
        ok = True
        non = [x for x in range(V) if x not in Sset]
        for a, b in zip(non, non[1:]):
            if sum(1 for x in S if a < x < b) % 2:
                ok = False
                break
        if ok:
            out.append(S)
    return out


def hull_facets_trig(V):
    t = 2 * np.pi * np.arange(V) / V
    P = np.column_stack([np.cos(t), np.sin(t), np.cos(2 * t), np.sin(2 * t)])
    h = ConvexHull(P)
    return sorted(tuple(sorted(map(int, s))) for s in h.simplices)


def recipe_facets(V):
    S = set()
    for i in range(V):
        for j in range(V):
            f = {i, (i + 1) % V, j, (j + 1) % V}
            if len(f) == 4:
                S.add(tuple(sorted(f)))
    return sorted(S)


def gamma(t):
    return [mp.cos(t), mp.sin(t), mp.cos(3 * t), mp.sin(3 * t)]


def run(V, rng, check_hull=True):
    m = (V - 1) // 2
    G = sorted(gale_facets(V))
    assert len(G) == V * (V - 3) // 2
    assert G == recipe_facets(V), "recipe facet list != Gale facet list"
    if check_hull:
        H = hull_facets_trig(V)
        # trig-moment hull is cyclic, but with vertex labels in cyclic order; Gale list is invariant
        # under cyclic shift for even d, so compare directly
        assert H == G, f"hull facets differ from Gale facets at V={V}"
    facets, o = coherent_orientation(G)
    C = [[(-1) ** i * x for x in gamma(mp.pi * i / V)] for i in range(V)]
    Ms = [mp.matrix([[C[j][r] for j in f] for r in range(4)]) for f in facets]
    ds = [mp.det(M) for M in Ms]
    mind = min(abs(d) for d in ds)
    assert mind > mp.mpf(10) ** -40
    s = [o[i] * (1 if ds[i] > 0 else -1) for i in range(len(facets))]
    same = abs(sum(s)) == len(s)
    Minv = [M ** -1 for M in Ms]

    def count(u):
        tot = 0
        minc = mp.inf
        for i in range(len(facets)):
            c = Minv[i] * u
            minc = min(minc, min(abs(x) for x in c))
            if all(x > 0 for x in c):
                tot += s[i]
        assert minc > mp.mpf(10) ** -30, "u too close to a triangle cone"
        return tot

    degs = set()
    for _ in range(5):
        u = mp.matrix([mp.mpf(float(x)) for x in rng.standard_normal(4)])
        degs.add(count(u))
    for tau in [mp.pi / V * mp.mpf(k) / 7 for k in (1, 3, 6)]:
        degs.add(count(mp.matrix(gamma(tau))))
    assert len(degs) == 1, degs
    deg = degs.pop()
    return dict(V=V, F=len(facets), same_sign=same, deg=deg, pred=(V - 1) * (V - 3) // 8,
                wrapsign=[s[i] for i, f in enumerate(facets) if f[0] == 0 and f[3] == V - 1][:1])


if __name__ == "__main__":
    rng = np.random.default_rng(2026)
    allok = True
    for V in range(5, 42, 2):
        r = run(V, rng)
        ok = r["same_sign"] and abs(r["deg"]) == r["pred"]
        allok &= ok
        print(f"V={V:2d} F={r['F']:4d} same_sign={r['same_sign']} deg={r['deg']:4d} "
              f"(V-1)(V-3)/8={r['pred']:4d}  {'OK' if ok else 'MISMATCH'}", flush=True)
    print("ALL OK" if allok else "FAILURES")
