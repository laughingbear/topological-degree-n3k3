"""Item 5(b): test Theorem M (identity of step 2, relation of step 4, edge bound of step 5, D_t bounds,
F = E - V, and 4 deg < 2F - E/2 - q/2) on random vertex maps on boundaries of random simplicial
4-polytopes with 6..14 vertices, and on the example files in topdeg4/ (JSON read directly).

Exact (integer) arithmetic: facet signs, deg, w_v.  Floating point: w_e (angle sum, checked near-integer),
mu(L_{t,v}) (closed-form 3-d solid angle), mu(cone t) and D_t (Genz orthant probabilities).
"""
import sys, json, itertools, math, time
import numpy as np
from scipy.spatial import ConvexHull
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from core import (VertexMap, det, colmat, sgn, solid_angle_3, orth_basis_perp, orthant_probs_4)
from scipy.stats import multivariate_normal

TOL_GENZ = 1e-6


def orthant(S, pat):
    Dg = np.diag(pat)
    S2 = Dg @ S @ Dg
    dd = 1 / np.sqrt(np.diag(S2))
    R = S2 * dd[:, None] * dd[None, :]          # correlation matrix (orthant prob is scale invariant)
    R = (R + R.T) / 2
    return float(multivariate_normal.cdf(np.zeros(4), mean=np.zeros(4), cov=R, allow_singular=True, abseps=TOL_GENZ,
                                         releps=TOL_GENZ, maxpts=2_000_000))


def facet_measures(M):
    """returns mu(cone), D_t, and (for checking) the 1-neg orthants"""
    Mi = np.linalg.inv(M)
    S = Mi @ Mi.T
    P = orthant(S, (1, 1, 1, 1))
    one = []
    for k in range(4):
        pat = [1] * 4
        pat[k] = -1
        one.append(orthant(S, pat))
    return P, sum(one), one


def analyse(vm, rng, name="", verbose=False):
    V, E, F = vm.V, vm.E, vm.F
    res = {"name": name, "V": V, "E": E, "F": F}
    deg = vm.degree(rng)
    if deg < 0:
        vm.flip()
        deg = -deg
    p = vm.s.count(1)
    q = vm.s.count(-1)
    res.update(deg=deg, p=p, q=q)
    # F = E - V
    res["F=E-V"] = (F == E - V)
    # link degrees, exact
    wv = [vm.link_degree(v, rng) for v in range(V)]
    # edges
    we = {}
    for e in vm.E_list:
        we[e] = vm.edge_winding(e)
    edge_bound_ok = all(w <= (pe - 1) / 2 for (w, pe, me) in we.values())
    res["edge_bound_ok"] = edge_bound_ok
    # per-vertex Gauss-Bonnet 4 w_v = 2 sum_{e ni v} w_e - (p_v - q_v)
    gb_ok = True
    for v in range(V):
        pv = sum(1 for i, f in enumerate(vm.facets) if v in f and vm.s[i] > 0)
        qv = sum(1 for i, f in enumerate(vm.facets) if v in f and vm.s[i] < 0)
        se = sum(w for e, (w, _, _) in we.items() if v in e)
        gb_ok &= (4 * wv[v] == 2 * se - (pv - qv))
    res["GB_vertex_ok"] = gb_ok
    sum_wv = sum(wv)
    sum_we = sum(w for (w, _, _) in we.values())
    res["GB_global_ok"] = (sum_wv == sum_we - (p - q))
    # measures
    Cf = [np.array(c, float) for c in vm.C]
    X = Y = 0.0
    vol = 0.0
    Dmin, Dmax = 1, 0
    per_facet_err = 0.0
    Lsum = [0.0] * V
    for i, f in enumerate(vm.facets):
        M = np.column_stack([Cf[j] for j in f])
        P, D, _ = facet_measures(M)
        vol += vm.s[i] * P
        Dmin, Dmax = min(Dmin, D), max(Dmax, D)
        if vm.s[i] > 0:
            X += D
        else:
            Y += D
        Ls = 0.0
        for v in f:
            B = orth_basis_perp(Cf[v])
            pr = [B.T @ Cf[j] for j in f if j != v]
            L = solid_angle_3(*pr)
            Ls += L
            Lsum[v] += vm.s[i] * L
        per_facet_err = max(per_facet_err, abs(Ls - (4 * P + D)))
    res["link_measure_err"] = max(abs(Lsum[v] - wv[v]) for v in range(V))
    res["vol_err"] = abs(vol - deg)
    res["per_facet_err"] = per_facet_err
    res["identity_err"] = abs(4 * deg - (sum_wv - X + Y))
    res.update(X=X, Y=Y, Dmin=Dmin, Dmax=Dmax, sum_wv=sum_wv, sum_we=sum_we)
    res["bound"] = 2 * F - E / 2 - q / 2
    res["bound_ok"] = 4 * deg < res["bound"]
    res["slack"] = res["bound"] - 4 * deg
    res["ok"] = (res["F=E-V"] and edge_bound_ok and gb_ok and res["GB_global_ok"] and res["bound_ok"]
                 and res["link_measure_err"] < 1e-9 and res["vol_err"] < 1e-4 and res["per_facet_err"] < 5e-5
                 and res["identity_err"] < 1e-3 and 0 < Dmin and Dmax < 0.5)
    return res


# ---------------------------------------------------------------- random polytopes
def random_polytope(V, rng, scale=1000):
    while True:
        X = rng.standard_normal((V, 4))
        X = X / np.linalg.norm(X, axis=1)[:, None]
        A = [[int(round(scale * x)) for x in row] for row in X]
        try:
            h = ConvexHull(np.array(A, float))
        except Exception:
            continue
        if len(h.vertices) != V:
            continue
        facets = sorted({tuple(sorted(map(int, s))) for s in h.simplices})
        # exact simpliciality: no 5 points of a facet hyperplane
        ok = True
        for f in facets:
            for z in range(V):
                if z in f:
                    continue
                if det([[1] + A[j] for j in list(f) + [z]]) == 0:
                    ok = False
        if ok:
            return A, facets


def hull_facets_exact(A):
    h = ConvexHull(np.array(A, float))
    assert len(h.vertices) == len(A), "not all points are vertices"
    facets = sorted({tuple(sorted(map(int, s))) for s in h.simplices})
    V = len(A)
    for f in facets:
        # all other points strictly on one side (exact)
        sides = {sgn(det([[1] + list(A[j]) for j in list(f) + [z]])) for z in range(V) if z not in f}
        assert 0 not in sides and len(sides) == 1, "hull facet not exact/simplicial"
    return facets


def random_C(V, rng, R=12):
    return [[int(x) for x in rng.integers(-R, R + 1, 4)] for _ in range(V)]


def family_C(kind, A, rng):
    """three families: 'uniform' random integers; 'near-identity' C_v = 3 A_v/1000*... + noise (deg +-1,
    mixed signs after noise); 'trig' C_v = R*(cos k t, sin k t, cos l t, sin l t) at random t_v, random
    k<l, rounded (tends to give large |deg|)."""
    V = len(A)
    if kind == "uniform":
        return random_C(V, rng)
    if kind == "near":
        lam = rng.uniform(0.02, 0.08)
        return [[int(round(lam * a + rng.integers(-40, 41))) for a in row] for row in A]
    k = int(rng.integers(1, 4)); l = int(rng.integers(k + 1, 8))
    out = []
    for _ in range(V):
        t = rng.uniform(0, 2 * np.pi)
        out.append([int(round(1000 * x)) for x in (np.cos(k * t), np.sin(k * t), np.cos(l * t), np.sin(l * t))])
    return out


def C_columns(Cjson):
    """JSON stores C as 4 rows x V columns"""
    C = np.array(Cjson)
    assert C.shape[0] == 4
    return [[int(x) for x in C[:, j]] for j in range(C.shape[1])]


def examples(T):
    out = []
    d = json.load(open(T + "/v10/vertexmap_deg9_nonpolytopal.json"))
    facets = [tuple(f) for f, _ in d["oriented_facets"]]
    out.append(("deg9 nonpolytopal sphere", facets, C_columns(d["C"]), d["oriented_facets"]))
    d = json.load(open(T + "/v10/system_deg8_N10_small.json"))
    out.append(("system deg8 N10", hull_facets_exact(d["A"]), C_columns(d["C"]), None))
    d = json.load(open(T + "/v10/D45_family.json"))
    Fh = hull_facets_exact(d["A"])
    for k, C in d["C_by_degree"].items():
        out.append((f"D45 family claimed deg {k}", Fh, C_columns(C), None))
    d = json.load(open(T + "/deg6_systems.json"))
    for k, sy in d.items():
        out.append((f"deg6 system {k}", hull_facets_exact(sy["A"]), C_columns(sy["C"]), None))
    return out


if __name__ == "__main__":
    T = sys.argv[1] if len(sys.argv) > 1 else None
    ntrials = int(sys.argv[2]) if len(sys.argv) > 2 else 220
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 12345
    rng = np.random.default_rng(seed)
    bad = []
    t0 = time.time()
    if T and T != "-":
        for name, facets, C, ofac in examples(T):
            vm = VertexMap(facets, C)
            if ofac is not None:  # compare supplied orientation with mine (up to global sign)
                mine = dict(zip(vm.facets, vm.o))
                rel = {mine[tuple(sorted(f))] * o for f, o in ofac}
                note = f" supplied orientation coherent: {len(rel) == 1}"
            else:
                note = ""
            r = analyse(vm, rng, name)
            print(f"{name}: V={r['V']} E={r['E']} F={r['F']} p={r['p']} q={r['q']} deg={r['deg']} "
                  f"sum_wv={r['sum_wv']} sum_we={r['sum_we']} X={r['X']:.4f} Y={r['Y']:.4f} "
                  f"idErr={r['identity_err']:.1e} Dmax={r['Dmax']:.3f} 4deg={4*r['deg']} < {r['bound']} "
                  f"{'OK' if r['ok'] else 'FAIL ' + str(r)}{note}", flush=True)
            if not r["ok"]:
                bad.append(r)
    stats = []
    nskip_cond = [0]
    trial = 0
    while trial < ntrials:
        V = int(rng.integers(6, 15))
        A, facets = random_polytope(V, rng)
        kind = ["uniform", "near", "trig"][trial % 3]
        C = family_C(kind, A, rng)
        try:
            vm = VertexMap(facets, C)
        except ValueError:
            continue
        cond = max(np.linalg.cond(np.array([vm.C[j] for j in f], float).T) for f in vm.facets)
        if cond > 1e6:   # badly conditioned facet basis: exact parts fine, Genz unreliable -> resample
            nskip_cond[0] += 1
            continue
        r = analyse(vm, rng, f"random {trial} {kind}")
        r["kind"] = kind
        stats.append(r)
        if not r["ok"]:
            bad.append(r)
            print("FAIL", r, flush=True)
        trial += 1
        if trial % 20 == 0:
            print(f"{trial} random maps done ({time.time()-t0:.0f}s); max identity err "
                  f"{max(s['identity_err'] for s in stats):.1e}", flush=True)
    import collections
    print("random maps:", len(stats), " V distribution:", dict(collections.Counter(s["V"] for s in stats)))
    print("  with q>0 and p>0 (mixed signs):", sum(1 for s in stats if s["q"] > 0 and s["p"] > 0))
    for kd in ("uniform", "near", "trig"):
        print("  family", kd, ":", sum(1 for s in stats if s.get("kind") == kd), "maps, deg distribution",
              dict(sorted(collections.Counter(s["deg"] for s in stats if s.get("kind") == kd).items())))
    print("  deg distribution:", dict(collections.Counter(s["deg"] for s in stats)))
    print("  max identity err %.2e, max vol err %.2e, max per-facet err %.2e, max link-measure err %.2e" % (
        max(s["identity_err"] for s in stats), max(s["vol_err"] for s in stats),
        max(s["per_facet_err"] for s in stats), max(s["link_measure_err"] for s in stats)))
    print("  min slack of bound:", min(s["slack"] for s in stats), " max Dmax:", max(s["Dmax"] for s in stats))
    print("resampled because a facet basis had condition number > 1e6:", nskip_cond[0])
    print("FAILURES:", len(bad))
