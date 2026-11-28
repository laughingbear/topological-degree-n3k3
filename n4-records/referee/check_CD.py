r"""Claims C and D (combinatorial part). Exact integer arithmetic throughout.
usage: check_CD.py V12_polytope_deg11.json system_N12_deg11.json"""
import itertools, random, sys
from common import *


from hull import hom, hull_facets


def report_polytope(name, P, claimed=None):
    N = len(P); H = hom(P)
    fac, degen = hull_facets(P)
    F = sorted(G for G, _ in fac)
    zero5 = sum(1 for S in itertools.combinations(range(N), 5) if det([H[i] for i in S]) == 0)
    f = check_closed_manifold(F)
    print(f"[{name}] hull: {len(F)} simplicial facets, degenerate supporting 4-sets: {len(degen)}, "
          f"vanishing 5x5 minors of (1,p): {zero5} of {len(list(itertools.combinations(range(N),5)))}, f={f}, "
          f"vertices on hull: {len({v for G in F for v in G})}")
    if claimed is not None:
        print(f"[{name}] hull facets == claimed facet list: {F == sorted(tuple(sorted(t)) for t in claimed)}")
    # Does the geometric side sign agree with the combinatorial coherent orientation?
    o = orientation(F); side = dict(fac)
    rel = {o[i] * side[G] for i, G in enumerate(F)}
    print(f"[{name}] side sign chi(t,w) = eps * o_t with one global eps: {len(rel) == 1} (eps = {rel})")
    return F, o, side


dC = load(sys.argv[1]); dD = load(sys.argv[2])

# ---------------- Claim C ----------------
P = dC["polytope_points"]; Cmap = dC["C_int"]
F, o, side = report_polytope("C", P, dC["facets"])
# sanity of the whole degree machinery: the vertex map v -> p_v - c (c interior) must have degree ±1
c4 = [sum(p[k] for p in P) for k in range(4)]            # 12 * centroid
Pc = [[12 * p[k] - c4[k] for k in range(4)] for p in P]
s_id = facet_signs(F, o, Pc)
print("[C] control: radial map of the polytope itself: signs", set(s_id), "degree",
      {degree(F, s_id, Pc, u) for u in ([1, 2, 3, 5], [-3, 1, 7, 2], [5, -2, -1, 9])})
s = facet_signs(F, o, Cmap)
rng = random.Random(7)
vals = set()
for _ in range(80):
    u = [rng.randrange(-10**9, 10**9) for _ in range(4)]
    vals.add(degree(F, s, Cmap, u))
deg = vals.pop(); assert not vals
print(f"[C] vertex map degree = {deg}, q = {sum(1 for x in s if x == -sgn(deg))}")

# ---------------- Claim D ----------------
A = dD["A"]; C = dD["C"]; N = len(A)
assert len(C) == 4 and all(len(r) == N for r in C)
print(f"[D] A range [{min(min(a) for a in A)}, {max(max(a) for a in A)}], C range "
      f"[{min(min(r) for r in C)}, {max(max(r) for r in C)}], distinct exponents: {len({tuple(a) for a in A}) == N}")
Ccols = [[C[r][i] for r in range(4)] for i in range(N)]
minA = [abs(det([[1] + A[i] for i in S])) for S in itertools.combinations(range(N), 5)]
minC = [abs(det([Ccols[i] for i in S])) for S in itertools.combinations(range(N), 4)]
print(f"[D] (U_A): {all(minA)} (min |minor| {min(minA)}, {len(minA)} minors);  (U_C): {all(minC)} "
      f"(min |minor| {min(minC)}, {len(minC)} minors)")
FD, oD, sideD = report_polytope("D", A)
FC = sorted(tuple(sorted(t)) for t in dC["facets"])
print("[D] facets of conv(A) == facet list of Claim C (same labelling):", FD == FC)
print("[D] coefficient columns of D proportional to C_int of Claim C? ",
      all(all(Ccols[i][a] * Cmap[i][b] == Ccols[i][b] * Cmap[i][a] for a in range(4) for b in range(4)) for i in range(N)))
# chirotope of C (D) vs chirotope of C_int (claim C), restricted to facets
print("[D] sign det C_Gamma equals sign det C_int,Gamma on all facets (up to global sign):",
      len({sgn(det([Ccols[i] for i in G])) * sgn(det([Cmap[i] for i in G])) for G in FD}) == 1)

def bihan(u):
    """deg F = (-1)^n sum_{Gamma, u in cone C_Gamma} sign det C_Gamma * sign det Ahat_{Gamma ∪ {i}}.
    Ahat has COLUMNS (1,a_i); det of the transpose is the same, so sideD is exactly eps_Gamma."""
    tot = 0; used = 0
    for G in FD:
        vecs = [Ccols[i] for i in G]
        if in_open_cone(vecs, u):
            tot += sgn(det(vecs)) * sideD[G]; used += 1
    return (-1) ** 4 * tot, used

def generic(u):
    return all(det([Ccols[i] for i in S] + [u]) != 0 for S in itertools.combinations(range(N), 3))

us = [[1, 2, 3, 5], [-3, 1, 7, 2], [5, -2, -1, 9], [1, 0, 0, 0], [0, 0, 0, -1]]
us += [[rng.randrange(-10**6, 10**6) for _ in range(4)] for _ in range(300)]
res = {}; cov = set(); ng = 0
for u in us:
    if not generic(u):
        ng += 1; continue
    d, used = bihan(u); res[d] = res.get(d, 0) + 1; cov.add(used)
print(f"[D] Bihan formula over {sum(res.values())} generic u ({ng} non-generic skipped): values {res}; "
      f"number of facets whose cone contains u ranges over {sorted(cov)}")
# every proper face in an open half-space (H)?  Under U_A,U_C it suffices for facets: 4 independent
# vectors always lie in an open half-space.  Recorded for completeness.
sD = facet_signs(FD, oD, Ccols)
dg = {degree(FD, sD, Ccols, u) for u in us if generic(u)}
print(f"[D] same thing via CLAIMS.md degree definition with recomputed orientation: {dg}; q = "
      f"{sum(1 for x in sD if x == -sgn(next(iter(dg))))}")
