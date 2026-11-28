r"""Controls for the Claim-B encoder.
1. Real polytopes (Claim C points; conv(A) of Claim D; cyclic polytope C(12,4) on the moment
   curve; a random simplicial polytope): the ACTUAL chirotope computed from coordinates must
   satisfy every unit condition (with a single eps) and every GP relation of the encoder, and the
   propagator must report no conflict and never derive a sign different from the true one.
2. Negative control: the same complexes with two vertex labels of the ORIENTATION input swapped
   are not tested here; instead see REVIEW.md for the Barnette-sphere remark.
"""
import itertools, random, sys
from common import *
from chiro import gp_terms, unit_conditions, propagate
from hull import hull_facets


def real_chirotope(P):
    H = [[1] + list(p) for p in P]
    return {S: sgn(det([H[i] for i in S])) for S in itertools.combinations(range(len(P)), 5)}


def control(name, P):
    N = len(P)
    fac, degen = hull_facets(P)
    F = sorted(G for G, _ in fac); assert not degen
    chi = real_chirotope(P); assert all(chi.values()), "points not in general position"
    o = orientation(F)
    ok_eps = []
    for eps in (1, -1):
        if all(chi[k] == v for k, v, _ in unit_conditions(F, o, eps)):
            ok_eps.append(eps)
    bad = 0; n = 0
    for s in itertools.combinations(range(N), 3):
        rest = [x for x in range(N) if x not in s]
        for q in itertools.combinations(rest, 4):
            vals = [c * chi[k1] * chi[k2] for c, k1, k2 in gp_terms(s, q)]; n += 1
            bad += vals[0] == vals[1] == vals[2]
    msg = []
    for eps in ok_eps:
        status, val, _ = propagate(F, o, eps)
        wrong = sum(1 for k, v in val.items() if chi[k] != v) if status == 'NO CONFLICT' else None
        msg.append(f"eps={eps:+d}: propagator {status}, derived {len(val) if wrong is not None else '-'} signs, "
                   f"{wrong} disagree with the true chirotope")
    print(f"[{name}] N={N} facets={len(F)}; unit conditions hold for eps in {ok_eps}; "
          f"GP relations violated by the true chirotope: {bad} of {n}; " + "; ".join(msg), flush=True)


if __name__ == "__main__":
    dC = load(sys.argv[1]); dD = load(sys.argv[2])
    control("Claim C polytope", dC["polytope_points"])
    control("conv(A) of Claim D", dD["A"])
    control("cyclic C(12,4)", [[t, t**2, t**3, t**4] for t in range(1, 13)])
    rng = random.Random(3)
    while True:
        P = [[rng.randrange(-50, 51) for _ in range(4)] for _ in range(40)]
        P = [p for p in P if True]
        # project to a sphere-ish set: keep points with large norm to get many vertices
        P = sorted(P, key=lambda p: -sum(x * x for x in p))[:11]
        try:
            fac, degen = hull_facets(P)
            if not degen and len({v for G, _ in fac for v in G}) == 11 and all(real_chirotope(P).values()):
                break
        except AssertionError:
            pass
    control("random 11-vertex polytope", P)
