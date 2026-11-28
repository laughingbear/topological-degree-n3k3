"""Exact convex hull facets of points in R^4 (referee code)."""
import itertools
from common import det, sgn


def hom(P):
    return [[1] + list(p) for p in P]


def hull_facets(P):
    """All 4-subsets G of the points such that every other point is strictly on one side of
    aff(G).  Also returns the list of 'degenerate supports' (G affinely independent, all other
    points weakly on one side with at least one ON the hyperplane) - must be empty for the
    boundary complex to be exactly the returned list.  Returns [(G, side)], side = sign of
    det[(1,p_a),(1,p_b),(1,p_c),(1,p_d),(1,p_w)] (rows), w off G, G increasing."""
    H = hom(P); N = len(P); out = []; degen = []
    for G in itertools.combinations(range(N), 4):
        sides = [sgn(det([H[i] for i in G] + [H[w]])) for w in range(N) if w not in G]
        if all(x == 0 for x in sides):
            continue
        nz = {x for x in sides if x}
        if len(nz) == 1:
            if 0 in sides:
                degen.append(G)
            else:
                out.append((G, nz.pop()))
    return out, degen


