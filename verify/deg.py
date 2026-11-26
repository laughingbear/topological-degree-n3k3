"""Topological degree of fewnomial maps via Bihan's Newton polytope formula (arXiv 2611.34336, Lemma 2.1 / Thm A).
deg F = degree of the PL map  boundary(P) -> S^{n-1},  a_s |-> C_s  (n = 3 here).
For a simplicial 3-polytope with outward-oriented facets (i,j,k), deg = (1/4pi) * sum of signed solid angles
of the cones (C_i, C_j, C_k).  Condition (H): for every face, the C's lie in an open half-space.
"""
import numpy as np
from scipy.spatial import ConvexHull

def hull_facets(pts):
    """Outward-oriented facets (index triples) of a simplicial 3-polytope."""
    h = ConvexHull(pts)
    c = pts.mean(0)
    F = []
    for s in h.simplices:
        i, j, k = s
        nrm = np.cross(pts[j] - pts[i], pts[k] - pts[i])
        if np.dot(nrm, pts[i] - c) < 0: j, k = k, j
        F.append((i, j, k))
    return F, sorted(set(h.vertices))

def solid(a, b, c):
    a, b, c = (v / np.linalg.norm(v) for v in (a, b, c))
    det = np.dot(a, np.cross(b, c))
    den = 1 + a @ b + b @ c + c @ a
    return 2 * np.arctan2(det, den)

def H_ok(F, C, tol=1e-9):
    """(H) on facets (implies edges/vertices): 0 not in the closed convex cone hull <=> exists w with w.C>0."""
    for (i, j, k) in F:
        M = np.array([C[i], C[j], C[k]])
        if abs(np.linalg.det(M)) < tol: return False
        # cone of 3 independent vectors is pointed: always lies in open half-space. Fine.
    # edges: not antiparallel; guaranteed by independence of facet triples.
    return True

def degree(F, C):
    s = sum(solid(C[i], C[j], C[k]) for (i, j, k) in F)
    d = s / (4 * np.pi)
    return d

def degree_thmA(F, C, A, u=None, rng=np.random.default_rng(1)):
    """Independent check: Bihan Theorem A count with a random u.  sign conventions -> compare |deg|."""
    if u is None: u = rng.normal(size=3)
    tot = 0
    for (i, j, k) in F:
        M = np.array([C[i], C[j], C[k]]).T
        lam = np.linalg.solve(M, u)
        if (lam > 0).all():
            tot += np.sign(np.linalg.det(M))   # facet outward-oriented => eps_Gamma constant sign
    return tot
