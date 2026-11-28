r"""3-sphere recognition by SHELLING (not bistellar flips).

Fact used: a shellable simplicial complex in which every codimension-1 face lies in exactly
two facets is a PL sphere (Danaraj-Klee 1974; Bjorner, 'Posets, regular CW complexes and
Bruhat order', Prop. 4.3 / Ziegler, Lectures on Polytopes, Ch. 8).

find_shelling: randomised greedy + restarts.  check_shelling: independent, direct from the
definition: for k >= 2, F_k ∩ (F_1 ∪ ... ∪ F_{k-1}) is a non-empty union of triangles of F_k,
i.e. every intersection F_j ∩ F_k (j<k) is contained in a triangle F_k \ {v} that itself lies
in some earlier facet.
Also: integer homology via Smith normal form (independent second evidence; necessary only).
"""
import itertools, random, sys
from common import load, check_closed_manifold


def check_shelling(order):
    order = [frozenset(f) for f in order]
    for k in range(1, len(order)):
        Fk = order[k]
        prev = order[:k]
        good = [Fk - {v} for v in Fk if any((Fk - {v}) <= Fj for Fj in prev)]
        if not good:
            return False
        for Fj in prev:
            I = Fj & Fk
            if not any(I <= g for g in good):
                return False
    return True


def find_shelling(facets, tries=2000, seed=0):
    rng = random.Random(seed)
    F = [frozenset(f) for f in facets]
    for _ in range(tries):
        order = [rng.choice(F)]
        rest = set(F) - {order[0]}
        # incremental: track faces present
        present = set()
        def addfaces(f):
            for r in range(1, 5):
                for s in itertools.combinations(sorted(f), r):
                    present.add(s)
        addfaces(order[0])
        ok = True
        while rest:
            cands = []
            for f in rest:
                fs = sorted(f)
                good = [tuple(x for x in fs if x != v) for v in fs]
                good = [g for g in good if g in present]
                if not good:
                    continue
                # every face of f that is present must lie in a good triangle
                fine = True
                for r in range(1, 4):
                    for s in itertools.combinations(fs, r):
                        if s in present and not any(set(s) <= set(g) for g in good):
                            fine = False; break
                    if not fine:
                        break
                if fine:
                    cands.append((len(good), f))
            if not cands:
                ok = False; break
            m = max(c[0] for c in cands)   # prefer facets glued along most triangles
            f = rng.choice([c[1] for c in cands if c[0] == m])
            order.append(f); rest.remove(f); addfaces(f)
        if ok:
            return [tuple(sorted(f)) for f in order]
    return None


# ---------- integer homology (Smith normal form over Z, exact) ----------
def boundary_matrix(rows_faces, cols_faces):
    idx = {f: i for i, f in enumerate(rows_faces)}
    M = [[0] * len(cols_faces) for _ in rows_faces]
    for j, c in enumerate(cols_faces):
        for i in range(len(c)):
            M[idx[c[:i] + c[i + 1:]]][j] = (-1) ** i
    return M


def smith_diag(M):
    """Nonzero invariant factors of an integer matrix (exact)."""
    M = [r[:] for r in M]
    m = len(M); n = len(M[0]) if M else 0
    diag = []; r0 = 0; c0 = 0
    while r0 < m and c0 < n:
        piv = None
        best = None
        for i in range(r0, m):
            for j in range(c0, n):
                if M[i][j] != 0 and (best is None or abs(M[i][j]) < best):
                    best = abs(M[i][j]); piv = (i, j)
        if piv is None:
            break
        i, j = piv
        M[r0], M[i] = M[i], M[r0]
        for row in M:
            row[c0], row[j] = row[j], row[c0]
        while True:
            done = True
            p = M[r0][c0]
            for i in range(r0 + 1, m):
                if M[i][c0] != 0:
                    q = M[i][c0] // p
                    if q:
                        M[i] = [a - q * b for a, b in zip(M[i], M[r0])]
                    if M[i][c0] != 0:
                        M[r0], M[i] = M[i], M[r0]; p = M[r0][c0]; done = False
            for j in range(c0 + 1, n):
                if M[r0][j] != 0:
                    q = M[r0][j] // p
                    if q:
                        for row in M:
                            row[j] -= q * row[c0]
                    if M[r0][j] != 0:
                        for row in M:
                            row[c0], row[j] = row[j], row[c0]
                        p = M[r0][c0]; done = False
            if done:
                # ensure divisibility of the remainder is not needed for torsion detection:
                # we only report |p|; torsion-free iff all |p| == 1 after full reduction
                break
        diag.append(abs(M[r0][c0])); r0 += 1; c0 += 1
    return diag


def homology(facets):
    facets = sorted(tuple(sorted(f)) for f in facets)
    sk = {k: sorted({s for f in facets for s in itertools.combinations(f, k)}) for k in (1, 2, 3, 4)}
    res = {}
    ranks = {}; tors = {}
    for k in (2, 3, 4):  # boundary from k-sets to (k-1)-sets
        d = smith_diag(boundary_matrix(sk[k - 1], sk[k]))
        ranks[k] = len(d); tors[k] = [x for x in d if x != 1]
    # H_i (i = dim) : dim C_i = len(sk[i+1]); rank d_i = ranks[i+1]
    b0 = len(sk[1]) - ranks[2]
    b1 = len(sk[2]) - ranks[2] - ranks[3]
    b2 = len(sk[3]) - ranks[3] - ranks[4]
    b3 = len(sk[4]) - ranks[4]
    return (b0, b1, b2, b3), (tors[2], tors[3], tors[4])


if __name__ == "__main__":
    for path in sys.argv[1:]:
        d = load(path)
        facets = [tuple(sorted(f)) for f in d["facets"]]
        f = check_closed_manifold(facets)
        betti, tors = homology(facets)
        sh = find_shelling(facets)
        ok = sh is not None and sorted(sh) == sorted(facets) and check_shelling(sh)
        print(path.split("/")[-1], "f =", f, "| closed 3-manifold: yes | Betti (Z) =", betti,
              "torsion(H0,H1,H2) =", tors, "| shelling found+verified:", ok)
        if ok:
            import json, os
            os.makedirs("out", exist_ok=True)
            json.dump(sh, open("out/shelling_" + path.split("/")[-1], "w"))
