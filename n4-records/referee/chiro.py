r"""Claim B, independent encoder + unit propagator with proof logging.

Necessary conditions for K = boundary of a simplicial 4-polytope P with vertices p_0..p_{V-1}
(after an arbitrarily small perturbation of the vertices, which does not change the boundary
complex of a simplicial polytope, no 5 points lie on a hyperplane):
  chi(i1..i5) = sign det[rows (1,p_i1),...,(1,p_i5)]  is alternating, never 0, and
  (U) for each facet t=(a<b<c<d) and each w not in t:   chi(a,b,c,d,w) = eps * o_t
      (all other vertices lie strictly on one side of aff(t); that these side signs form a
       coherent orientation is standard, and it is verified numerically on two real polytopes
       in check_CD.py; the coherent orientation of a connected pseudomanifold is unique up to sign)
  (GP) for each 3-set s and a<b<c<d outside s (three-term Grassmann-Pluecker identity
       [sab][scd] - [sac][sbd] + [sad][sbc] = 0 for real determinants):
       the three signs  chi(sab)chi(scd),  -chi(sac)chi(sbd),  chi(sad)chi(sbc)  are not all equal.
Only NECESSITY of (U),(GP) is used, so a contradiction is a proof of non-polytopality.
chi -> -chi maps solutions for eps to solutions for -eps (mirror image), so eps=+1 is w.l.o.g.;
we nevertheless run both.

Propagation rule (= unit propagation on the clausal form): if two of the three terms of a (GP)
are known and equal to v, the third term must be -v; if one of its two factors is known the other
is forced.  Conflict: all three terms known and equal.
"""
import itertools, json, sys
from collections import deque
from common import load, check_closed_manifold, orientation


def par_sorted(t):
    """(sorted tuple, sign of the sorting permutation)"""
    inv = sum(1 for i in range(len(t)) for j in range(i + 1, len(t)) if t[i] > t[j])
    return tuple(sorted(t)), (-1 if inv & 1 else 1)


def gp_terms(s, q):
    """s sorted 3-tuple, q=(a<b<c<d). Returns three terms (coef, key1, key2) meaning
    coef*chi[key1]*chi[key2] with keys sorted 5-tuples."""
    a, b, c, d = q
    out = []
    for coef, (x, y), (z, w) in ((1, (a, b), (c, d)), (-1, (a, c), (b, d)), (1, (a, d), (b, c))):
        k1, p1 = par_sorted(s + (x, y)); k2, p2 = par_sorted(s + (z, w))
        out.append((coef * p1 * p2, k1, k2))
    return out


def unit_conditions(facets, o, eps):
    V = 1 + max(max(f) for f in facets)
    units = []
    for t, ot in zip(facets, o):
        for w in range(V):
            if w in t:
                continue
            k, p = par_sorted(tuple(t) + (w,))
            units.append((k, p * eps * ot, (tuple(t), w)))
    return units


def propagate(facets, o, eps=1, log=None):
    """Returns ('CONFLICT', steps, conflict) or ('NO CONFLICT', val, None).
    steps: list of (key, value, reason); reason = ('U', t, w) or ('GP', s, q, index_of_forced_term)."""
    V = 1 + max(max(f) for f in facets)
    val = {}; steps = []; Q = deque()
    for k, v, (t, w) in unit_conditions(facets, o, eps):
        if k in val:
            if val[k] != v:
                return 'CONFLICT', steps, ('U', t, w, k, v)
            continue
        val[k] = v; steps.append((k, v, ('U', t, w))); Q.append(k)
    nrel = 0
    while Q:
        S = Q.popleft(); Sset = set(S)
        outside = [x for x in range(V) if x not in Sset]
        for s in itertools.combinations(S, 3):
            e, f = [x for x in S if x not in s]
            for g, h in itertools.combinations(outside, 2):
                q = tuple(sorted((e, f, g, h)))
                terms = gp_terms(s, q); nrel += 1
                tv = []
                for coef, k1, k2 in terms:
                    v1 = val.get(k1); v2 = val.get(k2)
                    tv.append(None if v1 is None or v2 is None else coef * v1 * v2)
                known = [x for x in tv if x is not None]
                if len(known) == 3:
                    if tv[0] == tv[1] == tv[2]:
                        return 'CONFLICT', steps, ('GP', s, q)
                elif len(known) == 2 and known[0] == known[1]:
                    i = tv.index(None); coef, k1, k2 = terms[i]
                    v1 = val.get(k1); v2 = val.get(k2)
                    if v1 is None and v2 is None:
                        continue
                    kk, other = (k1, v2) if v1 is None else (k2, v1)
                    newv = -known[0] * coef * other   # coef*newv*other = -known
                    val[kk] = newv; steps.append((kk, newv, ('GP', s, q, i))); Q.append(kk)
    return 'NO CONFLICT', val, None


def minimise(steps, conflict):
    """Keep only the steps the conflict depends on (backward closure)."""
    pos = {k: i for i, (k, v, r) in enumerate(steps)}
    need = set()
    def deps_of_rel(s, q, skip=None):
        ks = []
        for coef, k1, k2 in gp_terms(s, q):
            ks += [k1, k2]
        return [k for k in ks if k != skip]
    stack = list(deps_of_rel(conflict[1], conflict[2])) if conflict[0] == 'GP' else [conflict[3]]
    while stack:
        k = stack.pop()
        if k in need:
            continue
        need.add(k)
        r = steps[pos[k]][2]
        if r[0] == 'GP':
            stack += deps_of_rel(r[1], r[2], skip=k)
    return [st for st in steps if st[0] in need]


if __name__ == "__main__":
    import os, time
    os.makedirs("out", exist_ok=True)
    for path in sys.argv[1:]:
        d = load(path); facets = [tuple(sorted(f)) for f in d["facets"]]
        check_closed_manifold(facets); o = orientation(facets)
        name = os.path.basename(path)
        for eps in (1, -1):
            t0 = time.time()
            status, steps, conflict = propagate(facets, o, eps)
            if status == 'CONFLICT':
                proof = minimise(steps, conflict)
                json.dump({"facets": facets, "eps": eps, "steps": proof, "conflict": conflict},
                          open(f"out/proof_eps{eps:+d}_{name}", "w"))
                nu = sum(1 for p in proof if p[2][0] == 'U')
                print(f"{name} eps={eps:+d}: CONFLICT by unit propagation after {len(steps)} assignments; "
                      f"minimised proof: {nu} facet conditions + {len(proof)-nu} GP steps + final GP {conflict[1:]}"
                      f"  [{time.time()-t0:.0f}s]", flush=True)
            else:
                print(f"{name} eps={eps:+d}: no conflict by unit propagation ({len(steps)} of "
                      f"{len(list(itertools.combinations(range(1+max(max(f) for f in facets)),5)))} "
                      f"signs fixed)  [{time.time()-t0:.0f}s]", flush=True)
