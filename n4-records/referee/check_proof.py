r"""Stand-alone checker for the proofs written by chiro.py.  Shares NO code with chiro.py
(own parity, own orientation, own GP expansion), so that an encoder bug in chiro.py cannot
validate itself.  Usage: check_proof.py out/proof_*.json  [--print]

A proof is a list of steps (key, value, reason) on sorted 5-tuples plus a final conflict.
  reason ['U', t, w]         : facet t, vertex w not in t:  chi(t_0,t_1,t_2,t_3,w) = eps*o_t
  reason ['GP', s, q, i]     : in the GP relation (s; q) the two terms other than i are already
                               known and equal, the partner factor in term i is already known
  conflict ['GP', s, q]      : all three terms known and equal.
"""
import json, sys


def perm_sign(t):
    """sign of the permutation sorting t, by cycle-free bubble count (independent implementation)."""
    t = list(t); sg = 1
    for i in range(len(t)):
        for j in range(len(t) - 1 - i):
            if t[j] > t[j + 1]:
                t[j], t[j + 1] = t[j + 1], t[j]; sg = -sg
    return tuple(t), sg


def orient(facets):
    """o_t with o_t(-1)^i = -o_t'(-1)^j across every triangle; o=+1 on first facet."""
    o = [0] * len(facets); o[0] = 1; changed = True
    while changed:
        changed = False
        for x, t in enumerate(facets):
            if not o[x]:
                continue
            for y, t2 in enumerate(facets):
                common = set(t) & set(t2)
                if len(common) != 3:
                    continue
                i = [k for k in range(4) if t[k] not in common][0]
                j = [k for k in range(4) if t2[k] not in common][0]
                want = -o[x] * (-1) ** i * (-1) ** j
                if o[y] == 0:
                    o[y] = want; changed = True
                else:
                    assert o[y] == want, "not coherently orientable"
    assert all(o)
    return o


def terms(s, q):
    a, b, c, d = q
    assert a < b < c < d and list(s) == sorted(s) and not set(s) & set(q) and len(set(s)) == 3
    res = []
    for coef, P, R in ((+1, (a, b), (c, d)), (-1, (a, c), (b, d)), (+1, (a, d), (b, c))):
        k1, s1 = perm_sign(tuple(s) + P); k2, s2 = perm_sign(tuple(s) + R)
        res.append((coef * s1 * s2, k1, k2))
    return res


def check(path, show=False):
    P = json.load(open(path))
    facets = [tuple(f) for f in P["facets"]]; eps = P["eps"]
    assert all(list(f) == sorted(f) for f in facets)
    o = orient(facets); fidx = {f: i for i, f in enumerate(facets)}
    chi = {}
    lines = []
    for key, v, r in P["steps"]:
        key = tuple(key)
        assert key not in chi and v in (1, -1)
        if r[0] == 'U':
            t, w = tuple(r[1]), r[2]
            assert t in fidx and w not in t
            k, sg = perm_sign(t + (w,))
            assert k == key and v == sg * eps * o[fidx[t]], "bad facet condition"
            lines.append(f"chi{list(key)} = {v:+d}   [facet {list(t)} (o={o[fidx[t]]:+d}), vertex {w} on its inner side]")
        else:
            s, q, i = tuple(r[1]), tuple(r[2]), r[3]
            T = terms(s, q)
            others = [T[j] for j in range(3) if j != i]
            vals = [c * chi[k1] * chi[k2] for c, k1, k2 in others]   # KeyError if not yet derived
            assert vals[0] == vals[1], "GP step: other two terms differ"
            c, k1, k2 = T[i]
            assert key in (k1, k2)
            partner = k2 if key == k1 else k1
            assert partner != key
            assert c * v * chi[partner] == -vals[0], "GP step: wrong forced value"
            lines.append(f"chi{list(key)} = {v:+d}   [GP s={list(s)}, q={list(q)}: other two terms both {vals[0]:+d}; partner chi{list(partner)}={chi[partner]:+d}]")
        chi[key] = v
    cf = P["conflict"]
    if cf[0] == 'GP':
        T = terms(tuple(cf[1]), tuple(cf[2]))
        vals = [c * chi[k1] * chi[k2] for c, k1, k2 in T]
        assert vals[0] == vals[1] == vals[2], "final relation is not violated"
        lines.append(f"CONTRADICTION in GP s={cf[1]}, q={cf[2]}: the three terms "
                     + ", ".join(f"({c:+d})*chi{list(k1)}*chi{list(k2)}" for c, k1, k2 in T)
                     + f" all equal {vals[0]:+d}, but they must sum to zero with these signs.")
    else:
        raise AssertionError("unit-level conflict: inspect by hand")
    if show:
        print("\n".join(f"{n+1:3d}. {l}" for n, l in enumerate(lines)))
    return len(P["steps"])


if __name__ == "__main__":
    show = "--print" in sys.argv
    for p in sys.argv[1:]:
        if p.startswith("--"):
            continue
        n = check(p, show)
        print(f"{p}: proof VALID ({n} steps)")
