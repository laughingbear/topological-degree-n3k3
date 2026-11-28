r"""Prints, for a proof file whose steps are all facet conditions (V = 12, 16, 20 cases), a
self-contained hand-checkable version: the dual paths that fix the relative orientation of the
facets involved, and the violated Grassmann-Pluecker relation."""
import json, sys
from collections import deque
P = json.load(open(sys.argv[1])); F = [tuple(f) for f in P["facets"]]
used = []
for key, v, r in P["steps"]:
    assert r[0] == 'U'
    if tuple(r[1]) not in used: used.append(tuple(r[1]))
base = used[0]
def nbrs(t):
    for t2 in F:
        c = set(t) & set(t2)
        if len(c) == 3:
            i = [k for k in range(4) if t[k] not in c][0]; j = [k for k in range(4) if t2[k] not in c][0]
            yield t2, i, j
print(f"Facets used: {used}.  Normalise o{list(base)} = +1.")
print("Rule: if t' is obtained from t by deleting position i of t and position j of t' (0-based), o_t' = -(-1)^(i+j) o_t.")
for target in used[1:]:
    prev = {base: None}; Q = deque([base])
    while Q:
        t = Q.popleft()
        if t == target: break
        for t2, i, j in nbrs(t):
            if t2 not in prev: prev[t2] = (t, i, j); Q.append(t2)
    path = []; t = target
    while prev[t]: path.append((t, prev[t])); t = prev[t][0]
    o = 1; txt = f"  {list(base)}(+)"
    for t, (p, i, j) in reversed(path):
        o = -((-1) ** (i + j)) * o; txt += f" -[i={i},j={j}]-> {list(t)}({'+' if o > 0 else '-'})"
    print(txt)
print("Facet conditions (eps = %+d):" % P["eps"])
for key, v, r in P["steps"]:
    print(f"  chi{key} = {v:+d}   from facet {r[1]} and vertex {r[2]}")
print("Violated relation: s =", P["conflict"][1], " q =", P["conflict"][2])
