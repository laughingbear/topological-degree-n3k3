r"""E1 side check: which of the stated symmetries do the certificate complexes really have?
Two blocks a_i = i, b_i = n+i (V = 2n).  rot: i -> i+1 in both blocks.
refl_c: a_i -> a_{-i}, b_i -> b_{c-i}."""
import sys, os
from common import load
for path in sys.argv[1:]:
    d = load(path); F = {tuple(sorted(f)) for f in d["facets"]}; V = 1 + max(max(f) for f in F)
    out = [os.path.basename(path), d.get("note", "")[:70]]
    if V % 2 == 0:
        n = V // 2
        rot = lambda v: (v // n) * n + (v % n + 1) % n
        out.append(f"rot-invariant: {({tuple(sorted(map(rot, f))) for f in F} == F)}")
        cs = []
        for c in range(n):
            ref = lambda v: ((-v) % n) if v < n else n + ((c - (v - n)) % n)
            if {tuple(sorted(map(ref, f))) for f in F} == F: cs.append(c)
        out.append(f"reflection-invariant for c in {cs}")
    # full automorphism-ish info: vertex degree sequence
    deg = [sum(v in f for f in F) for v in range(V)]
    out.append(f"facets per vertex {sorted(set(deg))}")
    print(" | ".join(out))
