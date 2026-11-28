"""Claims A1-A6 (ii),(iii) and the degree part of C: orientation recomputed, exact degree."""
import random, sys
from common import *

for path in sys.argv[1:]:
    d = load(path)
    facets = [tuple(sorted(f)) for f in d["facets"]]
    C = d["C_int"]
    f = check_closed_manifold(facets)
    o = orientation(facets)
    stored = d.get("orientation")
    rel = None
    if stored:
        rels = {a * b for a, b in zip(o, stored)}
        rel = "equal" if rels == {1} else "negated" if rels == {-1} else "DIFFERENT (not ± ours)"
    s = facet_signs(facets, o, C)
    rng = random.Random(20261127)
    us = [[1, 2, 3, 5], [-3, 1, 7, 2], [5, -2, -1, 9], [-7, -11, 13, -17]]
    us += [[rng.randrange(-10**9, 10**9) for _ in range(4)] for _ in range(60)]
    vals = {}; ncover = []
    for u in us:
        try:
            dg = degree(facets, s, C, u)
        except ValueError:
            continue
        vals[dg] = vals.get(dg, 0) + 1
        ncover.append(sum(in_open_cone([C[v] for v in t], u) for t in facets))
    assert len(vals) == 1, f"degree depends on u: {vals}"
    deg = next(iter(vals))
    q = sum(1 for x in s if x == -sgn(deg))
    V, E, T, F = f
    lhs, rhs2 = 8 * abs(deg), 4 * F - E - q   # 4|deg| < 2F - E/2 - q/2  <=>  8|deg| < 4F - E - q
    cyc = (V - 1) * (V - 3) / 8 if V % 2 else (V - 2) * (V - 4) / 8
    print(f"{path.split('/')[-1]}: f={f} (file f={tuple(d['f'])}) deg={deg} (file {d['degree']}) "
          f"over {sum(vals.values())} generic u; cones covering u: min {min(ncover)} max {max(ncover)}; "
          f"q={q} (file {d.get('q')}); orientation vs stored: {rel}; "
          f"Thm M: {abs(deg)*4} < {rhs2/2} -> {lhs < rhs2}; cyclic {cyc}; max|C|={max(abs(x) for r in C for x in r)}")
