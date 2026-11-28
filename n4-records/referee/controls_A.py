r"""Controls for the degree code (is a degree 'right for the wrong reason'?).
(1) random integer vertex maps on each certified sphere: the signed count must not depend on u;
(2) the same count with ONE facet's orientation sign flipped must depend on u for the claimed
    maps (so the u-independence observed in check_A.py really tests the orientation);
(3) summing |s_t| over cones (unsigned cover count) vs signed count is reported in check_A.py."""
import random, sys, os
from common import *
rng = random.Random(5)
for path in sys.argv[1:]:
    d = load(path); F = [tuple(sorted(f)) for f in d["facets"]]; o = orientation(F); V = 1 + max(max(f) for f in F)
    us = [[rng.randrange(-10**6, 10**6) for _ in range(4)] for _ in range(40)]
    res = []
    for trial in range(5):
        C = [[rng.randrange(-99, 100) for _ in range(4)] for _ in range(V)]
        try:
            s = facet_signs(F, o, C); vals = {degree(F, s, C, u) for u in us}
        except (AssertionError, ValueError):
            continue
        res.append(sorted(vals))
    C = d["C_int"]; s = facet_signs(F, o, C)
    broken = 0
    for k in range(len(F)):
        s2 = list(s); s2[k] = -s2[k]
        broken += len({degree(F, s2, C, u) for u in us}) > 1
    print(f"{os.path.basename(path)}: random maps -> degree sets per map {res}; "
          f"flipping one facet sign breaks u-independence for {broken} of {len(F)} facets")
