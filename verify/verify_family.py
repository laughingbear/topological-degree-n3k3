"""Exact check of every system in data/sharp_family.json: (U_A), (U_C), and |deg F| = floor((3N-8)/4)
via Bihan's Theorem A (verify_exact.check), for several generic u."""
import json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from verify_exact import check
F = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'data/sharp_family.json'))
for N in sorted(F, key=int):
    A, C = F[N]['A'], F[N]['C']; n = int(N); target = (3 * n - 8) // 4; degs = set(); ua = uc = True
    for u in ([7, 11, -13], [1, 2, 3], [-5, 3, 8]):
        try: UA, UC, fac, d, _ = check(A, C, u)
        except ValueError: continue
        ua &= UA; uc &= UC; degs.add(d)
    ok = ua and uc and len(degs) == 1 and abs(degs.pop()) == target
    print(f'N={n:2d}  bound={target:2d}  (U_A)={ua} (U_C)={uc}  deg attains bound: {ok}')
