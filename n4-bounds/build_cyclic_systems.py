"""Explicit polynomial systems for Theorem C.
Odd N = 2m+1: exponents a_i = (i, i^2, i^3, i^4) (moment curve, so the Newton polytope is the cyclic polytope C(N,4)),
coefficients C_i = (-1)^i * (cos t, sin t, cos 3t, sin 3t), t = pi*i/N, rounded to integers after scaling by 10^6.
Even N: the system for N-1 plus one monomial in the interior of the Newton polytope (the degree depends only on the vertices).
Every system is checked exactly with verify_exact.py (Bihan, arXiv:2611.34336, Thm 1.3): (U_A), (U_C), and |deg F| = m(m-1)/2."""
import json, math, random
from verify_exact import structure, degree, UC
random.seed(7); out = {}
for N in range(5, 15):
    V = N if N % 2 else N - 1; m = (V - 1) // 2
    A = [[i, i**2, i**3, i**4] for i in range(V)]
    C = [[round(1e6 * (-1)**i * f(k * math.pi * i / V)) for i in range(V)] for f, k in ((math.cos, 1), (math.sin, 1), (math.cos, 3), (math.sin, 3))]
    if N % 2 == 0:
        while True:
            w = [random.randint(1, 9) for _ in range(V)]; s = sum(w)
            p = [round(sum(w[i] * A[i][r] for i in range(V)) / s) for r in range(4)]
            c = [random.randint(-10**6, 10**6) for _ in range(4)]
            A2 = A + [p]; C2 = [C[r] + [c[r]] for r in range(4)]
            UA, facets = structure(A2)
            if UA and UC(C2) and all(V not in G for G, _ in facets): A, C = A2, C2; break
    UA, facets = structure(A); assert UA and UC(C)
    degs = set()
    for _ in range(12):
        u = [random.randint(-50, 50) for _ in range(4)]
        try: degs.add(degree(facets, C, u)[0])
        except ValueError: pass
    assert len(degs) == 1 and abs(next(iter(degs))) == m * (m - 1) // 2, (N, degs)
    out[str(N)] = {"N": N, "n": 4, "k": N - 5, "deg": next(iter(degs)), "A": A, "C": C}
    print(f"N = {N:2d}: (U_A), (U_C) hold, facets {len(facets)}, deg F = {next(iter(degs))}, |deg| = {m*(m-1)//2}", flush=True)
json.dump(out, open("cyclic_systems.json", "w"))
