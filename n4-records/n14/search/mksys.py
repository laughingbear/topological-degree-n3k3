"""Turn a BEST block of sympoly output into an explicit system JSON (integer exponents A, integer coefficients C).
usage: mksys.py out.txt scaleA scaleC result.json [jitter]
(jitter: relative size of a random perturbation of the vectors before rounding; needed when the exactly symmetric
vectors have vanishing 4x4 minors, which (U_C) forbids)"""
import sys, json, math, re
def orbit(base, n, p, q):
    P = []
    for b in base:
        for i in range(n):
            t1 = 2 * math.pi * p * i / n; t2 = 2 * math.pi * q * i / n
            P.append([math.cos(t1) * b[0] - math.sin(t1) * b[1], math.sin(t1) * b[0] + math.cos(t1) * b[1],
                      math.cos(t2) * b[2] - math.sin(t2) * b[3], math.sin(t2) * b[2] + math.cos(t2) * b[3]])
    return P
L = open(sys.argv[1]).read().split("\n"); i = max(k for k, l in enumerate(L) if l.startswith("BEST"))
h = L[i].split(); kv = dict(zip(h[1::2], map(int, h[2::2]))); n, m = kv['n'], kv['m']
X = [list(map(float, L[i + 1 + j].split()[1:])) for j in range(m)]; Y = [list(map(float, L[i + 1 + m + j].split()[1:])) for j in range(m)]
A = orbit(X, n, kv['p1'], kv['q1']); C = orbit(Y, n, kv['p2'], kv['q2'])
sa = float(sys.argv[2]); sc = float(sys.argv[3])
if len(sys.argv) > 5:
    import random; random.seed(1); jit = float(sys.argv[5]); mc0 = max(abs(x) for r in C for x in r)
    C = [[x + jit * mc0 * random.gauss(0, 1) for x in r] for r in C]
    ma0 = max(abs(x) for r in A for x in r); A = [[x + jit * ma0 * random.gauss(0, 1) for x in r] for r in A]
mx = max(abs(x) for r in A for x in r); Ai = [[round(x / mx * sa) for x in r] for r in A]
mn = [min(r[k] for r in Ai) for k in range(4)]; Ai = [[r[k] - mn[k] for k in range(4)] for r in Ai]
mc = max(abs(x) for r in C for x in r); Ci = [[round(C[v][k] / mc * sc) for v in range(len(C))] for k in range(4)]
json.dump(dict(N=n * m, n=4, k=n * m - 5, deg=kv['deg'], A=Ai, C=Ci, source=dict(kv, X=X, Y=Y)), open(sys.argv[4], "w"))
print("wrote", sys.argv[4], "N", n * m, "claimed", kv['deg'], "max exponent", max(max(r) for r in Ai), "max |coef|", max(abs(x) for r in Ci for x in r))
