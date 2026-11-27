"""Positive roots of F(x) = sum_i C_i x^{a_i}: starting points from all (n+1)-monomial subsystems with a positive solution
(each has exactly one, by linear algebra in log coordinates), plus random starts; damped Newton in log coordinates with mpmath."""
import json, sys, itertools, numpy as np, mpmath as mp
S = json.load(open(sys.argv[1])); A = np.array(S['A'], float); C = np.array(S['C'], float); N = len(A); mp.mp.dps = 50
Am = [[mp.mpf(int(v)) for v in r] for r in S['A']]; Cm = [[mp.mpf(int(v)) for v in r] for r in S['C']]
def FJ(y):
    e = [sum(Am[i][j] * y[j] for j in range(4)) for i in range(N)]; m = max(e); w = [mp.e ** (v - m) for v in e]
    F = mp.matrix([sum(Cm[r][i] * w[i] for i in range(N)) for r in range(4)])
    J = mp.matrix(4, 4)
    for r in range(4):
        for j in range(4): J[r, j] = sum(Cm[r][i] * w[i] * Am[i][j] for i in range(N))
    return F, J
def newton(y0, iters=80):
    y = mp.matrix([mp.mpf(float(v)) for v in y0]); F, J = FJ(y); nf = mp.norm(F)
    for it in range(iters):
        try: d = mp.lu_solve(J, -F)
        except ZeroDivisionError: return None
        t = mp.mpf(1)
        while t > mp.mpf(10) ** -8:
            y2 = y + t * d; F2, J2 = FJ(y2); n2 = mp.norm(F2)
            if n2 < nf: break
            t /= 2
        else: return None
        y, F, J, nf = y2, F2, J2, n2
        if nf < mp.mpf(10) ** -35: return [float(v) for v in y]
    return None
starts = []
for T in itertools.combinations(range(N), 5):
    CT = C[:, T]; _, _, vt = np.linalg.svd(CT); k = vt[-1]
    if (k > 0).all() or (k < 0).all():
        M = np.concatenate([A[list(T)], np.ones((5, 1))], 1); sol = np.linalg.solve(M, np.log(np.abs(k))); starts.append(sol[:4])
rng = np.random.default_rng(0); base = list(starts)
for s in base:
    for _ in range(int(sys.argv[2])): starts.append(s + rng.normal(size=4) * rng.choice([0.05, 0.3, 1]))
roots = []
for y0 in starts:
    y = newton(y0)
    if y is None or max(abs(v) for v in y) > 80: continue
    if not any(max(abs(a - b) for a, b in zip(y, z)) < 1e-8 for z in roots): roots.append(y)
tot = 0
for y in sorted(roots):
    y = np.array(y); w = np.exp(A @ y - (A @ y).max()); s = int(np.sign(np.linalg.det((C * w) @ A))); tot += s
    print('log x =', np.array2string(y, precision=4), 'sign det J =', '+' if s > 0 else '-')
print('simplex subsystems with a positive root:', len(base), ' starts:', len(starts), ' roots found:', len(roots), ' sum of Jacobian signs:', tot)
