"""Referee's own root search. Two independent routes:
 (1) primal: Powell hybrid (scipy root) in log coordinates from many wide random starts, polished by 60-digit Newton (own loop);
 (2) Gale side: zeros in the polyhedron Delta = {y: D(y,1) > 0} of psi_j = sum_i B_ij log p_i(y), starts sampled inside Delta.
Reports sign det(dF/dx) at each root.  No completeness claim."""
import json, sys, numpy as np, mpmath as mp, sympy as sp
from scipy.optimize import root
fam = json.load(open(sys.argv[1])); A = np.array(fam['A'], float); C = np.array(fam['C_by_degree'][sys.argv[2]] if 'C_by_degree' in fam else fam['C'], float)
nst = int(sys.argv[3]); N, n = A.shape
mp.mp.dps = 60
Am = mp.matrix(A.tolist()); Cm = mp.matrix(C.tolist())
def Fy(y):
    e = A @ y; w = np.exp(e - e.max()); return C @ w
def Jy(y):
    e = A @ y; k = e.argmax(); w = np.exp(e - e[k]); return (C * w) @ (A - A[k])
def polish(y):
    y = mp.matrix(list(map(mp.mpf, y)))
    for it in range(60):
        w = [mp.e ** sum(Am[i, j] * y[j] for j in range(n)) for i in range(N)]
        F = mp.matrix([sum(Cm[r, i] * w[i] for i in range(N)) for r in range(n)])
        J = mp.matrix(n, n)
        for r in range(n):
            for j in range(n): J[r, j] = sum(Cm[r, i] * w[i] * Am[i, j] for i in range(N))
        scale = max(abs(Cm[r, i] * w[i]) for r in range(n) for i in range(N))
        if max(abs(v) for v in F) / scale < mp.mpf(10) ** -50: return y, mp.det(J), max(abs(v) for v in F) / scale
        y = y - mp.lu_solve(J, F)
    return None
rng = np.random.default_rng(12345); found = []
def add(y0, tag):
    r = polish(y0)
    if r is None: return
    y, dJ, res = r
    yf = np.array([float(v) for v in y])
    if any(np.abs(yf - z[0]).max() < 1e-8 for z in found): return
    found.append((yf, int(mp.sign(dJ)), float(res), tag))
# route 1
for t in range(nst):
    y0 = rng.normal(size=n) * rng.choice([0.5, 2, 5, 10])
    s = root(Fy, y0, jac=Jy, method='hybr', tol=1e-13)
    if s.success and np.abs(Fy(s.x)).max() < 1e-8: add(s.x, 'primal')
n1 = len(found)
# route 2: Gale dual
Ahat = sp.Matrix([[1] * N] + [[int(A[i, r]) for i in range(N)] for r in range(n)])
B = np.array(sp.Matrix.hstack(*Ahat.nullspace()).tolist(), float)           # N x k
D = np.array(sp.Matrix.hstack(*sp.Matrix(C.astype(int).tolist()).nullspace()).tolist(), float)   # N x (k+1)
k = B.shape[1]
def psi(v):                       # v in R^{k+1}, homogeneous; add normalisation sum(Dv)=1
    p = D @ v
    if (p <= 0).any(): return np.full(k + 1, 1e6)
    return np.append(B.T @ np.log(p), p.sum() - 1)
def dpsi(v):
    p = D @ v
    if (p <= 0).any(): return np.eye(k + 1)
    return np.vstack([B.T @ (D / p[:, None]), D.sum(0)])
# sample interior points of the cone K = {Dv >= 0}: random positive combos of feasible points found by rejection + hit and run
from scipy.optimize import linprog
lp = linprog(np.zeros(k + 2 - 1 + 1)[:k + 1].tolist() + [-1], A_ub=np.hstack([-D, np.ones((N, 1))]), b_ub=np.zeros(N), A_eq=[D.sum(0).tolist() + [0]], b_eq=[1], bounds=[(None, None)] * (k + 1) + [(None, 1)])
ngale = 0
if lp.status == 0 and lp.x[-1] > 1e-12:
    v = lp.x[:k + 1]
    for t in range(nst):
        for _ in range(8):       # hit-and-run inside {Dv>0, sum Dv = 1}
            d = rng.normal(size=k + 1); d -= D.sum(0) * (D.sum(0) @ d) / (D.sum(0) @ D.sum(0))
            p, q = D @ v, D @ d; lo = max([-p[i] / q[i] for i in range(N) if q[i] > 0], default=-1e9); hi = min([-p[i] / q[i] for i in range(N) if q[i] < 0], default=1e9)
            a = rng.uniform(lo, hi) if rng.random() < .5 else (lo + (hi - lo) * rng.beta(.3, .3)); v = v + a * d
        s = root(psi, v, jac=dpsi, method='hybr', tol=1e-13)
        if s.success and np.abs(psi(s.x)).max() < 1e-9 and (D @ s.x > 0).all():
            w = D @ s.x                         # proportional to x^{a_i}; solve for log x by least squares
            y = np.linalg.lstsq(np.hstack([np.ones((N, 1)), A]), np.log(w), rcond=None)[0][1:]
            ngale += 1; add(y, 'gale')
else: print('Delta empty (no positive vector in ker C): no positive roots')
tot = 0
for yf, s, res, tag in sorted(found, key=lambda z: tuple(z[0])):
    tot += s; print('x =', np.array2string(np.exp(yf), precision=6), 'sign det dF/dx =', '+' if s > 0 else '-', 'rel.residual %.1e' % res, tag)
print(f'label {sys.argv[2]}: roots found {len(found)} (primal route {n1}, Gale route converged {ngale} times, new from Gale {len(found) - n1}); sum of Jacobian signs = {tot}')
