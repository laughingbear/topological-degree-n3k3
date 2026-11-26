"""Independent numerical check: positive roots of F(x)=sum_i C_i x^{a_i} by multi-start Newton in log coordinates (mpmath, 40 digits);
sum of sign det Jacobian over the roots found.  (No completeness claim.)"""
import json, sys, numpy as np, mpmath as mp
from scipy.optimize import least_squares
S = json.load(open(sys.argv[1])); A = np.array(S['A'], float); C = np.array(S['C'], float) if 'C' in S else np.array(S['C_by_degree'][sys.argv[2]], float)
def G(y):
    e = A @ y; m = e.max(); w = np.exp(e - m); return C @ w
def J(y):
    e = A @ y; m = e.max(); w = np.exp(e - m); return (C * w) @ A - np.outer(C @ w, A[e.argmax()])
rng = np.random.default_rng(0); roots = []
mp.mp.dps = 40
Am = mp.matrix(A.tolist()); Cm = mp.matrix(C.tolist())
def Fm(*y):
    w = [mp.e ** sum(Am[i, j] * y[j] for j in range(4)) for i in range(9)]
    return [sum(Cm[r, i] * w[i] for i in range(9)) for r in range(4)]
for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 4000):
    y0 = rng.normal(size=4) * rng.choice([1, 3, 6])
    r = least_squares(G, y0, jac=J, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=200)
    if np.abs(r.fun).max() > 1e-9: continue
    y = r.x
    try:
        ym = mp.findroot(lambda a, b, c, d: Fm(a - 0, b, c, d), tuple(y), tol=1e-30, maxsteps=50)
    except Exception: continue
    y = np.array([float(v) for v in ym])
    if any(np.abs(y - z).max() < 1e-7 for z in roots): continue
    roots.append(y)
tot = 0
for y in sorted(roots, key=lambda z: tuple(z)):
    x = np.exp(y); w = np.exp(A @ y)
    Jx = (C * w) @ A / x          # dF/dx
    s = int(np.sign(np.linalg.det(Jx))); tot += s
    print('x =', np.array2string(x, precision=6), ' sign det J =', '+' if s > 0 else '-')
print('positive roots found:', len(roots), ' sum of Jacobian signs:', tot)
