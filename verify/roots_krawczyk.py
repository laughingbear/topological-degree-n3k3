"""Independent root solve for F(x)=sum_i C_i x^{a_i} on R^3_{>0}, plus Krawczyk certification.

Formulation: z = log x, Fhat(z) = C @ softmax(A z)  (= F(x)/sum_j x^{a_j}, same positive zeros).
Starts: (i) uniform z in a big box, (ii) moment-map sampling: random interior point q of P=conv(A),
solve grad logsumexp(Az) = q (convex), giving starts spread over the whole compactified orthant.
Certification: interval Krawczyk operator (mpmath.iv) on G(z) = C @ exp(A z - m) with a fixed shift m.
Jacobian sign: sign det DF(x) = sign det(C diag(x^a) A^T)  (diag(1/x) > 0 dropped).
"""
import sys, json
import numpy as np
from scipy.optimize import least_squares, minimize
import mpmath
from mpmath import iv, mpf

S = json.load(open(sys.argv[1])); A = np.array(S['A'], float); C = np.array(S['C'], float)
nstart = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
N, n = A.shape
def soft(z):
    m = A @ z; m -= m.max(); e = np.exp(m); return e / e.sum()
def Fh(z): return C @ soft(z)
def Jh(z):
    p = soft(z); Ab = p @ A
    return C @ (p[:, None] * (A - Ab))
rng = np.random.default_rng(12345)
roots = []
def add(z):
    for r in roots:
        if np.linalg.norm(r - z) < 1e-7: return False
    roots.append(z); return True
def polish(z0):
    r = least_squares(Fh, z0, jac=Jh, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
    z = r.x
    for _ in range(8):  # Newton polish
        try: z = z - np.linalg.solve(Jh(z), Fh(z))
        except np.linalg.LinAlgError: return None
    if np.linalg.norm(Fh(z)) < 1e-13 and np.all(np.isfinite(z)) and np.abs(z).max() < 60: return z
    return None
# vertices of P for moment-map sampling
from scipy.spatial import ConvexHull
H = ConvexHull(A); Vx = A[H.vertices]
hits_per_start = []
for t in range(nstart):
    if t % 2 == 0:
        z0 = rng.uniform(-6, 6, n)
    else:
        w = rng.dirichlet(np.full(len(Vx), 0.3)); q = w @ Vx
        f = lambda z: np.log(np.exp(A @ z - (A @ z).max()).sum()) + (A @ z).max() - q @ z
        z0 = minimize(f, np.zeros(n), method='BFGS').x
    z = polish(z0)
    if z is not None and add(z): print(f'  start {t}: new root #{len(roots)}', flush=True)
print(len(roots), 'distinct positive zeros found from', nstart, 'starts')

# Krawczyk certification
mpmath.mp.prec = iv.prec = 120
def G_iv(Z, m):
    E = [iv.exp(sum(int(A[i, r]) * Z[r] for r in range(n)) - m) for i in range(N)]
    return [sum(int(C[j, i]) * E[i] for i in range(N)) for j in range(n)], E
def certify(z, rad=1e-9):
    m = iv.mpf(float((A @ z).max()))
    Y_num = np.linalg.inv(C @ (np.exp(A @ z - (A @ z).max())[:, None] * A))  # inverse of dG at z (float)
    Y = [[iv.mpf(float(Y_num[a, b])) for b in range(n)] for a in range(n)]
    y = [mpf(float(v)) for v in z]
    X = [iv.mpf([y[r] - rad, y[r] + rad]) for r in range(n)]
    Gy, _ = G_iv([iv.mpf(v) for v in y], m)
    _, EX = G_iv(X, m)
    J = [[sum(int(C[j, i]) * EX[i] * int(A[i, r]) for i in range(N)) for r in range(n)] for j in range(n)]
    K = []
    for a in range(n):
        val = iv.mpf(y[a]) - sum(Y[a][b] * Gy[b] for b in range(n))
        for r in range(n):
            M = (1 if a == r else 0) - sum(Y[a][b] * J[b][r] for b in range(n))
            val += M * (X[r] - y[r])
        K.append(val)
    ok = all(K[r].a > X[r].a and K[r].b < X[r].b for r in range(n))
    return ok
tot = 0
for z in roots:
    x = np.exp(z)
    s = int(np.sign(np.linalg.det(C @ (np.exp(A @ z - (A @ z).max())[:, None] * A))))
    tot += s
    print('x =', np.array2string(x, precision=12), ' z =', np.array2string(z, precision=6), ' jac sign', s,
          ' Krawczyk unique root in box r=1e-9:', certify(z))
print('sum of Jacobian signs =', tot)
