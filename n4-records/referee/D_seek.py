r"""Claim D, step 1 (NUMERICAL, not a certificate): search for positive zeros of F in log
coordinates t = log x with multiprecision Newton.  Seeds: (a) for every 5-subset S of monomials
whose coefficient columns have a positive linear dependence lambda, the point t with
a_i.t + c = log(lambda_i) (i in S) - the 'Viro' approximate zeros; (b) random points.
Output: out/D_roots.json (decimal strings, 60 digits).  Certification is done in D_certify.py."""
import itertools, json, random, sys
import mpmath as mp
mp.mp.dps = 80
d = json.load(open(sys.argv[1])); A = d["A"]; C = d["C"]; N = 12
Am = mp.matrix(A); 

def FJ(t):
    e = [sum(A[i][k] * t[k] for k in range(4)) for i in range(N)]
    m0 = max(e); m = [mp.e ** (x - m0) for x in e]           # normalised monomials
    F = mp.matrix([sum(C[j][i] * m[i] for i in range(N)) for j in range(4)])
    J = mp.matrix(4, 4)
    for j in range(4):
        for k in range(4):
            J[j, k] = sum(C[j][i] * m[i] * A[i][k] for i in range(N))
    return F, J

def newton(t, it=200):
    t = mp.matrix(t)
    F, J = FJ(t); nf = mp.norm(F)
    for _ in range(it):
        try:
            dt = mp.lu_solve(J, -F)
        except ZeroDivisionError:
            return None
        lam = mp.mpf(1)
        while lam > mp.mpf(10) ** -8:
            t2 = t + lam * dt; F2, J2 = FJ(t2); n2 = mp.norm(F2)
            if n2 < nf: break
            lam /= 2
        else:
            return None
        t, F, J, nf = t2, F2, J2, n2
        if nf < mp.mpf(10) ** -70:
            return t
        if mp.norm(t, mp.inf) > 50: return None
    return None

seeds = []
for S in itertools.combinations(range(N), 5):
    M = mp.matrix([[C[j][i] for i in S[:4]] for j in range(4)])
    lam = mp.lu_solve(M, mp.matrix([-C[j][S[4]] for j in range(4)]))
    lam = list(lam) + [mp.mpf(1)]
    if all(x > 0 for x in lam):
        B = mp.matrix([[A[i][k] for k in range(4)] + [1] for i in S])
        sol = mp.lu_solve(B, mp.matrix([mp.log(x) for x in lam]))
        seeds.append([sol[k] for k in range(4)])
npos = len(seeds)
rng = random.Random(1)
for sc in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
    for _ in range(int(sys.argv[2])):
        seeds.append([mp.mpf(rng.uniform(-sc, sc)) for _ in range(4)])
roots = []; hits = []
for n, s0 in enumerate(seeds):
    r = newton(s0)
    if r is None: continue
    for idx, q in enumerate(roots):
        if mp.norm(r - q, mp.inf) < mp.mpf(10) ** -30:
            hits[idx] += 1; break
    else:
        roots.append(r); hits.append(1)
print("5-subsets with positive dependence (Viro seeds):", npos, "| total seeds:", len(seeds))
tot = 0
for r, h in sorted(zip(roots, hits), key=lambda z: [float(x) for x in z[0]]):
    F, J = FJ(r); sg = int(mp.sign(mp.det(J))); tot += sg
    print("t =", [mp.nstr(x, 12) for x in r], "sign det J =", sg, "hits", h)
print("distinct zeros found:", len(roots), "signed count:", tot)
json.dump([[mp.nstr(x, 60) for x in r] for r in roots], open("out/D_roots.json", "w"))
