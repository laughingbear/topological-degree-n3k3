r"""Merge candidate zeros (60-digit strings from D_seek.py, float64 points from D_seek_fast.py /
D_exclude.py), polish each with multiprecision Newton in t = log x, dedupe, write
out/D_roots_all.json for D_certify.py.  (Numerical; certification happens in D_certify.py.)"""
import json, sys, glob
import numpy as np, mpmath as mp
mp.mp.dps = 80
d = json.load(open(sys.argv[1])); A = d["A"]; C = d["C"]; N = 12
def FJ(t):
    e = [sum(A[i][k] * t[k] for k in range(4)) for i in range(N)]
    m0 = max(e); m = [mp.e ** (x - m0) for x in e]
    F = mp.matrix([sum(C[j][i] * m[i] for i in range(N)) for j in range(4)])
    J = mp.matrix(4, 4)
    for j in range(4):
        for k in range(4):
            J[j, k] = sum(C[j][i] * m[i] * A[i][k] for i in range(N))
    return F, J
cands = [[mp.mpf(s) for s in r] for r in json.load(open("out/D_roots.json"))]
cands += [[mp.mpf(float(x)) for x in r] for r in np.load("out/D_candidates.npy")]
for f in glob.glob("out/D_excl_part*.json"):
    cands += [[mp.mpf(x) for x in r] for r in json.load(open(f))]
roots = []
for c in cands:
    t = mp.matrix(c)
    for _ in range(12):
        F, J = FJ(t); t = t + mp.lu_solve(J, -F)
    F, J = FJ(t)
    assert mp.norm(F) < mp.mpf(10) ** -60, mp.norm(F)
    if not any(mp.norm(t - r, mp.inf) < mp.mpf(10) ** -30 for r in roots):
        roots.append(t)
roots.sort(key=lambda r: [float(x) for x in r])
json.dump([[mp.nstr(x, 60) for x in r] for r in roots], open("out/D_roots_all.json", "w"), indent=0)
print("distinct polished zeros:", len(roots))
