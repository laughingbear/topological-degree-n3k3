r"""Claim D, discovery only (float64, NOT a certificate): massive multi-start capped Newton in
log coordinates, vectorised with numpy.  Prints clusters of converged points."""
import json, sys, numpy as np
d = json.load(open(sys.argv[1])); A = np.array(d["A"], float); C = np.array(d["C"], float)
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
def FJ(T):
    E = T @ A.T; m0 = E.max(axis=1, keepdims=True); M = np.exp(E - m0)          # (B,12)
    F = M @ C.T                                                               # (B,4)
    J = np.einsum('bi,ji,ik->bjk', M, C, A)                                   # (B,4,4)
    return F, J
roots = []
nb = int(sys.argv[2])
for scale in (0.05, 0.15, 0.4, 1, 2, 4, 7, 12):
    for cap in (0.3, 0.03):
        T = rng.uniform(-scale, scale, size=(nb, 4))
        for it in range(400):
            F, J = FJ(T)
            dj = np.linalg.det(J); ok = np.abs(dj) > 1e-300
            J[~ok] = np.eye(4)
            dT = -np.linalg.solve(J, F[..., None])[..., 0]
            nrm = np.linalg.norm(dT, axis=1, keepdims=True)
            c = cap * max(scale, 0.05) if it < 300 else 1e9
            dT *= np.minimum(1, c / np.maximum(nrm, 1e-300))
            T = T + dT
            T = np.clip(T, -40, 40)
        F, J = FJ(T)
        good = (np.linalg.norm(F, axis=1) < 1e-9) & (np.abs(T).max(axis=1) < 39)
        for t in T[good]:
            for r in roots:
                if np.abs(r[0] - t).max() < 1e-6 * (1 + np.abs(t).max()):
                    r[1] += 1; break
            else:
                roots.append([t, 1])
tot = 0
for t, h in sorted(roots, key=lambda r: tuple(r[0])):
    F, J = FJ(t[None]); s = int(np.sign(np.linalg.det(J[0]))); tot += s
    print(np.array2string(t, precision=10, floatmode='fixed'), "sign", s, "hits", h, "cond", f"{np.linalg.cond(J[0]):.1e}")
print("clusters:", len(roots), "signed:", tot)
np.save("out/D_candidates.npy", np.array([r[0] for r in roots]))
