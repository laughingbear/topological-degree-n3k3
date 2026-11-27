"""Item 5(c): active search for counterexamples.
 (i)  hill-climb on g = 4 deg - (2F - E/2 - q/2) over vertex maps (moves: Gaussian perturbation of one C_v,
      sign flip C_v -> -C_v), on random simplicial 4-polytopes and on cyclic polytopes C(V,4);
      also tracks the same-sign refinement  4 deg < sum_e floor((m_e-1)/2) - F  whenever q = 0;
 (ii) maximise D_t over 4x4 bases (Nelder-Mead + random restarts) to probe D_t < 1/2.
Fast float degree for the search; every record-setting map is re-verified with exact integer arithmetic.
"""
import sys, itertools, math, time, json
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from core import coherent_orientation, edges_of, VertexMap
from randmaps import random_polytope
from cyclic import gale_facets


class FastK:
    def __init__(self, facets):
        self.facets, self.o = coherent_orientation(facets)
        self.o = np.array(self.o)
        self.idx = np.array(self.facets)
        self.V = int(self.idx.max()) + 1
        self.E_list = edges_of(self.facets)
        self.E, self.F = len(self.E_list), len(self.facets)
        me = {e: 0 for e in self.E_list}
        for f in self.facets:
            for e in itertools.combinations(f, 2):
                me[e] += 1
        self.refine = sum((m - 1) // 2 for m in me.values()) - self.F

    def evaluate(self, C, us):
        M = C[self.idx].transpose(0, 2, 1)  # F x 4 x 4, columns = vectors
        d = np.linalg.det(M)
        if np.min(np.abs(d)) < 1e-9:
            return None
        s = self.o * np.sign(d)
        c = np.linalg.solve(M[None, :, :, :], us[:, None, :, None])[..., 0]  # U x F x 4
        if np.min(np.abs(c)) < 1e-10:
            return None
        cnt = ((c > 0).all(-1) * s[None, :]).sum(1)
        if len(set(cnt.tolist())) != 1:
            return None
        deg = int(cnt[0])
        if deg < 0:
            s = -s
            deg = -deg
        q = int((s < 0).sum())
        return deg, q


def thmC(V):
    i = np.arange(V); t = np.pi * i / V
    return ((-1.0) ** i)[:, None] * np.column_stack([np.cos(t), np.sin(t), np.cos(3 * t), np.sin(3 * t)])


def climb(K, rng, iters=3000, C0=None, temp=0.0):
    us = rng.standard_normal((3, 4))
    C = rng.standard_normal((K.V, 4)) if C0 is None else C0.copy()
    r = K.evaluate(C, us)
    while r is None:
        C = rng.standard_normal((K.V, 4))
        r = K.evaluate(C, us)
    bound = lambda q: 2 * K.F - K.E / 2 - q / 2
    g = 4 * r[0] - bound(r[1])
    best = (g, r, C.copy())
    for it in range(iters):
        C2 = C.copy()
        v = rng.integers(K.V)
        if rng.random() < 0.15:
            C2[v] = -C2[v]
        else:
            C2[v] += rng.standard_normal(4) * rng.choice([0.05, 0.3, 1.0])
        r2 = K.evaluate(C2, us)
        if r2 is None:
            continue
        g2 = 4 * r2[0] - bound(r2[1])
        if g2 >= g or (temp > 0 and rng.random() < math.exp((g2 - g) / temp)):
            C, g, r = C2, g2, r2
            if g > best[0]:
                best = (g, r, C.copy())
    return best


def verify_exact(K, C, rng):
    Ci = [[int(round(x * 10**6)) for x in row] for row in C]
    vm = VertexMap(K.facets, Ci)
    deg = vm.degree(rng)
    if deg < 0:
        vm.flip(); deg = -deg
    q = vm.s.count(-1)
    return deg, q, 4 * deg < 2 * vm.F - vm.E / 2 - q / 2


def D_of(Mflat):
    from scipy.stats import multivariate_normal
    M = Mflat.reshape(4, 4)
    if abs(np.linalg.det(M)) < 1e-8:
        return 0.0
    Mi = np.linalg.inv(M)
    S = Mi @ Mi.T
    tot = 0.0
    for k in range(4):
        pat = np.ones(4); pat[k] = -1
        S2 = S * pat[:, None] * pat[None, :]
        dd = 1 / np.sqrt(np.diag(S2))
        R = S2 * dd[:, None] * dd[None, :]
        R = (R + R.T) / 2
        try:
            tot += multivariate_normal.cdf(np.zeros(4), mean=np.zeros(4), cov=R, allow_singular=True,
                                           abseps=1e-7, releps=1e-7)
        except np.linalg.LinAlgError:
            return 0.0
    return tot


if __name__ == "__main__":
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 600
    rng = np.random.default_rng(777)
    t0 = time.time()
    Ks = []
    for V in range(7, 14, 2):
        Ks.append((f"cyclic C({V},4)", FastK(gale_facets(V))))
    for V in (8, 10, 12):
        Ks.append((f"cyclic C({V},4)", FastK(gale_facets(V))))
    for _ in range(14):
        V = int(rng.integers(6, 13))
        Ks.append((f"random V={V}", FastK(random_polytope(V, rng)[1])))
    report = []
    full_results = []
    if len(sys.argv) > 2 and sys.argv[2] == "Donly":
        Ks = []
    per = budget * 0.85 / max(len(Ks), 1)
    for name, K in Ks:
        tK = time.time()
        bestg, bestr, bestC = -1e9, None, None
        refine_viol = 0
        runs = 0
        seedC = thmC(K.V) if (name.startswith("cyclic") and K.V % 2) else None
        while time.time() - tK < per:
            if seedC is not None and runs % 2 == 0:
                C0 = seedC + 0.05 * rng.standard_normal(seedC.shape)
            else:
                C0 = None
            g, r, C = climb(K, rng, iters=1500, C0=C0, temp=[0.0, 0.5, 1.5][runs % 3])
            runs += 1
            if r[1] == 0 and 4 * r[0] >= K.refine:
                refine_viol += 1
            if g > bestg:
                bestg, bestr, bestC = g, r, C
        ex = verify_exact(K, bestC, rng)
        # full check of identity / Gauss-Bonnet / edge bound / D_t on the record map
        from randmaps import analyse
        vmx = VertexMap(K.facets, [[int(round(x * 10**6)) for x in row] for row in bestC])
        rr = analyse(vmx, rng, name)
        full_results.append(rr)
        line = (f"{name}: V={K.V} E={K.E} F={K.F} restarts={runs} best deg={bestr[0]} q={bestr[1]} "
                f"4deg-bound={bestg:+.1f} (exact recheck deg={ex[0]} q={ex[1]} bound holds={ex[2]}) "
                f"same-sign refinement violations={refine_viol} | full analyse ok={rr['ok']} idErr={rr['identity_err']:.1e} "
                f"sum_wv={rr['sum_wv']} sum_we={rr['sum_we']} Dmax={rr['Dmax']:.3f}")
        print(line, flush=True)
        report.append(line)
    if Ks: print("record maps fully analysed:", len(full_results), " failures:", sum(1 for r in full_results if not r["ok"]),
          " degrees:", [r["deg"] for r in full_results], flush=True)
    # D_t: (a) random bases, (b) Nelder-Mead from the best random ones, (c) the explicit family
    # c = M^{-1}u with Cov(c) = a a^T + eps I, a = (1,1,1,-3)  (D -> 1/2 as eps -> 0, never reaching it)
    from scipy.optimize import minimize
    rs = [(D_of(x), x) for x in rng.standard_normal((300, 16))]
    rs.sort(key=lambda z: -z[0])
    print(f"max D_t over 300 random Gaussian bases: {rs[0][0]:.6f}", flush=True)
    bestD = rs[0][0]
    for D0, x0 in rs[:3]:
        res = minimize(lambda x: -D_of(x), x0, method="Nelder-Mead", options=dict(maxiter=300))
        bestD = max(bestD, -res.fun)
    print(f"max D_t after Nelder-Mead polishing: {bestD:.6f}  (claim: < 1/2)", flush=True)
    a = np.array([1, 1, 1, -3.0])
    for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
        S = np.outer(a, a) + eps * np.eye(4)
        L = np.linalg.cholesky(S)          # c = L z  =>  M = L^{-1}
        print(f"  family eps={eps:g}: D = {D_of(np.linalg.inv(L).ravel()):.6f}", flush=True)
