r"""Claim D, step 3: COMPLETENESS search for the positive zeros (FLOATING POINT - evidence, not proof).

Work in t = log x.  A zero needs  C m = 0  with m_i = exp(a_i.t) > 0.

Three kinds of cells cover R^4:  the cube [-R0,R0]^4 (R0 = 1/16), subdivided into cubes;  conic
shells {rho*nu : nu in a patch of the boundary of [-1,1]^4, R0 <= rho <= R}, subdivided in angle and
geometrically in rho;  and the far field rho >= R.
Cube B = c + [-r,r]^4 (for a shell cell: exponent ranges from the corner values, same tests):  m_i in [l_i,h_i] = exp(a_i.c -+ r|a_i|_1).
  (LP)  B has no zero if some w in R^4 has  sum_i min((w.C_i) l_i, (w.C_i) h_i) > 0.
        (linear program in (w,z); solved with HiGHS; the returned w is re-evaluated directly and
        must clear a relative margin 1e-7.)
  (K)   for small r: Krawczyk operator in t with interval Jacobian
        J(B)_jk in sum_i C_ji a_ik [l_i,h_i]:   K(B) ∩ B = ∅  -> no zero;
        K(B) ⊂ int B -> exactly one zero (recorded);  otherwise subdivide.
Far field (max|t_k| >= R): cells {rho*nu : nu in a patch P of the boundary of [-1,1]^4, rho >= R}.
  With g_i^-+ = min/max over P of a_i.nu, G = max_i g_i^-, i* a maximiser:  relative to exp(rho G),
  (all monomials are divided by the monomial k that is largest at the cell centre:  m_k = 1,
  m_i <= exp(R g_i^+) if g_i^+ < 0,  m_i in [0,inf) otherwise; g^-+ computed for a_i - a_k.)
  Same LP with these bounds (h = inf  ->  constraint w.C_i >= 0).
Every cell is excluded, certified to contain one zero, or subdivided; the run ends when the
stack is empty.  Floating point: results are numerical evidence only.
usage: D_exclude.py system.json R part nparts [comma-separated job indices]   (73 jobs for R = 8)
"""
import json, sys, time, itertools
import numpy as np
from scipy.optimize import linprog

d = json.load(open(sys.argv[1])); A = np.array(d["A"], float); C = np.array(d["C"], float)
R = float(sys.argv[2]); part = int(sys.argv[3]); nparts = int(sys.argv[4])
N = 12; S1 = A.sum(axis=1)
MARGIN = 1e-7
PROG = [0.0]; t0 = time.time()
SELF = [0, 0]                      # self-test: (points tested, failures)
_rng = np.random.default_rng(12345)
def selftest(w, k, pts):
    """direct evaluation: phi_w(t) = sum_i (w.C_i) exp((a_i-a_k).t) must be > 0 at every sample point
    of a cell that was excluded with the vector w."""
    p = C.T @ w; D = A - A[k]
    for t in pts:
        e = D @ t
        with np.errstate(over='ignore', invalid='ignore'):
            v = np.sum(np.where(p == 0, 0.0, p * np.exp(np.minimum(e, 700))))
        if np.any((e > 700) & (p < 0)): v = -1.0
        SELF[0] += 1; SELF[1] += not (v > 0)

stats = dict(lp=0, lp_excl=0, quick=0, kr_excl=0, kr_one=0, split=0, far_excl=0, far_split=0)

def gval(w, l, h, inf):
    p = C.T @ w
    v = np.where(inf, np.where(p >= 0, p * l, -np.inf), np.minimum(p * l, p * h))
    return v.sum()

def lp_w(l, h, inf):
    """maximise sum z, z_i <= (w.C_i) l_i, z_i <= (w.C_i) h_i (finite h) or w.C_i >= 0 (h = inf)."""
    stats['lp'] += 1
    if stats['lp'] % 20000 == 0:
        print(f"# progress {time.time()-t0:.0f}s stats {stats} roots {len(roots)} done-fraction-of-current-job {PROG[0]:.6f}", flush=True)
    nv = 4 + N; c = np.zeros(nv); c[4:] = -1
    rows = []; rhs = []
    for i in range(N):
        row = np.zeros(nv); row[:4] = -C[:, i] * l[i]; row[4 + i] = 1; rows.append(row); rhs.append(0)
        row = np.zeros(nv)
        if inf[i]:
            row[:4] = -C[:, i]; rows.append(row); rhs.append(-1e-6 * np.abs(C[:, i]).sum())   # w.C_i >= small margin
        else:
            row[:4] = -C[:, i] * h[i]; row[4 + i] = 1; rows.append(row); rhs.append(0)
    res = linprog(c, A_ub=np.array(rows), b_ub=np.array(rhs), bounds=[(-1, 1)] * 4 + [(None, None)] * N, method="highs")
    if res.status != 0:
        return None
    return res.x[:4]

def excluded(l, h, inf, whint):
    # valid relaxations that keep the LP well scaled (ranges are relative to a monomial equal to 1):
    # a huge upper bound is replaced by +infinity, a tiny lower bound by 0.
    inf = inf | (h > 1e12); l = np.where(l < 1e-12, 0.0, l); h = np.where(inf, 1.0, h)
    scale = 1.0
    for w in ([whint] if whint is not None else []):
        if gval(w, l, h, inf) > MARGIN * scale * np.abs(C.T @ w).sum():
            stats['quick'] += 1; return True, w
    w = lp_w(l / scale, np.where(inf, 1, h / scale), inf)
    if w is not None and gval(w, l, h, inf) > MARGIN * scale * np.abs(C.T @ w).sum():
        stats['lp_excl'] += 1; return True, w
    return False, (w if w is not None else whint)

def krawczyk(c, r):
    e = A @ c; M = e.max(); m = np.exp(e - M); lo = m * np.exp(-r * S1); hi = m * np.exp(r * S1)
    F = C @ m
    J = (C * m) @ A
    try:
        Y = np.linalg.inv(J)
    except np.linalg.LinAlgError:
        return 'split', None
    CA = C[:, :, None] * A[None, :, :]                          # (4,12,4): C_ji a_ik
    Jlo = np.where(CA >= 0, CA * lo[None, :, None], CA * hi[None, :, None]).sum(axis=1)
    Jhi = np.where(CA >= 0, CA * hi[None, :, None], CA * lo[None, :, None]).sum(axis=1)
    # E = I - Y J(B): interval product
    Ylo = np.minimum(Y[:, :, None] * Jlo[None], Y[:, :, None] * Jhi[None]).sum(axis=1)
    Yhi = np.maximum(Y[:, :, None] * Jlo[None], Y[:, :, None] * Jhi[None]).sum(axis=1)
    Elo = np.eye(4) - Yhi; Ehi = np.eye(4) - Ylo
    rad = np.maximum(np.abs(Elo), np.abs(Ehi)).sum(axis=1) * r * (1 + 1e-9) + 1e-13
    mid = -Y @ F                                               # K = c + mid + [-rad, rad]
    if np.any(np.abs(mid) - rad > r * (1 + 1e-9)):
        return 'empty', None
    if np.all(np.abs(mid) + rad < r * (1 - 1e-9)):
        return 'one', c + mid
    return 'split', None

roots = []
def inner(c, r, whint):
    r_top = r; PROG[0] = 0.0
    stack = [(c, r, whint)]
    while stack:
        c, r, whint = stack.pop()
        e = A @ c; k = int(e.argmax()); D = A - A[k]; e = D @ c; sp = r * np.abs(D).sum(axis=1)
        l = np.exp(e - sp); h = np.exp(e + sp)             # monomials divided by monomial k
        ok, w = excluded(l, h, np.zeros(N, bool), whint)
        if ok:
            selftest(w, k, c + r * _rng.uniform(-1, 1, size=(3, 4)))
            PROG[0] += (r / r_top) ** 4; continue
        if r <= 0.02:
            st, z = krawczyk(c, r)
            if st == 'empty':
                stats['kr_excl'] += 1; PROG[0] += (r / r_top) ** 4; continue
            if st == 'one':
                stats['kr_one'] += 1; roots.append(z); PROG[0] += (r / r_top) ** 4; continue
        stats['split'] += 1
        if r < 1e-9:
            print("UNRESOLVED tiny box at", c, flush=True); continue
        for s in itertools.product((-0.5, 0.5), repeat=4):
            stack.append((c + r * np.array(s), r / 2, w))

def far(face, sign, c3, r3, whint):
    """patch of the face {nu_face = sign}, other three coords in c3 +- r3."""
    stack = [(c3, r3, whint)]
    others = [k for k in range(4) if k != face]
    while stack:
        c3, r3, whint = stack.pop()
        nu = np.zeros(4); nu[face] = sign; nu[others] = c3
        k = int((A @ nu).argmax()); D = A - A[k]
        g0 = D @ nu; sp = r3 * np.abs(D[:, others]).sum(axis=1)
        gm = g0 - sp; gp = g0 + sp                      # exponents of m_i/m_k are rho*[gm,gp], rho >= R
        inf = gp >= 0; inf[k] = False                    # m_i/m_k in [0, inf)  ->  need w.C_i >= 0
        l = np.zeros(N); l[k] = 1.0
        h = np.where(inf, 1.0, np.exp(np.maximum(R * np.minimum(gp, 0.0), -700))); h[k] = 1.0
        ok, w = excluded(l, h, inf, whint)
        if ok:
            pts = []
            for _ in range(3):
                t = np.zeros(4); t[face] = sign; t[others] = c3 + r3 * _rng.uniform(-1, 1, size=3)
                pts.append(t * R * np.exp(_rng.uniform(0, 3)))
            selftest(w, k, pts)
            stats['far_excl'] += 1; continue
        stats['far_split'] += 1
        if r3 < 1e-6:
            print("UNRESOLVED far patch", face, sign, c3, flush=True); continue
        for s in itertools.product((-0.5, 0.5), repeat=3):
            stack.append((c3 + r3 * np.array(s), r3 / 2, w))


def shell(face, sign, c3, r3, rho1, rho2, whint):
    """cells {rho*nu : nu in patch(c3 +- r3) of face, rho in [rho1, rho2]}"""
    others = [k for k in range(4) if k != face]
    top = np.log(rho2 / rho1) * 8.0
    stack = [(c3, r3, rho1, rho2, whint)]
    PROG[0] = 0.0
    while stack:
        c3, r3, p1, p2, whint = stack.pop()
        nu = np.zeros(4); nu[face] = sign; nu[others] = c3
        k = int((A @ nu).argmax()); D = A - A[k]
        g0 = D @ nu; sp = r3 * np.abs(D[:, others]).sum(axis=1)
        gm = g0 - sp; gp = g0 + sp
        elo = np.minimum(p1 * gm, p2 * gm); ehi = np.maximum(p1 * gp, p2 * gp)
        ok, w = excluded(np.exp(np.maximum(elo, -700)), np.exp(np.minimum(ehi, 700)), np.zeros(N, bool), whint)
        vol = np.log(p2 / p1) * (2 * r3) ** 3 / top
        if ok:
            pts = []
            for _ in range(3):
                t = np.zeros(4); t[face] = sign; t[others] = c3 + r3 * _rng.uniform(-1, 1, size=3)
                pts.append(t * np.exp(_rng.uniform(np.log(p1), np.log(p2))))
            selftest(w, k, pts)
            PROG[0] += vol; continue
        # bounding cube in t-space
        lo = np.zeros(4); hi = np.zeros(4)
        lo[face] = min(sign * p1, sign * p2); hi[face] = max(sign * p1, sign * p2)
        for k, o in enumerate(others):
            cands = [p * (c3[k] + e * r3) for p in (p1, p2) for e in (-1, 1)]
            lo[o] = min(cands); hi[o] = max(cands)
        rc = (hi - lo).max() / 2
        if rc <= 0.02:
            st, z = krawczyk((lo + hi) / 2, rc)
            if st == 'empty':
                stats['kr_excl'] += 1; PROG[0] += vol; continue
            if st == 'one':
                stats['kr_one'] += 1; roots.append(z); PROG[0] += vol; continue
        stats['split'] += 1
        if rc < 1e-9:
            print("UNRESOLVED tiny cell near", (lo + hi) / 2, flush=True); continue
        if (p2 - p1) > 2 * r3 * p2:
            pm = np.sqrt(p1 * p2)
            stack.append((c3, r3, p1, pm, w)); stack.append((c3, r3, pm, p2, w))
        else:
            for sg in itertools.product((-0.5, 0.5), repeat=3):
                stack.append((c3 + r3 * np.array(sg), r3 / 2, p1, p2, w))

R0 = 0.0625
jobs = [('in', np.zeros(4))]
for face in range(4):
    for sign in (-1, 1):
        for sg in itertools.product((-0.5, 0.5), repeat=3):
            jobs.append(('shell', face, sign, np.array(sg)))
        jobs.append(('far', face, sign))
only = set(map(int, sys.argv[5].split(','))) if len(sys.argv) > 5 else None
for idx, job in enumerate(jobs):
    if (idx % nparts != part) if only is None else (idx not in only):
        continue
    if job[0] == 'in':
        inner(job[1], R0, None)
    elif job[0] == 'shell':
        shell(job[1], job[2], job[3], 0.5, R0, R, None)
    else:
        far(job[1], job[2], np.zeros(3), 1.0, None)
    print(f"# job {idx} {job[0]} done at {time.time()-t0:.0f}s  stats {stats}  roots so far {len(roots)}  selftest points/failures {SELF}", flush=True)
uniq = []
for z in roots:
    if not any(np.abs(z - u).max() < 1e-7 for u in uniq):
        uniq.append(z)
tot = 0
for z in sorted(uniq, key=tuple):
    e = A @ z; m = np.exp(e - e.max()); J = (C * m) @ A; s = int(np.sign(np.linalg.det(J))); tot += s
    print("ZERO t =", np.array2string(z, precision=12, floatmode='fixed'), "sign", s, flush=True)
print(f"DONE part {part}/{nparts}: zeros {len(uniq)} signed {tot} stats {stats} selftest points/failures {SELF} time {time.time()-t0:.0f}s")
json.dump([list(map(float, z)) for z in uniq], open(f"out/D_excl_part{part}{'_' + sys.argv[5].replace(',', '-') if only else ''}.json", "w"))
