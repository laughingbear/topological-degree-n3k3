"""Biquadratic final polynomials (Bokowski / Richter-Gebert) for the 4 uniform rank-5 chirotopes compatible with sphere 101272.
For sigma (3-set) and a<b<c<d:  [sab][scd] - [sac][sbd] + [sad][sbc] = 0.  With the signs known, the three signed terms
t1, t2, t3 sum to 0, so the term whose sign is in the minority is larger in absolute value than each of the other two:
log|.| inequalities, linear in y_S = log|chi_S|.  If the resulting strict linear system is infeasible, the chirotope is not realizable.
An infeasibility certificate is a nonnegative integer combination of the rows that sums to the zero vector; it is found by LP
and then verified in exact integer arithmetic."""
import itertools, json, numpy as np, sys
from fractions import Fraction as Fr
from scipy.optimize import linprog
from polysat import gp5, parity
V = 10
X, cl, chi = gp5(V); NS = len(X)
D = json.load(open('chi27.json')); free = D['free']; units = D['units']
def rows_for(sign):                       # sign: dict var -> +-1
    def sg(t):
        l = chi(t); return sign[abs(l)] * (1 if l > 0 else -1)
    rows = []
    for sig in itertools.combinations(range(V), 3):
        rest = [z for z in range(V) if z not in sig]
        for a, b, c, d in itertools.combinations(rest, 4):
            T = [((sig + (a, b)), (sig + (c, d)), 1), ((sig + (a, c)), (sig + (b, d)), -1), ((sig + (a, d)), (sig + (b, c)), 1)]
            s = [sg(p) * sg(q) * k for p, q, k in T]
            assert abs(sum(s)) == 1, 'chirotope violates GP'
            minority = s.index(-1) if sum(s) == 1 else s.index(1)
            for j in range(3):
                if j == minority: continue
                r = {}
                for S_, coef in ((T[minority][0], 1), (T[minority][1], 1), (T[j][0], -1), (T[j][1], -1)):
                    v = X[tuple(sorted(S_))] - 1; r[v] = r.get(v, 0) + coef
                rows.append(r)            # sum r[v]*y_v > 0
    return rows
results = []
for mi, model in enumerate(D['models']):
    sign = {abs(l): (1 if l > 0 else -1) for l in units}; sign.update({abs(l): (1 if l > 0 else -1) for l in model})
    rows = rows_for(sign); M = np.zeros((len(rows), NS))
    for i, r in enumerate(rows):
        for v, c in r.items(): M[i, v] = c
    # primal feasibility: M y >= 1
    res = linprog(np.zeros(NS), A_ub=-M, b_ub=-np.ones(len(rows)), bounds=[(None, None)] * NS, method='highs')
    print('chirotope', mi, 'rows', len(rows), 'LP status:', res.status, res.message[:60], flush=True)
    cert = None
    if res.status == 2:                    # infeasible -> find lambda >= 0, lambda^T M = 0, sum lambda = 1  (Farkas)
        r2 = linprog(np.zeros(len(rows)), A_eq=np.vstack([M.T, np.ones((1, len(rows)))]), b_eq=np.concatenate([np.zeros(NS), [1]]), bounds=[(0, None)] * len(rows), method='highs-ds')
        if r2.status == 0:
            lam = r2.x; supp = [i for i in range(len(rows)) if lam[i] > 1e-9]
            # exact: solve for rational weights on the support with sum = 1 via integer nullspace search (small support) 
            import sympy
            Ms = sympy.Matrix([[int(M[i, v]) for i in supp] for v in range(NS)])
            ns = Ms.nullspace()
            ok = False
            for vec in ns:
                den = sympy.ilcm(*[sympy.fraction(x)[1] for x in vec]); w = [int(x * den) for x in vec]
                if all(x <= 0 for x in w): w = [-x for x in w]
                if all(x >= 0 for x in w) and any(w):
                    tot = [sum(w[k] * int(M[supp[k], v]) for k in range(len(supp))) for v in range(NS)]
                    ok = all(t == 0 for t in tot); cert = {'rows': [rows[i] for i in supp], 'weights': w}; break
            print('   Farkas support', len(supp), 'nullspace dim', len(ns), 'exact integer certificate verified:', ok, flush=True)
            if not ok: cert = None
    results.append({'model': mi, 'lp_status': int(res.status), 'certificate': cert})
json.dump(results, open('bfp27.json', 'w'))
print('non-realizable with exact certificate:', sum(1 for r in results if r['certificate']), 'of', len(results))
