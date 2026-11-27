"""Exact checks that need no SAT solver:
 (1) system_deg8.json and D45_family.json: (U_A), (U_C), |deg F| by Bihan's Theorem 1.3;
 (2) vertexmap_deg9_nonpolytopal.json: every 4 columns independent, degree 9, all facets of one sign;
 (3) bfp27.json: the four biquadratic final polynomials, re-derived from the chirotopes in chi27.json;
 (4) realizations_22.json: 10 vertices in general position (exact)."""
import json, itertools, random
from fractions import Fraction as Fr
from verify_exact import structure, degree, UC, det, cols, sgn
from polysat import gp5
random.seed(7); ok = True
def degs(fac, C, k=30):
    out = []
    while len(out) < k:
        try: out.append(degree(fac, C, [random.randint(-99, 99) for _ in range(4)])[0])
        except ValueError: pass
    return sorted(set(out))
S = json.load(open('system_deg8.json')); UA, fac = structure(S['A']); d = degs(fac, S['C'])
print('(1) system_deg8: (U_A)', UA, '(U_C)', UC(S['C']), 'facets', len(fac), 'deg F', d); ok &= UA and UC(S['C']) and d == [8]
F = json.load(open('D45_family.json')); UA, fac = structure(F['A'])
for k, C in F['C_by_degree'].items():
    d = degs(fac, C); good = UA and UC(C) and len(d) == 1 and abs(d[0]) == int(k); ok &= good
    print('    family target', k, 'deg F', d, 'OK' if good else 'FAIL')
M = json.load(open('vertexmap_deg9_nonpolytopal.json')); C = M['C']; of = [(tuple(G), e) for G, e in M['oriented_facets']]
ind = all(det(cols(C, T)) != 0 for T in itertools.combinations(range(10), 4)); sg = {e * sgn(det(cols(C, G))) for G, e in of}
ds = set()
for _ in range(40):
    u = [random.randint(-60, 60) for _ in range(4)]; t = 0; gen = True
    for G, e in of:
        CG = cols(C, G); dd = det(CG); lam = [det([[CG[r][j] if j != c else u[r] for j in range(4)] for r in range(4)]) / dd for c in range(4)]
        if any(l == 0 for l in lam): gen = False
        if all(l > 0 for l in lam): t += e * sgn(dd)
    if gen: ds.add(t)
print('(2) degree-9 vertex map: columns independent', ind, 'degrees', ds, 'facet signs', sg); ok &= ind and ds == {9} and len(sg) == 1
X, cl, chi = gp5(10); D = json.load(open('chi27.json')); B = json.load(open('bfp27.json'))
for mi, model in enumerate(D['models']):
    sign = {abs(l): (1 if l > 0 else -1) for l in D['units']}; sign.update({abs(l): (1 if l > 0 else -1) for l in model})
    sg_ = lambda t: sign[abs(chi(t))] * (1 if chi(t) > 0 else -1)
    valid = set()
    for sig in itertools.combinations(range(10), 3):
        rest = [z for z in range(10) if z not in sig]
        for a, b, c, d in itertools.combinations(rest, 4):
            T = [((sig + (a, b)), (sig + (c, d)), 1), ((sig + (a, c)), (sig + (b, d)), -1), ((sig + (a, d)), (sig + (b, c)), 1)]
            s = [sg_(p) * sg_(q) * k for p, q, k in T]; assert abs(sum(s)) == 1
            mn = s.index(-1) if sum(s) == 1 else s.index(1)
            for j in range(3):
                if j == mn: continue
                r = {}
                for S_, co in ((T[mn][0], 1), (T[mn][1], 1), (T[j][0], -1), (T[j][1], -1)):
                    v = X[tuple(sorted(S_))] - 1; r[v] = r.get(v, 0) + co
                valid.add(tuple(sorted((v, c_) for v, c_ in r.items() if c_)))
    cert = B[mi]['certificate']; rows = [tuple(sorted((int(v), c_) for v, c_ in r.items() if c_)) for r in cert['rows']]; w = cert['weights']
    tot = {}
    for r, wt in zip(rows, w):
        for v, c_ in r: tot[v] = tot.get(v, 0) + wt * c_
    good = all(r in valid for r in rows) and all(x >= 0 for x in w) and any(w) and all(t == 0 for t in tot.values()); ok &= good
    print('(3) chirotope', mi, 'inequalities', len(valid), 'certificate rows', len(rows), 'OK' if good else 'FAIL')
R = json.load(open('realizations_22.json')); g = 0
for k, A in R.items():
    UA, fac = structure(A); g += UA and len({v for G, e in fac for v in G}) == 10
print('(4) realizations with 10 vertices in general position:', g, 'of', len(R)); ok &= g == len(R)
print('ALL OK' if ok else 'FAILURE')
