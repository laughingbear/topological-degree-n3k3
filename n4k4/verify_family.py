"""Exact check of data/D44_family.json: one exponent set A (n = 4, N = 9) and seven coefficient matrices with |deg F| = 0..6,
and of the nine degree-6 systems in data/deg6_systems.json."""
import json, random
from verify_exact import structure, degree, UC

def degs(fac, C, k=40):
    random.seed(7); out = []
    while len(out) < k:
        try: out.append(degree(fac, C, [random.randint(-99, 99) for _ in range(4)])[0])
        except ValueError: pass
    return sorted(set(out))

S = json.load(open('data/D44_family.json')); A = S['A']; UA, fac = structure(A)
print('family: (U_A) =', UA, 'facets =', len(fac)); ok = UA
for k, C in S['C_by_degree'].items():
    d = degs(fac, C); good = UC(C) and len(d) == 1 and abs(d[0]) == int(k); ok &= good
    print('  target', k, '(U_C) =', UC(C), 'deg F over 40 generic u:', d, 'OK' if good else 'FAIL')
for k, s in json.load(open('data/deg6_systems.json')).items():
    UA, fac = structure(s['A']); d = degs(fac, s['C']); good = UA and UC(s['C']) and len(d) == 1 and abs(d[0]) == 6; ok &= good
    print('  sphere type', k, 'facets', len(fac), 'deg F:', d, 'OK' if good else 'FAIL')
print('ALL OK' if ok else 'FAILURE')
