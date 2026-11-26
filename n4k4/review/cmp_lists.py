import json, sys, collections
from sph import *
L = parse_lutz('/home/claude/ext/3_manifolds_all.txt', '/home/claude/ext/3_manifolds_all_type.txt')
print('Lutz entries', len(L), collections.Counter((V, t) for V, F, t in L.values()))
auth = json.load(open('/home/claude/author/spheres.json'))
out = {}
for V in range(5, 10):
    lz = {n: F for n, (v, F, t) in L.items() if v == V and t == 'S^3'}
    cl = {}
    for n, F in lz.items():
        ok, fv = manifold_checks(F); assert ok, n
        orient(F); cl.setdefault(canon(F), []).append(n)
    ca = {}
    for i, F in enumerate(auth[str(V)]):
        ok, fv = manifold_checks(F); assert ok, (V, i)
        ca.setdefault(canon(F), []).append(i)
    neigh = sum(1 for F in lz.values() if manifold_checks(F)[1][1] == V * (V - 1) // 2)
    print(f'V={V}: Lutz S^3 {len(lz)} ({len(cl)} distinct canonical forms); author {len(auth[str(V)])} ({len(ca)} distinct); '
          f'common {len(set(cl) & set(ca))}; only Lutz {len(set(cl) - set(ca))}; only author {len(set(ca) - set(cl))}; neighborly in Lutz list {neigh}', flush=True)
    out[V] = {n: ca.get(canon(F), [None])[0] for n, F in lz.items()}
    # the non-sphere must not be in author's list
    for n, (v, F, t) in L.items():
        if v == V and t != 'S^3': print('   non-sphere', n, t, 'in author list:', canon(F) in ca)
json.dump({str(V): out[V] for V in out}, open('lutz_to_author.json', 'w'))
