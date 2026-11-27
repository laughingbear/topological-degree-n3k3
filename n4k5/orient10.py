import json, itertools, random, sys
L = json.load(open('spheres10.json'))
sub = list(itertools.combinations(range(10), 4)); sid = {s: i for i, s in enumerate(sub)}
def orient(F):
    F = [tuple(f) for f in F]; tri = {}
    for f in F:
        for i in range(4): tri.setdefault(f[:i] + f[i+1:], []).append((f, i))
    eps = {F[0]: 1}; st = [F[0]]
    while st:
        f = st.pop()
        for i in range(4):
            for g, j in tri[f[:i] + f[i+1:]]:
                if g == f: continue
                e = -eps[f] * (-1) ** i * (-1) ** j
                if g in eps: assert eps[g] == e
                else: eps[g] = e; st.append(g)
    assert len(eps) == len(F)
    return [(sid[f], eps[f]) for f in F]
with open('pool10_all.txt', 'w') as out:
    out.write(f'{len(L)}\n')
    for F in L:
        o = orient(F); out.write(str(len(o)) + ' ' + ' '.join(f'{a} {e}' for a, e in o) + '\n')
from collections import Counter
print(sorted(Counter(len(F) for F in L).items()))
