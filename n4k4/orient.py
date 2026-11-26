"""Orient every sphere; write pool file for the C searcher: per sphere, list of (subset index among C(9,4), sign)."""
import json, itertools, sys
S = json.load(open('spheres.json'))
sub = list(itertools.combinations(range(9), 4)); sid = {s: i for i, s in enumerate(sub)}
def orient(F):
    F = [tuple(f) for f in F]
    tri = {}
    for f in F:
        for i in range(4):
            tri.setdefault(f[:i] + f[i+1:], []).append((f, i))
    eps = {F[0]: 1}; st = [F[0]]
    while st:
        f = st.pop()
        for i in range(4):
            t = f[:i] + f[i+1:]
            for g, j in tri[t]:
                if g == f: continue
                e = -eps[f] * (-1) ** i * (-1) ** j
                if g in eps: assert eps[g] == e
                else: eps[g] = e; st.append(g)
    assert len(eps) == len(F)
    return [(sid[f], eps[f]) for f in F]
for V in range(5, 10):
    with open(f'pool{V}.txt', 'w') as out:
        L = S[str(V)]; out.write(f'{len(L)}\n')
        for F in L:
            o = orient(F); out.write(str(len(o)) + ' ' + ' '.join(f'{a} {e}' for a, e in o) + '\n')
    print(V, len(L), max(len(F) for F in L))
