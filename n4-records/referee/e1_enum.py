r"""E1 cross-check, independent of ours/enum.c and ours/genprob.py: count unions of facet orbits
on V = 2n vertices (two n-gons; rotation only, or dihedral with parameter c) in which every
triangle lies in 0 or 2 facets.  Plain include/exclude backtracking over facet orbits in a fixed
order (a different search strategy from enum.c).  Reports: all non-empty solutions, and those
that are strongly connected (dual graph connected)."""
import itertools, sys
n = int(sys.argv[1]); dih = int(sys.argv[2]); c = int(sys.argv[3]) if len(sys.argv) > 3 else 0
V = 2 * n
def rot(v, s): return (v // n) * n + (v % n + s) % n
def ref(v, s): return ((-v + s) % n) if v < n else n + ((c - (v - n) + s) % n)
G = [lambda v, s=s: rot(v, s) for s in range(n)] + ([lambda v, s=s: ref(v, s) for s in range(n)] if dih else [])
def orbit(f): return frozenset(tuple(sorted(g(v) for v in f)) for g in G)
forb = sorted({orbit(f) for f in itertools.combinations(range(V), 4)}, key=lambda o: min(o))
tri_cnt = []                       # per orbit: dict triangle -> number of facets of the orbit containing it
for o in forb:
    dct = {}
    for f in o:
        for t in itertools.combinations(f, 3): dct[t] = dct.get(t, 0) + 1
    tri_cnt.append(dct)
ok_orb = [i for i, dct in enumerate(tri_cnt) if max(dct.values()) <= 2]
# representative triangles: one per triangle orbit suffices
def torb(t): return min(tuple(sorted(g(v) for v in t)) for g in G)
reps = sorted({torb(t) for t in itertools.combinations(range(V), 3)})
inc = [{t: tri_cnt[i].get(t, 0) for t in reps if tri_cnt[i].get(t, 0)} for i in ok_orb]
last = {}                          # last orbit position touching each rep triangle
for p, dct in enumerate(inc):
    for t in dct: last[t] = p
closing = [[t for t in dct if last[t] == p] for p, dct in enumerate(inc)]
cnt = {t: 0 for t in reps}; chosen = []; total = 0; conn = 0; gconn = 0
def connected(sel):
    facets = [f for p in sel for f in forb[ok_orb[p]]]
    by = {}
    for i, f in enumerate(facets):
        for t in itertools.combinations(f, 3): by.setdefault(t, []).append(i)
    adj = {i: [] for i in range(len(facets))}
    for l in by.values():
        for a in l:
            for b in l:
                if a != b: adj[a].append(b)
    comp = {}; nc = 0
    for s0 in range(len(facets)):
        if s0 in comp: continue
        comp[s0] = nc; st = [s0]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in comp: comp[y] = nc; st.append(y)
        nc += 1
    if nc == 1: return 1, 1
    # do the components form a single orbit under the group?
    idx = {f: i for i, f in enumerate(facets)}
    reach = {comp[0]}; st = [comp[0]]
    members = {}
    for i, k in comp.items(): members.setdefault(k, i)
    while st:
        k = st.pop(); f = facets[members[k]]
        for g in G:
            k2 = comp[idx[tuple(sorted(g(v) for v in f))]]
            if k2 not in reach: reach.add(k2); st.append(k2)
    return 0, int(len(reach) == nc)
def rec(p):
    global total, conn, gconn
    if p == len(inc):
        if chosen:
            total += 1; a, b = connected(chosen); conn += a; gconn += b
        return
    for take in (0, 1):
        if take:
            if any(cnt[t] + m > 2 for t, m in inc[p].items()): continue
            for t, m in inc[p].items(): cnt[t] += m
            chosen.append(p)
        if all(cnt[t] in (0, 2) for t in closing[p]): rec(p + 1)
        if take:
            for t, m in inc[p].items(): cnt[t] -= m
            chosen.pop()
rec(0)
print(f"n={n} dih={dih} c={c}: facet orbits {len(forb)} (usable {len(ok_orb)}); closed unions: {total}; strongly connected: {conn}; components form one group orbit (what ours/enum.c counts as \"closed\"): {gconn}")
