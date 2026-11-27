"""Overnight search: hill-climb over simplicial 3-spheres with V vertices on
   g(K) = largest d such that 'deg >= d' is satisfiable in the chirotope relaxation (omsat_fast.build_fast).
usage: walker.py V start_g log_threshold seconds seed outprefix
Self-contained except for python-sat (bundled in pydeps/).  Append-only outputs: <outprefix>.hits.jsonl, <outprefix>.log"""
import sys, json, time, random, itertools, math
from omsat_fast import build_fast
from pysat.solvers import Glucose4
V, g0, thr, secs, seed, pref = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
random.seed(seed)
def log(*a):
    with open(pref + ".log", "a") as f: f.write(time.strftime("%H:%M:%S ") + " ".join(str(x) for x in a) + "\n")
def orient(facets):
    facets = [tuple(sorted(f)) for f in facets]; ridge = {}
    for i, f in enumerate(facets):
        for k in range(4): ridge.setdefault(f[:k] + f[k+1:], []).append((i, (-1) ** k))
    o = {0: 1}; stack = [0]
    while stack:
        i = stack.pop(); f = facets[i]
        for k in range(4):
            (a, sa), (b, sb) = ridge[f[:k] + f[k+1:]]
            j, sj, si = (b, sb, sa) if a == i else (a, sa, sb)
            want = -o[i] * si * sj
            if j in o: assert o[j] == want
            else: o[j] = want; stack.append(j)
    return facets, [o[i] for i in range(len(facets))]
def random_flip(F):
    """one random bistellar 2-3 or 3-2 move (None if the random choice is not flippable)"""
    tri = {}; edges = {}
    for f in F:
        for t in itertools.combinations(f, 3): tri.setdefault(t, []).append(f)
        for e in itertools.combinations(f, 2): edges.setdefault(e, []).append(f)
    if random.random() < 0.5:
        t = random.choice(list(tri)); f1, f2 = tri[t]
        d = next(x for x in f1 if x not in t); e = next(x for x in f2 if x not in t)
        if tuple(sorted((d, e))) in edges: return None
        a, b, c = t; G = set(F); G.discard(f1); G.discard(f2)
        for p, q in ((a, b), (b, c), (a, c)): G.add(tuple(sorted((p, q, d, e))))
        return G
    cand = [e for e, fs in edges.items() if len(fs) == 3]
    if not cand: return None
    e = random.choice(cand); fs = edges[e]; others = set().union(*fs) - set(e)
    t = tuple(sorted(others))
    if len(others) != 3 or t in tri: return None
    G = set(F)
    for f in fs: G.discard(f)
    for x in e: G.add(tuple(sorted(t + (x,))))
    return G
def cyclic(V):
    F = set()
    for i in range(V):
        for j in range(V):
            f = {i, (i+1) % V, j, (j+1) % V}
            if len(f) == 4: F.add(tuple(sorted(f)))
    return F
def stellar(F, f, v):
    G = set(F); G.discard(f)
    for t in itertools.combinations(f, 3): G.add(tuple(sorted(t + (v,))))
    return G
def hval(F):
    m = {}
    for f in F:
        for e in itertools.combinations(f, 2): m[e] = m.get(e, 0) + 1
    return sum((x - 1) // 2 for x in m.values()) - len(F), len(m)
nsat = 0
def sat(F, d):
    """True / False / None (budget exhausted).  One orientation suffices: a reflection of R^4 reverses the degree."""
    global nsat; nsat += 1
    facets, o = orient(sorted(F))
    r = build_fast(V, [(f, e) for f, e in zip(facets, o)], d)
    if r[0] is None: return False
    gp, cl, X = r
    with Glucose4(bootstrap_with=gp) as s:
        for c in cl: s.add_clause(c)
        s.conf_budget(300000); return s.solve_limited()
def gval(F, gc):
    """exact g(F), searching from gc; budget-exhausted answers are treated as SAT but flagged"""
    flag = False; d = gc + 1
    r = sat(F, d)
    if r is not False:
        while r is not False:
            flag |= r is None; d += 1; r = sat(F, d)
        return d - 1, flag
    d = gc
    while d > 0:
        r = sat(F, d)
        if r is not False: return d, r is None
        d -= 1
    return 0, False
def starts():
    if V == 11:
        S10 = {tuple(f) for f in json.load(open("sphere_deg9.json"))["facets"]}
        return [cyclic(11)] * 3 + [stellar(S10, f, 10) for f in random.sample(sorted(S10), 2)]
    if V % 2 == 1: return [cyclic(V)]
    C = cyclic(V - 1); return [stellar(C, f, V - 1) for f in random.sample(sorted(C), 5)]
t0 = time.time(); S = starts(); seen = {}; hist = {}; best = 0; nlog = t0
cur = random.choice(S); gc, _ = gval(cur, g0); best = gc; log("start V", V, "g", gc, "F", len(cur))
while time.time() - t0 < secs:
    G = random_flip(cur)
    if G is None: continue
    key = frozenset(G)
    if key in seen: gg = seen[key]
    else:
        gg, flag = gval(G, gc); seen[key] = gg; hist[gg] = hist.get(gg, 0) + 1
        if len(seen) > 400000: seen.clear()
        if gg >= thr or gg > best:
            hv, E = hval(G); best = max(best, gg)
            with open(pref + ".hits.jsonl", "a") as f: f.write(json.dumps({"V": V, "g": gg, "budget_flag": flag, "F": len(G), "E": E, "h": hv, "facets": sorted(G)}) + "\n")
            log("RECORD/HIT g", gg, "flag", flag, "F", len(G), "E", E, "h", hv)
    if gg >= gc or random.random() < 0.01: cur, gc = G, gg
    if random.random() < 0.002: cur = random.choice(S); gc, _ = gval(cur, g0)
    if time.time() - nlog > 600:
        nlog = time.time(); log("progress: spheres", sum(hist.values()), "sat calls", nsat, "g histogram", sorted(hist.items()), "current g", gc)
log("DONE: spheres", sum(hist.values()), "sat calls", nsat, "g histogram", sorted(hist.items()), "best", best)
