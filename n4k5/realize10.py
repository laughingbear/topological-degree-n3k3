"""Try to realize an oriented simplicial 3-sphere as the boundary of a 4-polytope: find x in R^{V x 4} with
sigma * eps_F * det[1,x_F ; 1,x_w] > 0 for all facets F and vertices w not in F (multi-start least squares on hinge loss),
then round to integers and confirm with an exact convex hull."""
import itertools, json, sys, numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
def realize(V, ofac, tries=200, seed=0, verbose=False):
    rng = np.random.default_rng(seed)
    idx = np.array([list(G) + [w] for G, e in ofac for w in range(V) if w not in G]); sg = np.array([e for G, e in ofac for w in range(V) if w not in G], float)
    target = {tuple(G) for G, e in ofac}
    def margins(x):
        P = np.concatenate([np.ones((V, 1)), x.reshape(V, 4)], 1)[idx]          # (m,5,5)
        return np.linalg.det(P) * sg
    def loss(x, s):
        m = s * margins(x); return np.sum(np.maximum(0.05 - m, 0) ** 2) + 1e-3 * np.sum((x.reshape(V, 4).std(0) - 1) ** 2)
    for t in range(tries):
        x0 = rng.normal(size=V * 4)
        for s in (1, -1):
            r = minimize(loss, x0, args=(s,), method='L-BFGS-B', options={'maxiter': 3000})
            m = s * margins(r.x)
            if m.min() > 0:
                x = r.x.reshape(V, 4)
                for scale in (10, 20, 50, 100, 300, 1000):
                    A = np.round((x - x.mean(0)) * scale).astype(int)
                    try: h = ConvexHull(A.astype(float))
                    except Exception: continue
                    if len(h.vertices) == V and {tuple(sorted(f)) for f in h.simplices} == target and len(h.simplices) == len(target):
                        return A.tolist(), t
    return None, tries
if __name__ == '__main__':
    sub = list(itertools.combinations(range(10), 4)); SAT = json.load(open('sat8_idx.json'))
    lines = [l for l in open('pool10_sat8.txt').read().split('\n')[1:] if l.strip()]
    out = {}
    for k in json.loads(sys.argv[1]):
        z = list(map(int, lines[k].split())); ofac = [(sub[z[1 + 2 * j]], z[2 + 2 * j]) for j in range(z[0])]
        A, t = realize(10, ofac, tries=int(sys.argv[2]), seed=k)
        print(SAT[k], 'REALIZED after %d starts' % t if A else 'not found', flush=True)
        if A: out[SAT[k]] = A
    json.dump(out, open(sys.argv[3], 'w'))
