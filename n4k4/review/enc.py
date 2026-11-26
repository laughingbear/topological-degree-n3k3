"""Referee's own SAT encoder, written from the definitions (not from omsat.py).
Variables x_S (S a 4-subset, increasing): det(C_S) > 0.
(GP)  for every pair sigma and a<b<c<d outside sigma:  [s ab][s cd] - [s ac][s bd] + [s ad][s bc] = 0, brackets = determinants of
      the columns in the ORDER WRITTEN.  Encoded by truth table over the six underlying variables.
(DEG) for every vertex v:  sum_{G not containing v} eps_G * sign det(C_G) * [ -C_v in open cone(C_G) ]  >= d, where
      -C_v in cone(C_G)  <=>  the unique linear dependency of the five columns T = G + {v} (increasing t_0<..<t_4),
      sum_j (-1)^j det(C_{T - t_j}) C_{t_j} = 0, has all coefficients of one sign.
Cardinality: own one-directional unary counter (no pysat.card)."""
import itertools

def perm_sign(t):
    return (-1) ** sum(1 for i in range(len(t)) for j in range(i + 1, len(t)) if t[i] > t[j])

class Enc:
    def __init__(self, V, gp_middle_sign=-1):
        self.V = V; self.nv = 0; self.x = {}
        for S in itertools.combinations(range(V), 4): self.x[S] = self.new()
        self.gp = []
        for sig in itertools.combinations(range(V), 2):
            rest = [z for z in range(V) if z not in sig]
            for a, b, c, d in itertools.combinations(rest, 4):
                br = [sig + (a, b), sig + (c, d), sig + (a, c), sig + (b, d), sig + (a, d), sig + (b, c)]
                key = [tuple(sorted(t)) for t in br]; ps = [perm_sign(t) for t in br]
                assert len(set(key)) == 6
                for bits in itertools.product((1, -1), repeat=6):          # bits = sign of det of the SORTED subset
                    val = [bits[i] * ps[i] for i in range(6)]               # sign of the bracket as written
                    terms = (val[0] * val[1], gp_middle_sign * val[2] * val[3], val[4] * val[5])
                    if terms[0] == terms[1] == terms[2]:                    # impossible for real nonzero numbers summing to 0
                        self.gp.append([-self.x[key[i]] if bits[i] == 1 else self.x[key[i]] for i in range(6)])
    def new(self):
        self.nv += 1; return self.nv
    def formula(self, eps, d, verts=None):
        """eps: {sorted facet: +-1}. Returns (clauses, nvars) for 'degree counted at -C_v is >= d for all v'."""
        nv0 = self.nv; cl = list(self.gp)
        def AND(lits):
            y = self.new(); cl.extend([[-y, l] for l in lits]); cl.append([y] + [-l for l in lits]); return y
        def OR(lits):
            return -AND([-l for l in lits])
        for v in (range(self.V) if verts is None else verts):
            lits = []; m = 0
            for G, e in eps.items():
                if v in G: continue
                T = tuple(sorted(G + (v,)))
                coef = [self.x[T[:j] + T[j + 1:]] * (1 if j % 2 == 0 else -1) for j in range(5)]   # literal: coefficient j > 0
                inside = OR([AND(coef), AND([-l for l in coef])])
                s = self.x[G] * e                                                                  # literal: eps*sign det > 0
                lits += [AND([inside, s]), -AND([inside, -s])]; m += 1
            k = d + m                                   # need  #pos + #(not neg) >= d + m
            if k > len(lits): cl.append([]); continue
            if k <= 0: continue
            # r[i][j]: among first i literals at least j are true (only the direction r -> truth is encoded)
            prev = {0: None}                             # None = constant true ; missing = constant false
            for i, l in enumerate(lits, 1):
                cur = {0: None}
                for j in range(max(1, k - (len(lits) - i)), min(i, k) + 1):
                    r = self.new(); cur[j] = r
                    A = prev.get(j, False); B = prev.get(j - 1, False)       # r -> A or (l and B)
                    c1 = [-r, l]; c2 = [-r]
                    if A is None: continue                                    # A true: no constraint
                    if A is not False: c1.append(A); c2.append(A)
                    cl.append(c1)
                    if B is None: pass
                    elif B is False: cl.append(c2)
                    else: cl.append(c2 + [B])
                prev = cur
            cl.append([prev[k]])
        n = self.nv; self.nv = nv0
        return cl, n
