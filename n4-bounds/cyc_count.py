"""Theorem C, checked two ways without floating point in the decision:
(1) the sign rule  sign c_k = (-1)^(k + pos(k) + g + [k > l])  for u = gamma(tau), tau in (x_l, x_{l+1}),
    is compared with exact rational Cramer signs computed from the product formula;
(2) the number of facets with all c_k > 0 is C(m,2) for EVERY l in 0..V-2, V = 2m+1.
"""
from math import comb
def facets(V):
    F = set()
    for i in range(V):
        for j in range(V):
            f = {i, (i+1) % V, j, (j+1) % V}
            if len(f) == 4: F.add(tuple(sorted(f)))
    return sorted(F)
def rule_all_positive(f, l):
    g = sum(1 for k in f if k > l)
    return all((k + pos + g + (k > l)) % 2 == 0 for pos, k in enumerate(f, 1))
def product_sign(f, k, l):
    """sign of c_k from c_k*(-1)^k = prod_{j<k} sin(tau-x_j)/sin(x_k-x_j) * prod_{j>k} sin(x_j-tau)/sin(x_j-x_k),
    all x in [0,pi), x_l < tau < x_{l+1}: only the numerators can be negative."""
    s = (-1) ** k
    for j in f:
        if j < k and j > l: s = -s        # sin(tau - x_j) < 0
        if j > k and j <= l: s = -s       # sin(x_j - tau) < 0
    return s
for V in range(5, 80, 2):
    m = (V - 1) // 2; Fs = facets(V); assert len(Fs) == V * (V - 3) // 2
    for l in range(0, V - 1):
        cnt = 0
        for f in Fs:
            a = all(product_sign(f, k, l) > 0 for k in f); assert a == rule_all_positive(f, l)
            cnt += a
        assert cnt == comb(m, 2), (V, l, cnt)
print("OK: count = C(m,2) for all odd V in 5..79 and every l")
