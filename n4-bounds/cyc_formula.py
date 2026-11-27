"""Cyclic-polytope vertex maps: C_i = (-1)^i * g(pi*i/V), g(t) = (cos t, sin t, cos 3t, sin 3t), on the boundary of C(V,4), V odd.
Checks: every facet has the same sign, and deg = (V-1)(V-3)/8, by direct cone counting at random u."""
import numpy as np
from mixed_check import orient
from cyc import cyc_facets
rng = np.random.default_rng(5)
for V in range(5, 62, 2):
    facets, o = orient(cyc_facets(V)); i = np.arange(V); t = np.pi * i / V
    C = (-1.0) ** i * np.array([np.cos(t), np.sin(t), np.cos(3*t), np.sin(3*t)])
    M = np.array([C[:, list(f)] for f in facets]); d = np.linalg.det(M); s = np.sign(d) * np.array(o)
    assert abs(s.sum()) == len(facets), "mixed signs"
    Mi = np.linalg.inv(M); degs = {int(abs((s * ((Mi @ rng.standard_normal(4)) > 0).all(1)).sum())) for _ in range(20)}
    m = (V - 1) // 2
    assert degs == {m * (m - 1) // 2}, (V, degs)
    print(V, "facets", len(facets), "all same sign, deg =", m*(m-1)//2, " upper bound on this sphere: deg <", V*(V-3)/8, " general bound <", (3*V*V-11*V)/16)
