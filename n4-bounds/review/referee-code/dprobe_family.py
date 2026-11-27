"""D_t for the explicit family Cov(c) = a a^T + eps I, a = (1,1,1,-3): D_t -> 1/2 as eps -> 0, never equal."""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import numpy as np
from search import D_of
a = np.array([1, 1, 1, -3.0])
for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
    L = np.linalg.cholesky(np.outer(a, a) + eps * np.eye(4))
    print(f"eps={eps:g}: D = {D_of(np.linalg.inv(L).ravel()):.6f}", flush=True)
