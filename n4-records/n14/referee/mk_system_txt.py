#!/usr/bin/env python3
"""system.json -> plain text for S_inner (exponents shifted by the same integer vector b as S_outer.py)."""
import json, sys
s = json.load(open(sys.argv[1])); A = s['A']; C = s['C']; N = len(A)
b = [round(sum(a[k] for a in A) / N) for k in range(4)]
print(N)
for i in range(N):
    print(*[A[i][k] - b[k] for k in range(4)], *[C[r][i] for r in range(4)])
