"""Referee's helper routines. Written from CLAIMS.md only; exact integer arithmetic.
No code and no stored signs from ours/ are used here."""
import json
from fractions import Fraction
from itertools import combinations


def load(path):
    with open(path) as fh:
        return json.load(fh)


def det(M):
    """Exact determinant of a square integer/rational matrix (fraction-free Bareiss)."""
    n = len(M)
    M = [list(r) for r in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j])
                if isinstance(M[i][j], int):
                    M[i][j] //= prev
                else:
                    M[i][j] /= prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def sgn(x):
    return (x > 0) - (x < 0)


def det_cols(cols):
    """det of the matrix whose COLUMNS are the given vectors."""
    n = len(cols)
    return det([[cols[j][i] for j in range(n)] for i in range(n)])


def lifted(P, idx):
    """det of the matrix with columns (1, P[i]) for i in idx (in the given order)."""
    return det_cols([[1] + list(P[i]) for i in idx])
