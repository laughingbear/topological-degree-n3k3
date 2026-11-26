# The maximal topological degree for n = k = 3

**Main result (version 3):** for every k ≥ 1,

  **D(3,k) = {0, 1, …, ⌊(3k+4)/4⌋}**, and, by Bihan's H-duality, **D(n,3) = {0, 1, …, ⌊(3n+4)/4⌋}** for every n ≥ 1.

So the maximal topological degree for n = 3 with N monomials is exactly ⌊(3N−8)/4⌋.
- The upper bound is a Gauss–Bonnet argument.
- The lower bound comes from base cases N = 5…8 plus a 7-vertex "gadget" glued on repeatedly. Degree is additive under gluing along a triangle, and each gadget adds 4 monomials and raises the degree by 3.
- Explicit systems attaining the bound for N = 5…16, all exactly verified, are in [data/sharp_family.json](data/sharp_family.json). Check them with `python3 verify/verify_family.py`.
- The gadget is [data/gadget.json](data/gadget.json); check it with `python3 verify/gadget_exact.py data/gadget.json`.

History: D(3,3) = {0,1,2,3} (v1), then the N = 8, 9 addendum (v2), then all N (v3), all on 2026-11-26.

This answers the first open case in F. Bihan, *Descartes' rule of signs for the topological degree* ([arXiv:2611.34336](https://arxiv.org/abs/2611.34336), §6). The same question appeared earlier in [arXiv:2610.29808](https://arxiv.org/abs/2610.29808) v1, Oct 21, 2026.

**Note (3 pages):** [note/D33.pdf](note/D33.pdf) ([LaTeX source](note/D33.tex)).

## Statements

**Upper bound.** Let F : ℝ³₊ → ℝ³ be a fewnomial map with N monomials, where all maximal minors of the exponent matrix and of the coefficient matrix are nonzero. Let V be the number of vertices of the Newton polytope, and p, q the numbers of its facets whose image triangles are positively, respectively negatively, oriented. Then

  4·|deg F| ≤ max(p, q) + V − 4 ≤ 3V − 8 ≤ 3N − 8.

Consequences:
- For N = 7 this gives |deg F| ≤ 3, which improves Bihan's Theorem C (|deg F| ≤ 4).
- For every N ≥ 7 it is strictly stronger than Theorem C.
- Through Bihan's H-duality (Theorem D), D(n,3) ⊆ {0, …, ⌊(3n+4)/4⌋}.

The proof is a signed spherical Gauss–Bonnet count on Bihan's piecewise-linear model of the degree. Each vertex winding w_v satisfies 2·w_v ≤ (number of positive triangles at v) − 1, because every triangle angle is less than π.

**Lower bound.** An explicit system with n = 3, N = 7 and deg F = 3. Exponents:

```
a = (13,11,18) (25,10,1) (18,18,8) (5,18,1) (0,10,6) (5,0,0) (17,2,6)
```

Coefficients:

```
C = [ -2   2  -7   9   3  -4    2 ]
    [ -6   2  -1   3   9  -3    2 ]
    [  7  10  -7   3  -4  -9  -10 ]
```

The system is in [data/system_deg3.json](data/system_deg3.json). Systems with the same exponents and |deg F| = 0, 1, 2 are in [data/small_examples.json](data/small_examples.json). So every value in {0, 1, 2, 3} has an explicit, exactly verified example.

## Addendum (same day): the bound is attained for N = 8 and N = 9

The system in [data/system_N8_deg4.json](data/system_N8_deg4.json) has n = 3 and N = 8 monomials, with (U_A) and (U_C). Its degree is **deg F = −4**:
- The exact Theorem A computation gives −4.
- The signed-area formula gives −4, with all 12 triangles negatively oriented.
- Numerically there are 6 positive zeros, each Krawczyk-certified (one at x₁ ≈ 1.33·10⁴), with Jacobian signs +,−,−,−,−,−.

So **D(3,4) = D(3,5) = {0, …, 4}**, and by H-duality D(4,3) = D(5,3) = {0, …, 4}. The bound ⌊(3N−8)/4⌋ is therefore attained for every N ≤ 9. The first open case is N = 10, where it allows |deg F| = 5.

The boundary complex is the 8-vertex triangulation with valences (5,5,5,5,5,5,3,3). The vertex map is a 4-sheeted branched cover with vertex windings −2 at the six valence-5 vertices and −1 at the two valence-3 vertices, which is exactly the equality case of the proof. It was found by enumerating all 223,440 labelled triangulations of S² on 8 vertices and, for random vertex images, maximising the degree over all of them ([verify/search_n8.py](verify/search_n8.py)).

```sh
python3 verify/verify_exact.py data/system_N8_deg4.json
python3 verify/indep_verify.py data/system_N8_deg4.json
python3 verify/roots_krawczyk.py data/system_N8_deg4.json 20000
```

## Verification

```sh
pip install numpy scipy sympy mpmath

# exact (rational arithmetic, stdlib only): (U_A), (U_C), facets, Theorem A for several generic u
python3 verify/verify_exact.py data/system_deg3.json

# independent exact verifier (sympy), plus the signed-area cross-check
python3 verify/indep_verify.py data/system_deg3.json

# numerical zeros: 5 positive zeros, each certified by an interval Krawczyk test; Jacobian signs sum to 3
python3 verify/roots_krawczyk.py data/system_deg3.json 3000
```

Expected:
- **Theorem A:** `deg F = 3` for every u.
- **Signed area:** `3.000000000000`, with all 10 triangles positive.
- **Zeros:** 5 positive zeros with Jacobian signs +,+,+,−,+.

The exact Theorem A computation is the proof. The zero computation is independent corroboration; completeness of the zero list is not claimed.

## How it was found

Bihan's Lemma 2.1 identifies deg F with the degree of the vertex map ∂P → S², a_s ↦ C_s. Riemann–Hurwitz and the angle bound showed that a degree-3 map must be a 3-sheeted branched cover on the 7-vertex triangulation with valences (6,5,5,5,3,3,3). [verify/search_hexagon.py](verify/search_hexagon.py) finds one in seconds. The exponents are a realization of that triangulation, perturbed to satisfy (U_A).

## Credits

Laughing Bear, working with Claude (Anthropic). The search, the proof and the verification code were produced with Claude. An independent re-verification was done in a separate session, with separately written code (`indep_verify.py`, `roots_krawczyk.py`).

MIT License.
