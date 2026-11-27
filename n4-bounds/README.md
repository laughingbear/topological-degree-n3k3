# n = 4: new upper and lower bounds for every number of monomials

For polynomial maps F : ℝ⁴₊ → ℝ⁴ with N monomials (generic exponents and coefficients):

  **(N−1)(N−3)/8 ≤ max |deg F| < (3N² − 11N)/16**  (N odd; for even N the lower bound is (N−2)(N−4)/8).

- Bihan's bound (arXiv:2611.34336, Theorem 2) is |deg F| < N(N−3)/4 on the Newton-polytope side. The new upper bound is better for every N ≥ 9 (N = 9: ≤ 8 instead of ≤ 13).
- So the maximum grows like c·N² with **1/8 ≤ c ≤ 3/16**. Whether c = 1/8 is open.
- By Bihan's duality the same holds for k = 4 and any n.

**Note with full proofs (5 pages):** [note/n4bounds.pdf](note/n4bounds.pdf) ([LaTeX](note/n4bounds.tex)).

## The two theorems

Both are statements about *vertex maps*: a simplicial 3-sphere K (V vertices, E edges, F = E − V facets) with a vector C_v ∈ ℝ⁴ at each vertex, the four vectors of each facet independent. Bihan shows that the topological degree of a generic system is the degree of the vertex map on the boundary of its Newton polytope.

**Theorem M (upper bound, all sign patterns).** 4·deg < 2F − E/2 − q/2, where q is the number of negatively oriented facets. Hence |deg| < (3E − 4V)/8 ≤ (3V² − 11V)/16.

The proof uses an exact identity that holds in every dimension n:

  n·deg = Σ_v w_v − X + Y,

where w_v is the degree of the link map at the vertex v, and X, Y are the sums over the positive and negative facets of a "defect" D_t ∈ (0, 1/2): the measure of the points having exactly one negative coordinate in the basis of the facet t. For n = 4 the link degrees are then bounded by the spherical Gauss–Bonnet argument used for n = 3.

**Theorem C (lower bound).** Let V = 2m + 1 be odd, K the boundary of the cyclic polytope C(V,4), and

  C_i = (−1)^i · (cos t, sin t, cos 3t, sin 3t),  t = πi/V.

Then all facets have the same sign and |deg| = m(m−1)/2 = (V−1)(V−3)/8. The curve is Barvinok–Novik's symmetric moment curve; the construction is the 4-dimensional analogue of Bihan's equally-spaced-points construction for n = 2.

## Explicit systems

[data/cyclic_systems.json](data/cyclic_systems.json): integer systems for N = 5, …, 14 with exponents (i, i², i³, i⁴) and coefficients 10⁶·C_i rounded (even N: one extra interior monomial). Each satisfies (U_A) and (U_C) and is verified exactly with Bihan's Theorem 1.3:

| N | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|
| \|deg F\| | 1 | 1 | 3 | 3 | 6 | 6 | 10 | 10 | 15 | 15 |

Rebuild and re-verify (about 10 s): `python3 build_cyclic_systems.py`

## Checks

- `cyc_count.py`: the counting step of Theorem C, without floating point, for every odd V from 5 to 79.
- `cyc_formula.py`: same sign and degree (V−1)(V−3)/8 by direct cone counting, odd V from 5 to 61.
- `review/`: an independent review done in a separate session (report: [review/review.md](review/review.md)). **Both theorems confirmed.** The reviewer's own code (`review/referee-code/`) recomputes the cyclic degrees in 60-digit arithmetic for odd V ≤ 41, tests the identity and the bound on 251 vertex maps (random polytopes with 6–14 vertices, all our example systems, and maps from a counterexample search), and found no counterexample in 6,180 hill-climbing restarts. `review/CLAIMS-reviewed-v2.md` is the statement that was reviewed, with the reviewer's presentational fixes.
- `mixed_check.py`: our own first check of the identity. As the review points out, it tests the identity in a nearly circular way; the reviewer's `core.py` / `randmaps.py` are the real test. Kept for the record.

## Computer evidence on the constant (not proofs)

`search/` is the search run for 7 hours on 8 cores: a hill-climb over simplicial 3-spheres on the largest degree allowed by the chirotope SAT relaxation (the encoder of [n4k4/](../n4k4/README.md)). To rerun it you need python-sat, plus `omsat.py`, `omsat_fast.py` and `sphere_deg9.json` from [n4k5/](../n4k5/) in the same folder; `run.sh` as stored expects a bundled copy of python-sat and is kept as a record of the run.

| V | spheres evaluated | largest degree allowed | (V−1)(V−3)/8 |
|---|---|---|---|
| 11 | 179,995 | 10 | 10 |
| 12 | 32,810 | 12 | 12.4 |
| 13 | 11,026 | 15 | 15 |

No sphere beat the cyclic construction at V = 11 or 13. These are samples from random walks, not exhaustive searches, and no proofs were logged. An exact formula ⌊(V−1)(V−3)/8⌋ is false in general: for V = 10 the maximum is 8 on polytopes and 9 on spheres ([n4k5/](../n4k5/README.md)).

Found with Claude (Anthropic). Laughing Bear, 2026-11-27.
