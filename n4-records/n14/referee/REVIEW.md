# Referee report: 14 monomials, 4 variables, topological degree 16

Package reviewed: `topdeg-review-N14-2026-11-27.tar.gz` (CLAIMS.md dated 2026-11-27). Review date: 2026-11-27.

## Verdicts

| Claim | Verdict | Basis |
|---|---|---|
| **P1** convex and general position, facet list, f = (14, 70, 112, 56) | **CONFIRMED** | my code, exact integers |
| **P2** all facet determinants nonzero, \|deg\| = 16, q = 0 | **CONFIRMED** | my code, exact integers, two orientations, 6136 generic u |
| **P3** 64 < 77 | **CONFIRMED** (arithmetic only) | Theorem M itself was not checked; no reference to it is given |
| **S** (U_A), (U_C), 14 vertices / 56 facets, Bihan's formula = −16 for every generic u | **CONFIRMED** | my code, exact integers |
| **S, optional** independent root count | **DONE, complete, with one caveat** | 30 positive roots, Jacobian signs +7 / −23, signed sum **−16** |
| "one more than the best published construction (15)" | **COULD NOT VERIFY** | I could not find where 15 is published; see W1 |

I did not manage to break the object. The mathematical content of P and S holds. The weak part of the
package is the novelty framing, not the computation.

Independence: none of my confirming checks imports anything from `ours/` or uses the stored `orientation`.
The stored orientation is read once, at the end of `check_P.py`, for comparison only. I read
`ours/verify_exact.py`, `mksys.py`, `drive.py`, the top of `sph.py` and the header of `sympoly.c`, and I ran
`verify_exact.py` once as a comparison (it printed `(U_A) = True (U_C) = True facets = 56`, `deg F over 42 generic u: [-16]`).

## What I ran

All code is in the tarball; `sh run_all.sh <unpacked package dir>` re-runs everything (about 2 minutes;
python 3.13, numpy, mpmath, g++ 13). Logs of my run are in `logs/`.

### P1 (`check_P.py`, exact integer arithmetic, Bareiss determinants)

- All C(14,5) = 2002 maximal minors of Â are nonzero. Smallest absolute value: 14280. So no five points
  lie on a hyperplane, with a wide margin.
- Hull by brute force: a 4-set is a facet iff the other ten points are strictly on one side (signs of 5×5
  integer determinants). Result: 56 facets, identical to the claimed list (no facet missed, none invented).
- Every one of the 14 points is a vertex. f = (14, 70, 112, 56). Every triangle is in exactly two facets,
  V − E + T − F = 0, every vertex link is a connected closed surface with χ = 2, the dual graph is connected.
- The two certificate files agree with each other (A = polytope_points, columns of C = C_int).

### P2 (`check_P.py`, exact)

- All 1001 4×4 minors of C are nonzero (smallest absolute value 20); the 56 facet determinants have
  smallest absolute value 1441.
- Orientation computed two ways:
  (a) combinatorially, by propagating the coherence rule of CLAIMS §0 across triangles (consistent on all 112 triangles);
  (b) geometrically, o_t = sign det Â_{t ∪ {i}} with i last (well defined because the other points are on one side).
  I also checked directly that (b) satisfies the coherence rule. (a) = −(b) with my starting sign; one global sign, as it must be.
- Facet signs s_t: all 56 equal (+56/−0 with (a), −56/+0 with (b)). So q = 0.
- Degree at 6136 generic integer u (out of 6204 tried; a u is used only if det(C_a, C_b, C_c, u) ≠ 0 for all
  364 triples, which is stronger than the stated genericity condition). The u include perturbations of ±e_i,
  of ±C_v, of facet barycentres, points pushed just off each side of a wall of every facet cone, and random
  u at three scales. Every one gives deg = +16 with (a), −16 with (b). Every generic u lies in **exactly 16**
  facet cones, which is what q = 0 and |deg| = 16 require.
- The stored `orientation` equals my (a), i.e. it is the **negative** of the outward geometric orientation.
  This is harmless for a claim about |deg|, but see W5.

### P3

64 < 2·56 − 70/2 − 0 = 77 holds; the bound allows |deg| ≤ 19 for this f-vector, as stated. This is only a
consistency check of numbers I recomputed. I cannot say anything about Theorem M: the package gives no
statement of its hypotheses beyond the inequality and no reference.

### S, formula part (`check_P.py`, exact)

- (U_A): 2002 minors nonzero. (U_C): 1001 minors nonzero. conv(A) simplicial, 14 vertices, 56 facets.
- I read the statements in the two preprints. arXiv:2611.34336 Theorem 1.3 and arXiv:2610.29808 Theorem A
  read exactly as transcribed in CLAIMS §0: factor (−1)^n, columns of Γ in increasing order, the extra index
  last, u outside every cone spanned by n−1 columns, hypotheses (U_A) and (U_C). I checked the convention
  against the paper's own n = 1 case (deg = [c_N > 0] − [c_1 > 0]): it gives +1 for x − 1, as it should.
- With that convention the formula gives **−16** at all 6136 generic u. It equals the degree with the
  geometric orientation (b), as it must for n = 4.

### S, independent root count (the degree without Bihan's theorem)

Method, in logarithmic coordinates y = log x, G(y) = Σ C_i exp(⟨a_i − b, y⟩), b = (57, 60, 41, 46):

1. `S_outer.py` (**exact rational arithmetic**): no root with |y|_∞ ≥ 8. The directions are covered by
   3515 dyadic boxes; on each, a set S of at most 4 dominant monomials is identified, a row vector λ with
   λ·C_i = 1 on S is computed exactly, and λ·G > 0 follows from exp(g·T0) > Σ_{i∉S}|λ·C_i|, with exp bounded
   below by a partial sum of its series. (With this crude λ the argument fails for T0 ≤ 6; that says nothing
   about roots, only about the test.)
2. `S_inner.cpp` (**interval arithmetic in doubles, outward widening after every operation**): branch and
   bound on the cube |y|_∞ ≤ 8. A box is discarded only if an interval enclosure of some λ·G excludes 0 or
   the Krawczyk image misses the box; it is accepted as "exactly one root" only if the Krawczyk image is in
   the interior of the box and the interval determinant of the Jacobian excludes 0.
   Result: 2 187 259 boxes examined, **0 unresolved**, **30 boxes with exactly one root**.
3. `S_certify.py` (**exact rational interval arithmetic**, on the original polynomials with integer exponents):
   for each of the 30 roots, Krawczyk on a box of relative radius 2⁻³⁴ around a rational point proves
   existence and uniqueness, and an interval Laplace expansion proves the sign of det Jac F. The 30 boxes
   are pairwise disjoint (exact comparison). Signs agree with step 2.

**Result: F has exactly 30 positive solutions, all nondegenerate; 7 with positive and 23 with negative
Jacobian determinant; signed sum −16.** This agrees with Bihan's formula in value and in sign. The list
of roots is in `logs/S_certify.log`; they lie in 0.67 < x_i < 1.41.

What is and is not rigorous here:
- "At least 30 nondegenerate positive roots, with these signs" is proved in exact arithmetic (step 3).
- "No root outside the cube" is proved in exact arithmetic (step 1).
- "No other root inside the cube" rests on step 2, which is floating-point interval arithmetic, not exact
  arithmetic. Its assumptions: IEEE round-to-nearest for +, −, × (each result widened by one ulp), and
  glibc `exp` accurate to within the 4 ulps by which I widen it. Overflowed enclosures are mapped to the
  whole line. I consider this sound, but it does not meet the "exact arithmetic only" standard of rule 2,
  and it is a single implementation written today and not reviewed by anyone else.
- Given completeness, deg(F, ℝ⁴_{>0}, 0) = −16 follows from the definition of the degree (0 is a regular
  value and the zero set is compact). So **|deg F| = 16 no longer depends on the two recent preprints**; it
  depends on my step 2 instead. Conversely, this is one more successful test of Bihan's formula, sign included.

## Errors and weaknesses found

**W1 (main). The benchmark "15" is not supported by anything in the package, and the package contradicts itself about it.**
CLAIMS says the best published construction (cyclic polytopes) gives (V−2)(V−4)/8 = 15 at V = 14. No
reference is given. I could not find this formula or value in Bihan's two preprints, in Ryabichev, in
Apolonskaya–Musin or in Basak–Gupta–Trivedi. Meanwhile `ours/drive.py` computes the "cyclic value" as
`(V - 1) * (V - 3) // 8`, which is **17** at V = 14, and prints "BEATS CYCLIC" only above that. With that
number, degree 16 does not beat the cyclic value. One of the two formulas is wrong for even V, or they
mean different things; either way the sentence "one more than the best published construction" cannot be
checked from what I was given. Give the exact reference (paper, theorem number) and fix or explain `drive.py`.

**W2. The cited "neighbours" solve a different problem.** Ryabichev, Apolonskaya–Musin and
Basak–Gupta–Trivedi study *simplicial* maps to the boundary of the (n+1)-simplex (five target vertices for
n = 3), not vertex maps with 14 vectors in general position. Their values at 14 vertices are far below 15:
the formula of Wang quoted in Ryabichev's appendix gives degree q^{h−j}(q+1)^j = 2·3 = **6** for the periodic
colouring of the cyclic polytope with k = 14 (n = 3, q = 2, s = 4, j = 1; my evaluation of the quoted formula).
So these papers are neither prior art for 16 nor the source of 15. The right comparison class is Bihan's
D(n, k), and there the only thing I found for n, k ≥ 3 is an upper bound.

**W3. The claim was stated as conditional on unrefereed preprints** (2610.29808 is five weeks old, 2611.34336
two days old). The root count above removes that dependence if you accept interval arithmetic; say so, or
keep the caveat.

**W4. Unverifiable side statements.** "Theorem M (ours, published and reviewed)" and "the best value known
on any 14-vertex simplicial 3-sphere is 17" come with no reference and no data. I did not check either.

**W5. Orientation wording.** The stored `orientation` is the inward one (negative of sign det Â_{Γ∪{i}}),
while the header of `sympoly.c` says "outward orientation". With the outward orientation every facet sign
is −1 and deg = −16, matching Bihan's −16. Nothing is wrong numerically, but state which convention the
certificate uses.

**W6. Small inconsistencies in the certificate.**
- `note` says "exact, 60 generic u"; `verify_exact.py` uses 43 vectors and reports 42 generic ones.
- The rounding parameters are not recorded. I reproduced `polytope_points` and `C_int` exactly from
  `source` with exponent scale 60, coefficient scale 20 and no jitter. Record this.
- The first reference appears under two titles ("Topological degree of sparse polynomial maps" on alphaXiv,
  "Topological degree of fewnomial maps" in a search-engine copy of the PDF). Cite the version you used.
- C has three zero entries: C[0][5], C[2][1] and C[2][12] (0-based, as in the JSON). (U_C) is unaffected,
  but the four polynomials do not literally have the same support: the first has 13 monomials and the third
  has 12. Bihan's setting is phrased with a common support. Mention it or perturb.
- `drive.py` calls a `climb` binary that is not in the package, and `sph.py` has a dead `move32` method.
  Neither affects the claims.

**W7. Only one implementation of the hull exists on each side.** Mine and the authors' are both brute-force
side tests. They are simple enough that I trust them, and the pseudomanifold, Euler and link checks would
catch a wrong list, but a third check with an independent tool (polymake, Qhull with exact arithmetic)
would cost little.

Things I looked for and did not find: a point not in convex position; five points on a hyperplane (nearest
miss: |minor| = 14280); a missed or invented facet; a u-dependent degree; a degree that changes between the
combinatorial and the geometric orientation; a vanishing minor of C; a facet with the wrong sign; a sign
convention in Bihan's formula different from the one transcribed.

## E1. Symmetry

- The facet list is invariant under i → i+1 in each block: **yes** (exact). The 56 facets form 8 orbits of size 7.
- Exactly symmetric configuration (before rounding): same 56 facets, degree −16 with the outward orientation
  at 40 random u, all facet signs equal, no 4×4 minor of the symmetric vectors vanishing (smallest/largest
  ≈ 1.2·10⁻³), smallest 5×5 minor 3.6·10⁻³ against a largest of 1.7. **This is 60-digit floating point
  (`E_checks.py`, mpmath), not an exact check**: the true configuration has coordinates in ℚ(cos 2π/7, sin 2π/7)
  and the base points are only given as doubles. The margins are ten orders of magnitude above the working
  precision, so I believe it, but it is evidence, not proof.
- Consequence: no jitter was needed, and the rounded object is a faithful image of the symmetric one.

## E2. Maximality

Nothing to verify. For context: Bihan's Theorem 2 in arXiv:2611.34336 gives |deg F| < ½·min(M(4,14), M(9,14))
= ½·min(77, 252), i.e. |deg F| ≤ 38 for (n, k) = (4, 9), and the paper says sharpness is unknown for n, k ≥ 3.
Your Theorem M gives 19 for this particular f-vector. 16 is consistent with both. Bihan's Theorem 3
(D(n,k) = D(k,n)) would transfer the value 16 to 9 variables and 14 monomials; I did not check that.

## E3. Prior art

I found nothing with degree ≥ 16 for 14 monomials in 4 variables, or for a vertex map on a 14-vertex
polytopal 3-sphere. The search was brief (web search, alphaXiv, the reference lists of the papers below).
Absence of a hit is weak evidence in a field where the relevant papers are weeks old.

- F. Bihan, *Topological degree of sparse polynomial maps*, https://arxiv.org/abs/2610.29808. Theorem A is the facet formula used here. No lower-bound construction for n = 4.
- F. Bihan, *Descartes' rule of signs for the topological degree*, https://arxiv.org/abs/2611.34336. Theorem 1.3 restates the formula; Theorem 2 is the upper bound quoted in E2; sharp only for n = 2 or k = 2.
- A. Ryabichev, *On vertex-minimal simplicial maps to the sphere*, https://arxiv.org/abs/2512.01137. Simplicial maps to ∂Δ^{n+1}; appendix quotes Wang's formula for the periodic colouring of cyclic polytopes.
- Q. Wang, *Exact Periodic Degrees, Degree Spectra, and Sharp Thresholds for Cyclic-Polytope Colorings*, https://zenodo.org/records/22249508. Not read beyond the abstract.
- K. Apolonskaya, O. Musin, *Minimal simplicial spherical mappings with a given degree*, https://arxiv.org/abs/2511.10870. Small-degree values of λ(n, d).
- B. Basak, R. K. Gupta, A. Trivedi, *Simplicial degree d self-maps on n-spheres*, https://arxiv.org/abs/2409.00907. Vertex counts linear in d.
- F. Bihan, F. Santos, P.-J. Spaenlehauer, *A polyhedral method for sparse systems with many positive solutions*, https://arxiv.org/abs/1804.05683. Cyclic-polytope constructions for the *number* of positive solutions, not for the degree. A possible origin of the "cyclic polytope" benchmark, but I did not find the value 15 there.

## E4. Smaller numbers

Re-rounding the symmetric configuration at smaller scales (`E_checks.py`), with every candidate verified
exactly ((U_A), (U_C), 14 vertices, formula value at 300 generic u): the smallest exponent scale that works is
26 (against your 60), and with it the smallest coefficient scale is 17 (against 20). Exponent scales 4 to 25
fail for every coefficient scale from 2 to 20 (`E4_scan.py`).

Exponents in [0, 51], coefficients in [−17, 16], same 56 facets, q = 0, Bihan's formula = −16
(file `system_N14_deg16_small.json`):

```
A = [38,48,20,32] [16,51,7,18] [0,35,24,9] [2,13,29,27] [21,0,10,28] [42,7,13,9] [50,28,31,17]
    [27,32,0,12] [21,32,31,3] [18,27,33,35] [19,22,1,30] [24,19,13,0] [29,22,40,19] [30,27,16,40]
C columns = [3,1,5,-11] [1,3,0,12] [-1,3,-5,-11] [-3,0,9,7] [-2,-2,-12,-2] [0,-3,11,-3] [2,-2,-9,8]
    [-17,3,8,6] [-13,-11,-10,-2] [1,-17,10,-2] [14,-10,-8,6] [16,4,4,-9] [7,16,0,10] [-8,15,-4,-9]
```

Limits: this is a scan over two scale factors, not an optimisation; unimodular changes of the exponent
lattice or independent perturbations might do better. I did **not** count the roots of this smaller system:
my outer-region test failed for it up to T0 = 40, which I attribute to near-degenerate facets after coarse
rounding, not investigated further. For this system only the formula value is verified.

## Not done

- Theorem M, the value 17 on a non-polytopal sphere, and the cyclic value 15 were not verified (W1, W4).
- The proofs in Bihan's preprints were not checked; only the statements were compared with CLAIMS.
- The exactly symmetric configuration was checked numerically only (E1).
- `sympoly.c`, `sph.py` and `drive.py` were read in part, not audited; they are the search, not the claim.
- No second, independent implementation of my branch and bound exists.
