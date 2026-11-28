# Claims for cold review — a 14-monomial system in 4 variables with topological degree 16 (2026-11-27)

Everything needed is in this folder. You have no access to the project these came from. Please verify with your own code,
written from the definitions below; the scripts in `ours/` show what we did and are not to be trusted.

## 0. Definitions
- **Complex.** K is a simplicial complex given by its facets (4-element vertex sets, vertices 0..V−1). f = (V, E, T, F) counts
  vertices, edges, triangles, facets.
- **Coherent orientation.** A sign o_t = ±1 for every facet t (vertices in increasing order) such that for every triangle the two
  facets containing it induce opposite orientations: if the triangle is obtained by deleting position i from facet t and position j
  from facet t′, then o_t·(−1)^i = −o_t′·(−1)^j. Unique up to a global sign.
- **Vertex map.** A vector C_v ∈ ℝ⁴ for every vertex, such that the four vectors of every facet are linearly independent.
- **Sign of a facet.** s_t = o_t · sign det(C_a, C_b, C_c, C_d), a < b < c < d the vertices of t.
- **Degree.** For u ∈ ℝ⁴ not in the cone of any three of the vectors: deg = Σ_t s_t·[u ∈ open cone(C_a, C_b, C_c, C_d)].
  Independent of u (degree of the piecewise-linear map |K| → S³ sending each vertex to its vector, radially projected).
  Its sign flips with the global orientation; the claims are about |deg|.
- **q** = number of facets whose sign is opposite to the sign of deg.
- **Theorem M (ours, published and reviewed; sanity check only):** 4·|deg| < 2F − E/2 − q/2.
- **Published value.** The best construction published so far (cyclic polytopes) gives (V−2)(V−4)/8 for even V: 15 at V = 14.
  The best value known on any 14-vertex simplicial 3-sphere is 17 (on a sphere that is not polytopal).
- **Bihan's formula** (F. Bihan, arXiv:2611.34336, Theorem 1.3; also arXiv:2610.29808, Theorem A). A = {a_1, …, a_N} ⊂ ℤ⁴ exponent
  vectors, C a 4 × N real matrix, F(x) = Σ_i C_i x^{a_i} on the positive orthant (C_i = i-th column). Assume (U_A): all 5 × 5 minors
  of the matrix Â with columns (1, a_i) are nonzero, and (U_C): all 4 × 4 minors of C are nonzero. For u not in a cone spanned by
  3 columns of C:
  deg F = (−1)^n Σ_{facets Γ of conv(A), u ∈ cone(C_s, s ∈ Γ)} sign det C_Γ · sign det Â_{Γ ∪ {i_Γ}},
  n = 4, columns of Γ in increasing order, i_Γ any index off Γ placed last.
  These are recent preprints; the claim "deg F = ±16" below is conditional on that theorem.

## 1. Claim P: a 14-vertex simplicial 4-polytope with a vertex map of degree 16
`certs/V14_polytope_deg16.json`: `polytope_points` (14 integer points in ℝ⁴, coordinates in [0, 117]), `facets` (56), `C_int`
(one integer vector per vertex, entries in [−20, 20]), `orientation` (ours; recompute it).
Claims:
1. The 14 points are in convex position and in general position (no 5 on a hyperplane); the boundary complex of their convex hull
   is exactly `facets`; f = (14, 70, 112, 56).
2. The vertex map `C_int` has all facet determinants nonzero and |deg| = 16, with q = 0 (all 56 facets have the same sign).
3. Theorem M holds: 64 < 2·56 − 70/2 − 0 = 77.

## 2. Claim S: an explicit polynomial system with 14 monomials and topological degree 16
`certs/system_N14_deg16.json`: `A` (the same 14 points, as exponent vectors) and `C` (4 × 14, the same vectors as columns).
Claims: (U_A) and (U_C) hold; conv(A) is a simplicial polytope with 14 vertices and 56 facets; Bihan's formula gives the same value
for every generic u, namely −16 with the convention above (so |deg F| = 16).
Consequence claimed: the maximal topological degree for n = 4 variables and N = 14 monomials (k = 9) is at least 16; the value
from the published construction is 15.
Valuable if feasible: an independent count of the positive real solutions of F = 0 with the signs of their Jacobians (the signed
count should be ±16, so there are at least 16 positive solutions). Exponents go up to 117, so this may be hard; say clearly if it
was not feasible or not complete. A method that worked for a 12-monomial system of ours: branch and bound in logarithmic
coordinates with an interval (Krawczyk) test and an explicit bound on the region that can contain roots.

## 3. Lower-priority items (please comment if time allows)
- E1. Symmetry. The object was found as two orbits of Z₇ acting on ℝ⁴ by rot(2π/7) ⊕ rot(4π/7) (vertices 0..6 and 7..13, rotation
  i → i+1 in each block), with vectors equivariant for rot(2π/7) ⊕ rot(6π/7), then rounded to integers (`source` in the
  certificate holds the floating-point base points X and base vectors Y). Is the facet list invariant under the rotation? Does the
  exactly symmetric configuration (before rounding) have the same complex and the same degree? (The rounded one is what we claim.)
- E2. The degree 16 is the best found by a search over one symmetric family, not a proved maximum for N = 14. We make no
  maximality claim. The upper bounds we know for 14 vertices are far away (Theorem M allows 19 for this f-vector).
- E3. Novelty: is a value ≥ 16 for 14 monomials in 4 variables, or a vertex/simplicial map of degree ≥ 16 on a 14-vertex
  polytopal 3-sphere, already in the literature? Known neighbours: Bihan arXiv:2611.34336 and arXiv:2610.29808; Ryabichev
  arXiv:2512.01137; Apolonskaya–Musin arXiv:2511.10870; Basak–Gupta–Trivedi arXiv:2409.00907.
- E4. Smaller numbers: if you find a version with visibly smaller exponents or coefficients and the same degree, please give it.

## 4. What we would like back
A file `REVIEW.md` with, for each of P1, P2, P3 and S: CONFIRMED / REFUTED / COULD NOT VERIFY, what you actually ran, and any error
found, however small; then comments on E1–E4; plus your code. Please say plainly where a check of yours is not independent of ours.
