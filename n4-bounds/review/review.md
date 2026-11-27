# Referee report: "Two claims on the degree of vertex maps on simplicial 3-spheres"

**Material reviewed:** `topdeg-n4-review.tar.gz`, i.e. CLAIMS.md, four scripts and four JSON example files.
**Method:** I checked every proof step by hand. I wrote independent code from the definitions in CLAIMS.md §0 (enclosed as `referee-code.tar.gz`; nothing is imported from the attached scripts). I ran the attached scripts and read them for what they actually test. I searched the literature with arXiv full text and web search.

## Overall verdict

**Theorem M: CONFIRMED.** That covers the identity 4·deg = Σ_v w_v − X + Y, the inequality 4·deg < 2F − E/2 − q/2 for every sign pattern, and the corollary |deg| < (3E − 4V)/8 ≤ (3V² − 11V)/16.

**Theorem C: CONFIRMED.** On ∂C(V,4) with V = 2m+1, the vertex map C_i = (−1)^i γ(πi/V) is a vertex map, all its facets have the same sign, and |deg| = (V−1)(V−3)/8.

I found no mathematical error and no counterexample. The open points are presentational: two orientation conventions that the proofs use but do not spell out, and two small wording issues in the statements. Each is listed under its item below, with the fix. The independent computations agree with every claim. They also show that the attached `mixed_check.py` tests the step-2 identity in a nearly circular way; see item 5.

The asymptotic gap is V²/8 ≲ max|deg| ≲ 3V²/16, a factor 3/2. Bihan proves |deg| < ½·min(M(n,N), M(k,N)) (arXiv:2611.34336, Thm 2) and says he does not know whether this is sharp for n, k ≥ 3. For n = 4, the new bound beats his for every N ≥ 9; for N ≤ 8 his bound is as good or better, because of the Gale-side term M(k,N) with k ≤ 3. Some values:

| N | Bihan: max\|deg\| ≤ | new: max\|deg\| ≤ | Theorem C attains |
|---|---|---|---|
| 9 | 13 | 8 | 6 |
| 10 | 17 | 11 | 6 (using V = 9) |
| 11 | 21 | 15 | 10 |

 That transfer to polynomial systems uses Bihan's Theorem A and was not part of this review; see item 6.

---

## Item 1. Theorem M, step 2 (the identity): **CONFIRMED**

- **Pointwise count (step 1).** Correct. Write u = Σ c_k C_k in the basis C_t. Then u ∈ L_{t,v} exactly when every c_j with j ≠ v is positive. If all four c_k > 0, all four v count. If exactly one c_k < 0, only v = k counts. If two or more are negative, none does. So Σ_v μ(L_{t,v}) = 4μ(cone t) + D_t holds for each facet separately. I also verified this numerically on every facet of all 251 test maps (several thousand facets), to 5·10⁻⁵, which is the precision of the 4-d orthant integrals.
- **μ(L_{t,v}) is a normalised solid angle.** Correct. Let P be orthogonal projection onto C_v^⊥. Then u ∈ L_{t,v} ⇔ Pu ∈ cone(PC_j : j ∈ t∖v). (⇒: project. ⇐: u − Σ c_j C_j lies in ker P = ℝC_v.) The PC_j are independent because C_t is a basis. The standard Gaussian measure is a product, so μ(L_{t,v}) is the Gaussian measure in ℝ³ of that cone, which is its normalised solid angle.
- **Signs and orientation.** Correct, provided two conventions are used, which CLAIMS.md leaves implicit:
  - lk(v) is oriented by: if t = [v, σ] carries orientation ε, then σ gets ε (put v first);
  - C_v^⊥ is oriented so that (C_v, basis) is positive in ℝ⁴.

  With these, det(C_v, C_a, C_b, C_c) = det(C_v, PC_a, PC_b, PC_c) = |C_v|·det₃(PC_a, PC_b, PC_c). So the link triangle has sign o·sign det₃ = s_t, and the link orientation is coherent (standard). The link map is a vertex map on a 2-sphere. The integral formula one dimension lower then gives Σ_{t∋v} s_t μ(L_{t,v}) = w_v. Nothing here uses the facet signs being equal, so mixed signs are covered.
  - *Numerical check:* `core.py` computes w_v as a separate 3-dimensional problem, with exact integer arithmetic in C_v^⊥, induced orientation and projected vectors. `randmaps.py` checks Σ_{t∋v} s_t μ(L_{t,v}) = w_v with closed-form solid angles. Maximum error over all vertices of 251 maps: 5·10⁻¹⁵.
- **Σ_t s_t μ(cone t) = deg.** Correct. The map f is injective on each facet: two points of the linear simplex that are positive multiples of each other must coincide, because both have barycentric sum 1. Almost every u avoids the triangle cones, and for such u the signed count equals deg. Integrating over u gives the formula. Verified numerically to 8·10⁻⁶.
- **Suggested fix (presentation only).** State the two orientation conventions above explicitly. As written, "a facet t ∋ v gives the triangle t∖v with the same sign s_t" reads like a definition. It is really a consequence of those conventions, and it is what makes w_v an integer, the degree of a genuine map.

## Item 2. Theorem M, step 3, 0 < D_t < 1/2: **CONFIRMED**

The map u ↦ −u sends the 4 one-negative orthants onto the 4 three-negative orthants, so the two families both have measure D_t. The all-positive orthant has measure μ(cone t) > 0, and so does the all-negative one. Hence 2D_t ≤ 1 − 2μ(cone t) < 1. D_t > 0 because each open orthant has positive measure.

The constant 1/2 is sharp but never attained. Take Cov(c) = aaᵀ + εI with a = (1, 1, 1, −3). As ε → 0, D_t → 1/2: in the limit c is proportional to ±a, which has exactly one or exactly three negative coordinates. The search (item 5c) and the random runs reached D_t = 0.49995, and never 1/2. Step 6 needs only X ≥ 0 with X > 0 when q = 0, and Y < q/2 when q > 0. Both follow.

## Item 3. Theorem M, steps 4–6: **CONFIRMED**

- **Gauss–Bonnet, Σ_v w_v = Σ_e w_e − (p − q).**
  - *Area side.* 4π·w_v = Σ s·area is the integral formula on S² ⊂ C_v^⊥, and the area of a spherical triangle is its angle sum minus π (Girard).
  - *Angle side.* The angle at the image of link vertex a, in the triangle (a, b, c), is the angle between the projections of C_b and C_c onto span(C_v, C_a)^⊥. Around a, these consecutive angles, each multiplied by s_t, add up to 2π·w_{va}. This is the winding number of the projected link polygon of the edge va, with the induced orientation.
  - *Edge winding is well defined.* Viewed from v, the plane span(C_v, C_a)^⊥ is oriented by (C_v, C_a, ·); viewed from a, by (C_a, C_v, ·). The link cycle's orientation also reverses. The two flips cancel, so w_e does not depend on which endpoint is used. CLAIMS.md uses this without saying so. The cleanest statement is 2π·w_e = Σ_{t∋e} s_t·∠_t(e), which is manifestly symmetric.
  - *Summing.* Σ_v Σ_{e∋v} = 2Σ_e, Σ_v p_v = 4p and Σ_v q_v = 4q, which gives the relation.
  - *Numerical check:* the per-vertex relation 4w_v = 2Σ_{e∋v} w_e − (p_v − q_v) holds exactly, with w_v exact and w_e from angle sums, at every vertex of all 251 maps. The global relation does too.
- **Edge bound, w_e ≤ (p_e − 1)/2.** Correct in all cases. Every angle lies in (0, π), so 2πw_e < πp_e whenever p_e ≥ 1, because the negative facets only subtract. An integer w_e < p_e/2 satisfies w_e ≤ (p_e − 1)/2.
  - *p_e = 0* (all facets at e negative): 2πw_e < 0, so w_e ≤ −1 < −1/2. Fine.
  - *Mixed signs at e:* covered by the same inequality.
  - *Data:* the bound held on every edge of every map. Both edge types are common in the random families: in a 30-map sample, 108 of 1044 edges had p_e = 0 and 860 had facets of both signs.
- **Arithmetic (step 6).** Correct: Σ_e p_e = 6p, so Σ_e w_e ≤ 3p − E/2. Hence 4deg ≤ 2p + q − E/2 − X + Y < 2p + q − E/2 + q/2 = 2F − E/2 − q/2.
  - *Strictness:* if q = 0 it comes from X > 0 (since p = F > 0); if q > 0, from Y < q/2.
  - *Corollary:* drop q, use F = E − V to get 4|deg| < 3E/2 − 2V, i.e. |deg| < (3E − 4V)/8. Then E ≤ V(V−1)/2 gives (3V² − 11V)/16. For deg < 0, reversing the orientation negates deg and swaps p and q. Correct.
- **Wording.**
  1. The hypothesis "deg ≥ 0" in Theorem M is not used. The inequality holds for any coherent orientation, with q counted in that orientation; it is most useful when deg ≥ 0.
  2. In §0, "(V−1)(V−3)/8 ≤ max ≤ (3V²−11V)/16" should read "max < (3V² − 11V)/16", and the lower bound holds only for odd V. (For even V, the odd case V−1 gives a lower bound by monotonicity, provided adding a vertex cannot lower the maximum. That needs a short argument, such as a stellar subdivision of a negative facet, and is not in the file.)
- **Same-sign refinement and the Remark.** Both confirmed. With q = 0, step 5 gives w_e ≤ ⌊(m_e − 1)/2⌋, so 4deg < Σ_e⌊(m_e−1)/2⌋ − F. On ∂C(V,4) I checked the edge degrees for V = 9, 11, 21: the V edges {i, i+1} have m_e = V−2, the V edges {i, i+2} have m_e = 3, and all the others have m_e = 4. The refinement then gives deg < V(V−3)/8, so the gap to Theorem C is (V−3)/8.
  - A consequence the authors may want to state: for V ≤ 11, Theorem C is optimal among same-sign maps on ∂C(V,4). For example, V = 11 gives deg < 11, so deg ≤ 10 = (V−1)(V−3)/8.

## Item 4. Theorem C: **CONFIRMED**

- **Determinant lemma.** Correct. (cos t, sin t, cos 3t, sin 3t) is an invertible linear image of (e^{−3it}, e^{−it}, e^{it}, e^{3it}) = e^{−3it}(1, z, z², z³) with z = e^{2it}. Each factor z_b − z_a = 2i·e^{i(y_a+y_b)}·sin(y_b − y_a). Every y_a appears in 3 of the 6 pairs, so the exponential factors multiply to e^{3iΣy}, which cancels ∏e^{−3iy_a}. That leaves det = c·∏sin(y_b − y_a) with c a nonzero constant times (2i)⁶ = −64. c is real because both sides are real.
- **Vertex map.** Correct. 0 ≤ x_a < x_b < π, so every sine is positive and every 4×4 minor of C is nonzero.
- **Facet list.** Correct. It is the facet list of the cyclic 4-polytope, i.e. the Gale-evenness sets of size 4. Both families together give (V−3)(V−2)/2 + (V−3) = V(V−3)/2 facets. I generated the facets independently from Gale's evenness condition and from a Qhull hull of V points on the trigonometric moment curve; both lists coincide with the {i, i+1, j, j+1} recipe for every odd V from 5 to 41.
- **Signs.** Correct.
  - sign det C_f = sign(c)·(−1)^{Σf}: non-wrap facets have Σf = 2i+2j+2, even; wrap facets {0, j, j+1, V−1} have Σf = 2j+V, odd.
  - Orientation: for the moment-curve realisation, the 5×5 matrix is a Vandermonde matrix in the order (f_1, …, f_4, z), with sign (−1)^{#{k : f_k > z}}. For a non-wrap facet, z < i gives 4, i+1 < z < j gives 2, z > j+1 gives 0, so the orientation is +. For a wrap facet, 1 ≤ z ≤ j−1 gives 3 and j+2 ≤ z ≤ V−2 gives 1, so it is −. This is the outward-normal orientation, hence coherent.
  - So s_f = sign(c) for every facet. My code reproduces this with its own coherent-orientation routine and 60-digit determinants.
- **Cramer sign computation.** Correct. The formula (−1)^k c_k = ∏_{j<k}… · ∏_{j>k}… follows from the lemma. For each column, the (−1)^j factors of the other columns cancel between numerator and denominator, leaving (−1)^k.
  - With 0 = x_0 < τ < x_1: sin(τ − x_j) < 0 ⇔ j ≥ 1, and sin(x_j − τ) < 0 ⇔ j = 0.
  - *Non-wrap facets with i ≥ 1:* sorted positions r = 1..4 give parities (i, i, j, j). Counted iff i and j are both even, with 2 ≤ i < j ≤ 2m−2. That gives C(m−1, 2) facets.
  - *i = 0:* c_1 < 0, so not counted.
  - *Wrap facets:* signs (+, (−1)^j, (−1)^j, +). Counted iff j is even, j ∈ {2, …, 2m−2}: m−1 facets.
  - Total: C(m−1, 2) + (m−1) = C(m, 2) = (V−1)(V−3)/8.
- **u = γ(τ) is a regular value.** Correct. u lies in the span of C_a, C_b, C_c iff det(C_a, C_b, C_c, u) = 0, iff sin(τ − x_i) = 0 for some i ∈ {a, b, c}. That does not happen for τ ∈ (0, x_1). It also keeps u off the boundary of every facet cone.
- **Observation for the write-up.** γ(t+π) = −γ(t), so C_i = γ(2π(m+1)i/V). The vertex vectors are therefore V equally spaced points on the Barvinok–Novik symmetric moment curve SM₄, with cyclic-polytope vertex i assigned to point (m+1)i mod V, i.e. i·2⁻¹ mod V. This is the 4-dimensional analogue of Bihan's n = 2 construction C_r = (cos 2πjr/N, sin 2πjr/N) (arXiv:2611.34336, proof of Thm 2′). Saying so would make the construction less mysterious.

## Item 5. Independent computation

All code is in `referee-code.tar.gz`; logs are in `logs/`.

**(a) Cyclic construction (`cyclic.py`): CONFIRMED.** For every odd V from 5 to 41, all facets have the same sign and |deg| = (V−1)(V−3)/8: 1, 3, 6, 10, …, 171, 190. This is direct cone counting at 5 random u and at 3 values u = γ(τ), all agreeing, in 60-digit arithmetic, with a facet list from Gale evenness and my own orientation. The smallest |det| and the smallest Cramer coefficient were asserted to be far above working precision.

**(b) Random maps and examples (`randmaps.py`): CONFIRMED, no failures.**
- **Random maps.** 210 vertex maps on boundaries of random simplicial 4-polytopes with V = 6..14 (every value represented, with 17 to 28 maps each). The polytopes have integer coordinates, with simpliciality checked exactly. Three families: uniform random integer C; small multiples of the polytope coordinates plus noise; rounded points on random curves (cos kt, sin kt, cos lt, sin lt). 209 of the 210 have facets of both signs.
- **Example files.** All 20 maps: the nonpolytopal degree-9 sphere (its supplied orientation is coherent and agrees with mine), the degree-8 system, the D45 family with degrees 0–8 (each claimed degree reproduced exactly), and the nine degree-6 systems.
- **Record maps.** The 21 record maps from (c), with degrees up to 15.
- **Results, over all 251 maps:**

  | check | result |
  |---|---|
  | F = E − V | always |
  | deg, w_v | exact, integer arithmetic |
  | Σ_v w_v = Σ_e w_e − (p − q), and per vertex | exact, always |
  | w_e ≤ (p_e − 1)/2 | every edge |
  | Σ_{t∋v} s_t μ(L_{t,v}) = w_v | error ≤ 5·10⁻¹⁵ |
  | Σ s_t μ(cone t) = deg | error ≤ 8·10⁻⁶ |
  | per-facet decomposition | error ≤ 5·10⁻⁵ |
  | identity 4deg = Σw_v − X + Y | error ≤ 1.4·10⁻⁵ |
  | 0 < D_t < 1/2 | always; max D_t 0.49995 |
  | 4deg < 2F − E/2 − q/2 | always; smallest slack 4.0 among random maps |

- **Weakness.** Random maps mostly have degree 0–2. That is why the search maps in (c) were added to this test.

**(c) Counterexample search (`search.py`): no counterexample.**
- **Setup.** Hill-climbing with three temperatures (0, 0.5, 1.5) on g = 4deg − (2F − E/2 − q/2). Moves: Gaussian perturbation of one C_v, or the sign flip C_v → −C_v. 21 complexes: C(V,4) for V = 7..13, and 14 random polytopes with V = 6..12. On odd cyclic polytopes half the restarts were seeded near the Theorem C map. 6180 restarts of 1500 steps each.
- **Results.**
  - The best value of g was −4.5, on a V = 6 sphere (deg 1). Largest degrees found: 15 on C(13,4), equal to Theorem C; 7 on random 12-vertex spheres; 6 on a random 10-vertex sphere with q = 3.
  - Every record map was re-verified with exact arithmetic and passed the full set of checks in (b).
  - The same-sign refinement was never violated.
  - On odd C(V,4) the climber never exceeded the Theorem C degree (3, 6, 10, 15 for V = 7, 9, 11, 13). Those restarts were partly seeded from Theorem C, so this is weak evidence about optimality.
- **D_t probe.** Over 300 random Gaussian 4×4 bases, the largest D_t was 0.4999 (two runs: 0.499918 and 0.498693). In the explicit family Cov(c) = aaᵀ + εI with a = (1, 1, 1, −3), D_t is 0.420, 0.474, 0.492, 0.497 and 0.4997 for ε = 10⁻¹, 10⁻², 10⁻³, 10⁻⁴ and 10⁻⁶. So 1/2 is the supremum and is never reached, as the proof says. A Nelder–Mead polishing stage was too slow at this Genz tolerance, so I report the random-basis and family results only.

**The attached scripts.** I ran all three. Each prints its OK line (mixed_check: "ALL OK"; cyc_count: count = C(m,2) for odd V ≤ 79 and every interval for τ; cyc_formula: odd V ≤ 61). They do not all test what CLAIMS.md says they test:

- **`mixed_check.py` checks the step-2 identity in a nearly circular way.** It computes w_v as the pointwise signed count Σ_{t∋v} s_t[u ∈ L_{t,v}], i.e. via the same set L_{t,v} that step 1 uses, rather than as the degree of the link map in C_v^⊥. It then estimates X and Y by Monte Carlo on common samples.
  - By the step-1 pointwise identity, Σ_t s_t[exactly one c_k < 0](u) = Σ_v w_v − 4deg holds exactly for every generic sample u. So the Monte Carlo value of X − Y is exactly the integer Σ_v w_v − 4deg. That is why the log shows agreement to "0.000" with 4·10⁵ samples.
  - The check therefore confirms step 1 and the u-independence of the counts. It does not test the two substantive parts of step 2: that μ(L_{t,v}) is the solid angle of the projected cone, and that the signed sum is the degree of the link map with the induced orientation.
  - What it does test genuinely: the Gauss–Bonnet relation, since w_e by pointwise counting is a true winding count; D_t < 1/2 by Monte Carlo; and the bound.
  - Other weaknesses: its 30 random maps have V = 6..13, not 6..14 (`rng.integers(6, 14)` excludes 14); 26 of them have degree 0 and 4 have degree 1, so the bound is barely exercised; singular facets are skipped silently via `LinAlgError`; and the Monte Carlo volume Σ s_t μ(cone t) is printed but never asserted.
- **`cyc_count.py` checks the counting, not the Cramer formula.** Its docstring says it compares the sign rule with "exact rational Cramer signs". In fact `product_sign` reads the signs off the product formula of the proof, so it assumes the lemma and the Cramer formula and checks only the combinatorial count. It does show something more than the proof: the count is C(m,2) for τ in every interval (x_l, x_{l+1}), not only for l = 0.
- **`cyc_formula.py` tests what it says, with two caveats.** It uses float64, which is adequate up to V = 61 because the smallest sine is about π/V. It takes the facet list from the recipe without checking that the recipe gives the cyclic polytope.
- **`cyc.py`** is the search that found the construction. It is context only and was not needed.

## Item 6. Facts used without proof: all **standard and correct**

- **F = E − V.** Euler's formula V − E + T − F = 0 together with 2T = 4F (each triangle lies in 2 facets) gives T = 2F and F = E − V. This holds for every closed triangulated 3-manifold; Dehn–Sommerville is not needed. Verified on every complex.
- **Integral formula for the degree with mixed signs.** Standard: the degree is the signed preimage count at almost every regular value, then integrate over u. It needs f to be a homeomorphism of each facet onto Δ_t with local orientation s_t, which holds (item 1). The same applies one dimension down for w_v.
- **Girard's formula**, area = angle sum − π: classical.
- **Winding number as a sum of signed angles in (−π, π):** classical.
- **Coherence of the induced link orientation**, and the edge winding being independent of the endpoint: standard, but they should be stated (items 1 and 3).
- **E ≤ V(V−1)/2:** trivial.
- **Scope.** The proof of Theorem M never uses polytopality, and it would go through verbatim on any closed oriented combinatorial 3-manifold. The JSON example of a degree-9 map on a nonpolytopal sphere is consistent with that.
- **Not reviewed: the transfer to polynomial systems.** Bihan's theorem that the topological degree of the system equals the degree of the vertex map on ∂conv(A), under (H) or genericity, is quoted as context. I checked that his Theorem A (arXiv:2610.29808) says this, up to a global sign. For Theorem C to transfer, one also needs a cyclic-polytope support A with all maximal minors of (1, a_i) nonzero, which is easy, and the C above, all of whose maximal minors are nonzero (item 4).

## Item 7. Novelty

**Searched:** the full texts of Bihan arXiv:2611.34336 and 2610.29808 (both v1/v2 dated 25 Nov 2026), Ryabichev arXiv:2512.01137, and Basak–Gupta–Trivedi arXiv:2409.00907; Barvinok–Novik and Barvinok–Lee–Novik (symmetric moment curve), via web search; plus an arXiv topical search for degrees of PL or simplicial maps of spheres.

- **The identity n·deg = Σ_v w_v − X + Y.** Not found in any of these. Bihan works with signed counts at one direction u (and ±u for his half-count), not with averages over u, and has no link-degree identity. Ryabichev and BGT study simplicial maps to ∂Δ^{n+1}, where degree is a preimage count; there are no solid angles. The ingredients are classical: the angle-sum relations of Gram, Sommerville and Perles–Shephard type, and averaging the degree over u. I cannot rule out an equivalent identity in the literature on angle sums or on PL Gauss–Bonnet for maps; that literature was not searched exhaustively.
- **The bound |deg| < (3E − 4V)/8.** Not found. On the Newton side it is strictly better than Bihan's |deg| < F/2 = (E − V)/2, which he proves with a one-direction counting argument; he says he does not know whether his bound is sharp for n, k ≥ 3. For n = 4, Theorem M shows his bound is not sharp for any N ≥ 9; at N = 9, ≤ 13 becomes ≤ 8. At N = 8 his Gale-side bound (≤ 5) is still the better one.
- **The construction with (cos t, sin t, cos 3t, sin 3t) on cyclic polytopes.** Not found. The curve is Barvinok–Novik's symmetric moment curve, used there for bicyclic polytopes and face counts, not for degrees. The construction is the natural analogue of Bihan's n = 2 equally-spaced-points construction. The closest combinatorial relative I found is the periodic colouring φ_k of the cyclic polytope, cited in Ryabichev's appendix (Q. Wang's preprint, which I could not read); it gives degrees of order k²/25 for simplicial maps in dimension 3, a different setting. Novelty claims for the construction should cite Barvinok–Novik for the curve and Bihan for the n = 2 pattern.
- **Limitations.** The two Bihan papers are dated the day before this review, so I could not check citing literature. Wang's preprint was not available to me. MathSciNet/zbMATH were not searched.

## What I could not check

- The 4-d orthant probabilities (μ(cone t) and D_t) are numerical, from Genz's algorithm at about 10⁻⁶. The identity is therefore verified to about 10⁻⁵ per map, not exactly; everything else is exact or integer-rounded with large margin.
- My random maps have small degrees. The high-degree regime is covered only by the search maps (degrees up to 15, V ≤ 13) and the cyclic family.
- The polynomial-system context (Bihan's theorems) was taken as given.
