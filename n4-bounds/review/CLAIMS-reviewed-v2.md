# Two theorems on the degree of vertex maps on simplicial 3-spheres

**Cold-reviewed 2026-11-26 (Fable, fresh chat): both theorems CONFIRMED** (`topdeg-n4-coldreview.md`). This version (v2) includes the reviewer's presentational fixes.

## 0. Setting and definitions

- K: a simplicial 3-sphere (every 3-simplex is a "facet"), with V vertices, E edges, F facets. Euler and Dehn–Sommerville give F = E − V.
- A **vertex map** is a choice of a vector C_v ∈ ℝ⁴ for each vertex such that the four vectors of every facet are linearly independent. It defines a continuous map f: |K| → S³: on each facet extend linearly and project radially. The facet t is mapped onto the spherical simplex Δ_t spanned by its four normalised vectors, which is the unit sphere intersected with the open cone cone(t) = {Σ c_k C_k : c_k > 0}, up to boundary.
- Fix an orientation of K: a sign o_t for each facet with its vertices sorted. The **sign** of a facet is s_t = o_t · sign det(C_{t1}, …, C_{t4}). p and q are the numbers of facets with s_t = +1 and −1.
- **deg** is the Brouwer degree of f. For u outside the cones of all triangles, deg = Σ_{t : u ∈ cone(t)} s_t.
- **Link degree w_v.** For a vertex v, project the vectors of its neighbours to the hyperplane C_v^⊥ ≅ ℝ³. This is a vertex map on the 2-sphere lk(v) (one dimension lower); w_v is its degree, with these conventions: lk(v) is oriented by putting v first (if t = [v, σ] has orientation ε then σ gets ε), and C_v^⊥ is oriented so that (C_v, basis) is positive in ℝ⁴. Then det(C_v, C_a, C_b, C_c) = |C_v|·det₃ of the projections, so the triangle t∖v has sign s_t and the link orientation is coherent.
- **Edge winding w_e.** For an edge e = ab, project the remaining vertices of the facets around e to span(C_a, C_b)^⊥ ≅ ℝ². The link of e is a cycle; w_e is the winding number of the resulting closed polygon around 0; equivalently 2π·w_e = Σ_{t∋e} s_t·(dihedral angle of t at e), which shows that w_e does not depend on the endpoint from which e is viewed. p_e is the number of positive facets containing e, m_e the number of all facets containing e.
- μ is the rotation-invariant probability measure on the unit sphere S³ (equivalently the standard Gaussian measure on ℝ⁴, since all sets below are cones).

Why it matters (context only, not part of the review): for polynomial systems of n = 4 equations in 4 variables with N monomials, F. Bihan (arXiv:2611.34336) shows that the topological degree is the degree of such a vertex map on the boundary of the Newton polytope (V ≤ N), and proves |deg| < F/2 ≤ V(V−3)/4. The theorems below give (V−1)(V−3)/8 ≤ max for odd V and max < (3V²−11V)/16 for all V. The new upper bound beats Bihan's for N ≥ 9 (N = 9: ≤ 8 instead of ≤ 13); for N ≤ 8 his Gale-side bound is as good or better. The curve γ below is the Barvinok–Novik symmetric moment curve, and the construction is the 4-dimensional analogue of Bihan's n = 2 construction with equally spaced points on a circle.

## 1. Theorem M (upper bound, no assumption on signs)

**Claim.** For every vertex map on a simplicial 3-sphere and either orientation (q counted in that orientation; the statement is useful when deg ≥ 0):

    4·deg < 2F − E/2 − q/2.

Consequently, for every vertex map, |deg| < (3E − 4V)/8 ≤ (3V² − 11V)/16.

**Proof.**

*Step 1 (pointwise count).* Let t be a facet and u a point with u = Σ_{k∈t} c_k C_k, all c_k ≠ 0. For v ∈ t let L_{t,v} = {u : c_j > 0 for all j ∈ t, j ≠ v}. Then

    #{v ∈ t : u ∈ L_{t,v}} = 4·[all c_k > 0] + [exactly one c_k < 0].

(If all c_k > 0, all four v count. If exactly one c_k < 0, only v = k counts. If two or more are negative, none.)

*Step 2 (integrate).* Multiply by s_t, sum over t and integrate against μ:

    Σ_t s_t Σ_{v∈t} μ(L_{t,v}) = 4·Σ_t s_t μ(cone(t)) + Σ_t s_t D_t,

where D_t := μ{exactly one coordinate of u in the basis C_t is negative}.
- Σ_t s_t μ(cone(t)) = deg (integrate the regular-value formula for the degree over u).
- L_{t,v} = ℝ·C_v + cone(the other three vectors of t). Its μ-measure equals the normalised measure, in C_v^⊥, of the cone spanned by the projections of the other three vectors. Summing over t ∋ v with the signs s_t gives, by the same integral formula one dimension lower, the degree w_v of the link map.

Hence the **identity**

    4·deg = Σ_v w_v − X + Y,    X = Σ_{s_t = +1} D_t,   Y = Σ_{s_t = −1} D_t.

*Step 3 (0 < D_t < 1/2).* The 16 open orthants of the basis C_t partition ℝ⁴ up to a null set. The map u ↦ −u sends orthants with exactly one negative coordinate to orthants with exactly three, so these two families have the same measure D_t. The all-positive orthant has positive measure. So 2·D_t < 1, and D_t > 0.

*Step 4 (links; spherical Gauss–Bonnet).* For a vertex map on a simplicial 2-sphere (here: lk(v), with p_v positive and q_v negative triangles), the signed area formula gives 4π·w_v = Σ_{triangles} s·(angle sum − π). The signed angle sum at a link vertex (an edge e ∋ v of K) is 2π·w_e. Hence 4·w_v = 2·Σ_{e∋v} w_e − (p_v − q_v). Summing over v (each edge has 2 endpoints, each facet 4 vertices):

    Σ_v w_v = Σ_e w_e − (p − q).

*Step 5 (edges).* 2π·w_e = (sum of the angles of the positive facets at e) − (sum for the negative facets), and every angle lies in (0, π). So w_e < p_e/2, hence w_e ≤ (p_e − 1)/2 as w_e is an integer. Since Σ_e p_e = 6p,

    Σ_e w_e ≤ 3p − E/2.

*Step 6.* 4·deg = Σ_e w_e − p + q − X + Y ≤ 2p + q − E/2 − X + Y < 2p + q − E/2 + q/2 = 2F − E/2 − q/2, using X > 0, Y < q/2 and p = F − q.

For the consequence: if deg < 0 reverse the orientation. Then use q ≥ 0, F = E − V and E ≤ V(V−1)/2. ∎

**Remark (same-sign refinement).** If q = 0, step 5 gives 4·deg < Σ_e ⌊(m_e − 1)/2⌋ − F.

## 2. Theorem C (lower bound: cyclic polytopes)

Let V = 2m + 1 ≥ 5 be odd. K = boundary complex of the cyclic polytope C(V,4): vertices 0, …, V−1 and facets {i, i+1, j, j+1} (indices mod V) with four distinct elements; F = V(V−3)/2. Let

    γ(t) = (cos t, sin t, cos 3t, sin 3t),   x_i = πi/V,   C_i = (−1)^i · γ(x_i).

**Claim.** This is a vertex map, all facets have the same sign, and |deg| = m(m−1)/2 = (V−1)(V−3)/8.

**Proof.**

*Lemma.* det(γ(y_1), γ(y_2), γ(y_3), γ(y_4)) = c·∏_{a<b} sin(y_b − y_a) for a nonzero constant c.
(γ(t) is an invertible linear image of e^{−3it}(1, z, z², z³) with z = e^{2it}; the Vandermonde determinant ∏_{a<b}(z_b − z_a) with z_b − z_a = 2i·e^{i(y_a+y_b)} sin(y_b − y_a); the exponential factors cancel against ∏ e^{−3iy_a}.)

*Vertex map.* For a < b the difference x_b − x_a lies in (0, π), so all 4×4 minors of C are nonzero.

*Signs.* Split the facets into non-wrap facets {i, i+1, j, j+1} with 0 ≤ i, i+2 ≤ j, j+1 ≤ V−1, and wrap facets {0, j, j+1, V−1} with 1 ≤ j ≤ V−3.
- For a sorted facet f_1 < … < f_4: sign det(C_f) = sign(c)·(−1)^{f_1+f_2+f_3+f_4}. This is sign(c) for non-wrap facets (sum 2i+2j+2) and −sign(c) for wrap facets (sum 2j+V, V odd).
- Orientation: realise K on the moment curve a_i = (τ_i, τ_i², τ_i³, τ_i⁴), τ_0 < … < τ_{V−1}. For a sorted facet and any other vertex z, sign det((1,a_{f_1}), …, (1,a_{f_4}), (1,a_z)) = (−1)^{#{k : f_k > z}} is independent of z and defines a coherent orientation o_f. It is + for non-wrap facets and − for wrap facets.
- So s_f = sign(c) for every facet.

*Degree.* All signs agree, so |deg| = #{facets f : u ∈ cone(f)} for any u outside the cones of triangles. Take u = γ(τ) with 0 < τ < x_1. (By the Lemma, u is in the span of three of the C_i only if sin(τ − x_i) = 0 for some i.)
For a sorted facet f write u = Σ_k c_k C_k. By Cramer's rule and the Lemma,

    (−1)^k · c_k = ∏_{j∈f, j<k} sin(τ − x_j)/sin(x_k − x_j) · ∏_{j∈f, j>k} sin(x_j − τ)/sin(x_j − x_k).

All denominators are positive. With 0 = x_0 < τ < x_1: sin(τ − x_j) < 0 iff j ≥ 1, and sin(x_j − τ) < 0 iff j = 0.
- *Non-wrap facet with i ≥ 1.* For the vertex in sorted position r (r = 1..4), exactly r − 1 numerator factors are negative. So sign c_k = (−1)^{k + r − 1}: the four signs are (−1)^i, (−1)^{i+1+1}, (−1)^{j+2}, (−1)^{j+1+3}, i.e. (−1)^i, (−1)^i, (−1)^j, (−1)^j. All positive iff i and j are both even. Number of such facets: pairs of even numbers 2 ≤ i < j ≤ V−3 = 2m−2 (then j ≥ i+2 automatically): C(m−1, 2).
- *Non-wrap facet with i = 0* (f = {0, 1, j, j+1}). k = 0: no negative factor, c_0 > 0. k = 1: the factor sin(τ − x_0) > 0, no negative factor, so sign c_1 = (−1)^1 < 0. Not counted.
- *Wrap facet {0, j, j+1, V−1}.* k = 0: c_0 > 0. k = j: one factor j' < k (namely 0, positive), none negative among j' > k: sign (−1)^j. k = j+1: factors below are 0 (positive) and j (negative): sign (−1)^{j+1+1} = (−1)^j. k = V−1: factors below are 0 (+), j (−), j+1 (−): sign (−1)^{V−1} = +. All positive iff j is even: j ∈ {2, 4, …, 2m−2}, that is m − 1 facets.

Total: C(m−1, 2) + (m − 1) = C(m, 2). ∎

**Remark.** On this K the same-sign refinement of Theorem M gives deg < V(V−3)/8 (V edges of degree V−2, V edges of degree 3, all others of degree 4), so the construction is within (V−3)/8 of optimal for ∂C(V,4).

## 3. Files

- `mixed_check.py` — Theorem M numerically: computes deg, all w_v, all w_e, D_t by Monte Carlo with common samples; asserts the identity of step 2, the relation of step 4, D_t < 1/2 and the bound, on the examples in `topdeg4/` (a degree-9 map on a 10-vertex sphere, a degree-8 polynomial system with two negative facets, a family with degrees 0..8, nine degree-6 systems) and on 30 random maps on random polytopes. For the files that give a system (A, C), K is the boundary of conv(A) (rows of A) and C has one column per vertex.
- `cyc_formula.py` — Theorem C in floating point: same sign and deg = (V−1)(V−3)/8 for odd V = 5..61.
- `cyc_count.py` — Theorem C combinatorially: the Cramer sign rule and the count C(m,2) for every position of τ, odd V = 5..79.
- `cyc.py` — the scan over frequency pairs that found the construction (context).

Run: `python3 mixed_check.py; python3 cyc_formula.py; python3 cyc_count.py` (numpy, scipy).
