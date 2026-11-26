# Cold review: D(4,4) = {0,…,6} (record `topdeg-D44-2026-11-26.md`)

Reviewer: Claude (Fable), fresh session, 2026-11-26. Stance: assume the claim is wrong and try to break it.
All "own code" below was written from the definitions and from Bihan's text, not from the author's files.

## Verdict

**I could not break the claim.** Every item checks out; the corrections below are about wording, attribution and missing control code, not about the mathematics or the computations.

| # | Item | Verdict | Basis |
|---|---|---|---|
| 1 | Reduction deg F = ± deg(K,C) for n = 4 | **CONFIRMED** (conditional on Bihan's two preprints, see 1.4) | Statements read in Bihan's paper; signs corroborated by root counts on 7 systems |
| 2 | Lower bound, 0..6 attained | **CONFIRMED** | Own integer-arithmetic verifier; own two-route root search |
| 3 | Upper bound, mathematics (a–e) | **CONFIRMED** | Re-derived each step |
| 4 | Upper bound, implementation | **CONFIRMED** | Line-by-line read; own encoder; identical results sphere by sphere; 771 own controls |
| 5 | Completeness of the sphere list | **CONFIRMED** | Author's list equals Lutz's published list, class by class |

Still open (already known to the author): no DRAT certificate has been checked. Mitigation found here: two independent encodings × four solvers all agree.

---

## 1. Reduction

**1.1 Statements checked against arXiv:2611.34336v1 (full text), not the paraphrase.**
- Lemma 1.1 (= [D, Lemmas 3.1, 3.2]): under (H), for a triangulation T of ∂P with vertices in A, deg F is the degree of Π_T/|Π_T| : ∂P → S^{n−1}; for u outside every cone on n−1 columns C_s (s a vertex of T), a simplex τ meets the preimage of u iff u ∈ cone(C_τ), in one interior point, with local degree (−1)^n · sign det C_τ · sign det Â_{τ∪{i_τ}}.
- Theorem 1.3 (= [D, Theorem A]): under (U_A) and (U_C), deg F = (−1)^n Σ_{Γ facet, u ∈ cone(C_Γ)} sign det C_Γ · ε_Γ. For n = 4 the prefactor is +1.
- §1.1: (U_A) and (U_C) imply (H). Theorem 2 gives |deg F| < ½·M(4,9) = 13.5, i.e. ≤ 13, as the record says (stated explicitly at the end of Bihan's §5).
- D(n,k) is defined with exactly N = n+k+1 monomials and nonzero maximal minors of both matrices; real exponents allowed. The record uses the same definition.

**1.2 Interior points.** Under (U_A) no five points of A are on a hyperplane, so each facet of P contains exactly four points of A, all vertices. Theorem 1.3 therefore only ever sees columns of hull vertices. The citation of [D, Thm C] is correct (it is the vertex reduction, Prop. 1.5 of arXiv:2610.29808v2) but not needed here.

**1.3 Fewer than 9 hull vertices.** P is then a simplicial 4-polytope with 5 ≤ V ≤ 8 vertices, its boundary is a simplicial 3-sphere on V vertices, and the V columns still satisfy (U_C). The SAT runs for V = 5..8 cover this. The genericity of u = −C_v also holds with respect to all 9 columns (same argument as 3a).

**1.4 What the result rests on.** The reduction is taken from two preprints (arXiv:2610.29808v2 and 2611.34336v1, both dated November 2026, unrefereed, and by their own acknowledgements AI-assisted). I read the proof chain (Prop. 1.3 → Lemma 3.1 → Lemma 3.2 in [D]) and found nothing wrong: it is the standard "degree = degree of the boundary map, homotope to the PL map under (H)" argument. Independent evidence: for all seven systems in `D44_family.json` the sum of Jacobian signs over the positive roots I found equals the formula's value **including the sign** (0, −1, …, −6). The paper should say plainly that the theorem is conditional on Bihan's Theorem A.

**1.5 ε is a coherent orientation.** Checked in exact arithmetic for all 16 polytopes: Bihan's ε_Γ satisfies the boundary-orientation rule on every ridge.

## 2. Lower bound

Own verifier `exact.py`: Python integers only, Laplace-expansion determinants, hull facets found by exact side tests. No import from `verify_exact.py`.

| C (label) | (U_A) | (U_C) | hull vertices | facets | deg F (≈57–60 generic u, and at u = −C_v for all 9 v) | facet signs +/− |
|---|---|---|---|---|---|---|
| 0 | yes | yes | 9 | 22 | 0 | 8/14 |
| 1 | yes | yes | 9 | 22 | −1 | 6/16 |
| 2 | yes | yes | 9 | 22 | −2 | 6/16 |
| 3 | yes | yes | 9 | 22 | −3 | 5/17 |
| 4 | yes | yes | 9 | 22 | −4 | 4/18 |
| 5 | yes | yes | 9 | 22 | −5 | 3/19 |
| 6 (§2 of record) | yes | yes | 9 | 22 | **−6** | 2/20 |

- Every 3-subset of a facet lies in exactly two facets (simplicial 3-sphere, 22 facets, f-vector consistent).
- The nine systems in `deg6_systems.json`: all have (U_A), (U_C), 9 hull vertices, 22–27 facets, |deg| = 6 exactly; they are nine pairwise different combinatorial types (author indices 15, 19, 60, 101, 211, 244, 353, 702, 895); type 211 is neighborly with all 27 facets negative. All as stated.

**Roots, own code (`myroots.py`).** Two routes: (i) Powell-hybrid in log coordinates from 6,000 wide random starts, polished by my own 60-digit Newton; (ii) the Gale dual system on the polyhedron Δ with hit-and-run starts inside Δ. Degree-6 system: **8 positive roots, Jacobian signs seven − and one +, sum −6**; relative residuals ≤ 2·10⁻⁵¹; the Gale route found no further root. The roots agree with the record, including the one near x ≈ (4.99·10⁴, 1.54·10⁻⁴, 118, 1.90). For labels 0..5 I found 0, 1, 4, 5, 8, 7 roots with sign sums 0, −1, −2, −3, −4, −5. No completeness claim for the root lists (none is needed: the sign sums match the exact degree).

## 3. Upper bound, mathematics

**a. Genericity of u = −C_v.** If −C_v = Σ λ_i C_{s_i} with λ ≥ 0 over three columns, then either v ∉ {s_i} and four columns are dependent, or v ∈ {s_i} and three columns are dependent (the coefficient of C_v becomes 1 + λ > 0). Both contradict (U_C). The same argument with four columns shows no facet containing v has −C_v in its cone. Correct.

**b. Cramer.** −C_v = Σ λ_i C_{g_i} has λ_i = −det(C_G with column i replaced by C_v)/det(C_G), so λ_i > 0 iff the replaced determinant has the opposite sign. Under (U_C) no λ_i is 0, so the cone condition is open. Correct. (My encoder uses the equivalent circuit form instead: the dependency Σ_j (−1)^j det(C_{T∖t_j}) C_{t_j} = 0 of the five columns T = G ∪ {v} has all coefficients of one sign.)

**c. Grassmann–Plücker.** For fixed σ = (x,y) and any a,b,c,d, [σab][σcd] − [σac][σbd] + [σad][σbc] = 0 with brackets as determinants in the order written (it is the 2×2 Plücker relation in R⁴/span σ). The three terms are nonzero under (U_C), so they cannot share a sign. The minus sign belongs on the middle term; unsorted tuples must be reduced with the sign of the sorting permutation, which `chi()`/`parity()` in `omsat.py` does. Sensitivity test: with a plus on the middle term my formula rejects the true chirotope of 60 out of 60 random real matrices, so a sign error here would not go unnoticed.

**d. Reflection.** For R with det R < 0, RC satisfies (U_C), cone membership is preserved, every determinant changes sign, so deg(K,RC) = −deg(K,C). Hence "no C with deg ≥ d" implies "no C with |deg| ≥ d". The remark that C → −C does not help for n = 4 is right (determinants keep their sign, and the degree does not depend on u).

**e. Relaxation.** For real C with (U_C) and deg(K,C) ≥ d, set x_S = [det C_S > 0]. Then (GP) holds by c; for every v the facet count at u = −C_v equals deg(K,C) by a (the count is independent of generic u because (K,C) is a map from a closed oriented pseudomanifold to R⁴∖0: independent columns keep every facet image off the origin); the auxiliary variables are defined by equivalences, so they extend uniquely; the cardinality constraint is then true. So the formula is satisfiable whenever such a C exists. Nothing in the formula constrains x beyond GP and the degree count. **Sound.** The formula also admits non-realizable chirotopes, so SAT answers prove nothing, as the record says.

## 4. Upper bound, implementation

**4.1 Reading `omsat.py` (build) and `orient.py`.**
- `chi(t)`: literal for det > 0 of the ordered tuple, via sorted subset and inversion parity. Correct.
- GP loop: 36 pairs × 35 quadruples = 1,260 relations, 16 clauses each; the `neg` flag is on the middle term only. Correct.
- `e ↔ (a xor b)`, `y ↔ ∧e_i`, `p ↔ y ∧ s`, `n ↔ y ∧ ¬s`: all four clause sets are full equivalences. `s` carries ε. Correct.
- Count: #p − #n ≥ d rewritten as #p + #(¬n) ≥ d + m over 2m literals; the empty clause when d > m is right.
- `orient.py`: ε_g = −ε_f·(−1)^i·(−1)^j across a shared triangle is the correct coherence rule for increasing-order simplices; consistency is asserted; all facets are reached.
- No defect found.

**4.2 Own encoder (`enc.py`, `sph.py`, `run.py`).** Differences from the author's: spheres taken from Lutz's file, not from `spheres.py`; orientation by propagating ordered tetrahedra; GP clauses from a truth table; cone test in circuit form; own unary counter (no `pysat.card`).

| V | spheres | bound | own encoder (CaDiCaL 1.5.3) | record |
|---|---|---|---|---|
| 5 | 1 | ≥2 / ≥1 | UNSAT / SAT | same |
| 6 | 2 | ≥2 / ≥1 | UNSAT all / SAT both | same |
| 7 | 5 | ≥4 | UNSAT all 5 | same |
| 7 | 5 | ≥3 | SAT for exactly 1 (Lutz 3_7_4 = author 0) | same |
| 8 | 39 | ≥5 | UNSAT all 39 | same |
| 8 | 39 | ≥4 | SAT for exactly 2 (Lutz 3_8_21, 3_8_32 = author 7, 20) | same |
| 9 | 1,296 | **≥7** | **UNSAT all 1,296** | same |
| 9 | 1,296 | ≥6 | SAT for 24 | same 24 |
| 9 | 1,296 | ≥5 | SAT for 597 | same 597 |

- **No sphere differs.** I also ran the author's `omsat.py` and matched the SAT sets through canonical forms: identical sets for (9,7), (9,6), (9,5), (8,4), (7,3).
- "≥7" on my encoding is also UNSAT for all 1,343 spheres (V = 5..9) with Glucose 4 and MiniSat 2.2, and for all 1,296 nine-vertex spheres with Lingeling. The author's `crosscheck.py` reproduces: 1,343 instances, none SAT, with both Glucose and MiniSat; max 2,655 variables and 27,099 clauses as stated.
- All nine `deg6_systems.json` types and the family's type 702 are among the 24 spheres where "≥6" is SAT, as they must be.

**4.3 Own positive controls** (χ fixed to the true signs; expect SAT at d = |deg|, UNSAT at d+1):
- 600 random Gaussian C on random spheres (200 each for V = 7, 8, 9): 0 mismatches. Degrees reached: 0–2 only.
- 155 hill-climbed C (own climber) to exercise higher degrees: 0 mismatches. Degrees up to 4 at V = 8 and up to 5 at V = 9.
- 16 exact systems (nine degree-6 systems and the family 0..6), chirotope from exact integers: all SAT at |deg|, UNSAT at |deg|+1.
- The author's `control.py 150` reproduces: 450 controls, 0 mismatches.

## 5. Completeness of the sphere list

- **Is "the counts match" sufficient?** Logically yes, given three things: every listed complex is a 3-sphere (flips preserve the PL type), the listed complexes are pairwise non-isomorphic (canonical forms), and the literature count is right. Flip-connectivity is then not even needed. But it leans on the count and on nauty.
- **Direct comparison, done here.** I downloaded Lutz's `3_manifolds_all.txt` and `3_manifolds_all_type.txt` (Manifold Page, TU Berlin mirror): 1,344 complexes with 5–9 vertices, of which 1, 2, 5, 39, 1,296 are typed S³ and one 9-vertex complex is the twisted S²-bundle. With my own nauty-free canonical form (minimum over all facet/ordering seeds of a breadth-first relabelling), the author's `spheres.json` and Lutz's S³ lists are **identical as sets of isomorphism classes for every V = 5..9** (0 only-Lutz, 0 only-author). The non-sphere is not in the author's list.
- Extra checks on every complex: each triangle in two facets, every vertex link a 2-sphere, Euler characteristic 0, orientable. Neighborly counts 1, 4, 50 for V = 7, 8, 9 match the literature (50 neighborly among the 1,296: Altshuler–Steinberg, as quoted in Datta's survey arXiv:math/0701735, Thm 5.16, which also records the proper-bistellar connectivity the enumeration uses).
- Remaining trust: Lutz's list itself (it agrees with Altshuler–Steinberg's earlier classification). My own SAT runs used Lutz's file directly, so the upper bound no longer depends on `spheres.py` at all.

## Novelty

Nothing found on D(4,4) or on the maximal degree for n = k = 4 (alphaXiv full-text search, web search). Bihan's paper (v1, 25 Nov 2026) says the maximal degree is unknown for n, k ≥ 3, calls n = k = 3 the first open case, and gives only the bound 13 for n = k = 4. The companion paper arXiv:2610.29808v2 contains no such computation. The paper is one day old, so check again before posting.

## Statements stronger than the evidence

1. **"V ≤ 8 reproduces Bihan's known D(4,1), D(4,2), D(4,3)."** D(4,1) and D(4,2) are Bihan's (Theorem 2′). D(4,3) is not: Bihan states the maximum is unknown for n, k ≥ 3. D(4,3) = {0..4} is the author's own n = 3 theorem combined with Bihan's Theorem 3. The V = 8 run is still a real cross-check, of the author's Gauss–Bonnet proof by an unrelated method.
2. **"Bihan names n = k = 4 as the next open case."** He names n = k = 3 as the first open case and mentions n = k = 4 only through the value 13.
3. **Controls.** `control.py` as supplied contains only the random controls, and those only reach degree ≤ 2 (the author's own histogram: 0, 1, and two cases of 2). The control on "the nine exact degree-6 systems" is not in the supplied code. Random controls alone do not support "the encoding computes the degree correctly" near the bound; the degree-6 controls do, and I reproduced them independently.
4. **"Upper bound … proved."** It is computer-assisted with unchecked UNSAT answers, and conditional on Bihan's Theorem A (unrefereed). Both caveats belong next to the claim, not only in "not done yet".
5. Not checked by me and not needed for the claim: "~230M steps", "11–14 spheres", "1,064 of 1,142 polytopal types", "at V = 9 gluing known pieces gives only 5 or 4". (Two 3M-step runs of `search.c` gave best 6, on 5–6 spheres.)

## Section 5 (gluing) and Section 6

- §5: the additivity argument is correct at the vertex-map level (the two removed facets contribute +1 and −1 at any u in their common cone, and the orientation reversal makes the union coherent). Rerunning the author's `glue4.py` on the family's degree-6 system gives −12, −18, −24 on 14, 19, 24 vertices with all 4-subsets independent. I did not re-implement it. The step to "6m ∈ D(4, 5m−1)" still needs explicit polytopes with (U_A), as the record says.
- §6: I see no error. The inequality needs each spherical tetrahedron inside an open hemisphere (true for independent columns), the same-sign hypothesis, and (U_C) for the link maps (inherited). It gives ≤ 8 at V = 9, consistent with 6.

## Required corrections

1. Fix the attribution of D(4,3) and the "next open case" sentence (items 1, 2 above).
2. Add the degree-6 control code to `control.py` (or a new file) and add higher-degree controls; report the degree histogram of the controls.
3. Replace "counts match the literature" by the direct comparison with Lutz's files; cite the file names and the Datta survey for flip-connectivity.
4. State the theorem as conditional on Bihan's Theorem A / Theorem 1.3, with the root counts as corroboration.
5. State next to the claim that the UNSAT results are not yet certificate-checked; list the solvers and that two independent encodings agree.
6. In §3 point 3 add "nonzero under (U_C)" to the three-term argument; in §1 note that under (U_A) the interior points drop out directly from Theorem 1.3.
7. Before publishing: DRAT/LRAT certificates for the 1,343 "≥7" instances (unchanged from the author's list).

## Reviewer's files

`exact.py` (exact verifier), `myroots.py` (roots, two routes), `sph.py` (Lutz parser, canonical form, orientation, manifold checks), `cmp_lists.py` (list comparison), `enc.py` (encoder), `run.py` (SAT runs), `controls.py`, `controls2.py` (controls). External data: Lutz, `3_manifolds_all.txt`, `3_manifolds_all_type.txt`.
