# Cold review: D(4,5) = D(5,4) = {0,…,8} (record `topdeg4/v10/topdeg-D45-2026-11-26.md`)

Reviewer: Claude (Fable 5.2), fresh session, 2026-11-26. Stance: assume the claims are wrong and try to break them.
All "own code" below was written from the definitions and from Bihan's text. The author's scripts were read; only `spheres.py`/`enum10.py` (to regenerate the list) and `roots_trop.py` (to reproduce a number) were run.

## Verdict

**I could not break claim A or claim B.** The one item the record flags as its weak point (root count 14 / sign sum 6) does not reproduce: the bundled system has exactly 10 positive roots with sign sum +8. The corrections below are about a stale statement, attribution and wording.

| # | Item | Verdict | Basis |
|---|---|---|---|
| 1a | Degree-8 system and family 0..8: (U_A), (U_C), simplicial hull, degree | **CONFIRMED** | Own integer-arithmetic verifier |
| 1b | Root count | **CONFIRMED: 10 roots, sign sum +8.** The record's "14 roots, sum 6" is **FLAWED** (not reproducible) | Own branch-and-bound with Krawczyk test, no unresolved box; author's own script also gives 10 / +8 |
| 1c | Own degree-8 systems on other types | **CONFIRMED** (not required, done anyway) | Degree 8 on all 22 polytopal types; complete root count matches on 9 of them |
| 2a | Sphere 101272: 3-sphere, stated properties | **CONFIRMED** | Own checks; it is Lutz's `manifold_3_10_249002` |
| 2b | Reduction to chirotopes, global sign, 3-term GP | **CONFIRMED** | Re-derived |
| 2c | Exactly 4 compatible chirotopes | **CONFIRMED** | Own encoder, five solvers; same four as the author's |
| 2d | Four biquadratic final polynomials | **CONFIRMED** | Rows regenerated from the chirotope; exact integer check |
| 2e | Already known to be non-polytopal? | **CONFIRMED: yes, known** | Not among Firsching's 162,004 polytopal spheres |
| 2f | Validation numbers of the rank-5 test | **CONFIRMED** | 37/2 at V = 8; 1,143/153 at V = 9 |
| 3 | Degree-9 vertex map on sphere 101272 | **CONFIRMED** | Exact, 200 random u and all ten u = −C_v |
| 4a | Sphere list for V = 10 | **CONFIRMED** | Equal to Lutz's list class by class |
| 4b | Sweep "≥ 8", "≥ 9", "≥ 10" | **CONFIRMED** | Own encoder on all 247,882 spheres of Lutz's list; identical to the author's file sphere by sphere |
| 4c | Positive controls at V = 10 | **CONFIRMED** | 598 own controls, 0 mismatches, degrees 0 to 9 |
| 4d | 22 polytopal / 16 non-polytopal | **CONFIRMED** | Exact hulls; own rank-5 test; Firsching's list |
| 5 | Logic of claim A | **CONFIRMED** (conditional and computer-assisted, see 5) | Re-derived |

Still open: no proof certificate (DRAT/LRAT) for any V = 10 UNSAT answer; my root search uses ordinary floating point with explicit slack, not outward-rounded intervals; Bihan's two preprints are unrefereed.

---

## 1. The degree-8 system

**1a. Exact check.** Own verifier `exact.py` (Python integers, Bareiss determinants, hull facets by exact side tests, Bihan's Theorem 1.3 with ε_Γ as he defines it). I read Theorem 1.3, Lemma 1.1, Theorem 2, Theorem 2′, Theorem 3 and the definition of D(n,k) in arXiv:2611.34336v1; the record quotes them correctly.

| C (label) | (U_A) | (U_C) | hull vertices | facets | deg F (54–60 random generic u, and at u = −C_v for all ten v) | facet signs +/− |
|---|---|---|---|---|---|---|
| 0 | yes | yes | 10 | 29 | 0 | 9/20 |
| 1 | yes | yes | 10 | 29 | 1 | 20/9 |
| 2 | yes | yes | 10 | 29 | −2 | 10/19 |
| 3 | yes | yes | 10 | 29 | 3 | 22/7 |
| 4 | yes | yes | 10 | 29 | 4 | 21/8 |
| 5 | yes | yes | 10 | 29 | 5 | 24/5 |
| 6 | yes | yes | 10 | 29 | 6 | 25/4 |
| 7 | yes | yes | 10 | 29 | 7 | 26/3 |
| 8 (§5 of record) | yes | yes | 10 | 29 | **8** | 27/2 |

The hull is a simplicial 3-sphere of the type of author index 82118 = Lutz `manifold_3_10_248466`. Label 8 equals `system_deg8_N10_small.json`.

**1b. Root count: settled, and it agrees with the formula.**

Method (`bnb.cpp`, own). Work in y = log x. Branch and bound on the cube [−R, R]⁴ by bisection. A box is discarded only by a test that proves it has no root:
- gap test: at a root the largest monomial is at most 6κ times the fifth largest, where κ = max |(C_T⁻¹C_j)_r| over all 4-sets T (computed exactly; κ = 144 here). Proof: with T the four largest monomials, m_T = −C_T⁻¹ Σ_{j∉T} m_j C_j.
- form test: some combination g = wC of the equations has Σ g_i m_i of one sign on the box, by the natural interval extension or the mean-value form. The w used are the 4 rows, the 120 combinations that eliminate three monomials, and the rows of J(y₀)⁻¹.
- Krawczyk operator image disjoint from the box.

A box is closed as "one root" only when the Krawczyk image of the doubled box lies in its interior (existence and uniqueness there). Boxes are bisected until one of these happens.

Completeness of the region: the gap inequality with (U_A) puts every root in the union of 1,260 explicit parallelepipeds; this gives |log x_k| ≤ 390.5 for this system. I used R = 394.4 and, as a second run, R = 1000 split into three slabs.

Result for `system_deg8_N10_small.json`: 7,316,903 boxes, **0 unresolved, 10 roots**, and the same 10 in the second run.

| sign det J | x₁ | x₂ | x₃ | x₄ |
|---|---|---|---|---|
| + | 61.568781 | 1.2987929 | 44.069953 | 0.066049984 |
| + | 4.8623292 | 19.354202 | 0.037323859 | 1.4019541 |
| + | 1.048859 | 0.88457213 | 1.9991066 | 0.20463232 |
| + | 1.3190655 | 0.86007838 | 0.88234585 | 1.1697061 |
| + | 1.8152427 | 0.13761143 | 0.14442304 | 1.2706193 |
| + | 0.45613304 | 1.128726 | 2.040299 | 2.6227826 |
| + | 0.42321851 | 2.1397338 | 1.649363 | 1.5245866 |
| − | 0.94701853 | 1.012609 | 1.0753688 | 0.88909773 |
| + | 0.88654779 | 1.0447774 | 1.092965 | 0.86213811 |
| + | 0.67029418 | 0.52682446 | 1.4433182 | 0.42625089 |

**Sum of Jacobian signs = 9 − 1 = +8 = deg F**, sign included. Each root was polished at 60 digits (residual below 10⁻⁵⁸); the normalised Jacobian determinant is at least 3·10⁻⁴ in absolute value, so no sign is in doubt.

Checks on the method:
- On the reviewed N = 9 system it returns exactly the 8 known roots, sum −6, including the one at x ≈ (4.99·10⁴, 1.5·10⁻⁴, 118, 1.9).
- On all nine family systems (R computed per system, 420 to 549): roots / sign sum = 0/0, 3/1, 4/−2, 9/3, 6/4, 9/5, 12/6, 9/7, 10/8. Every sum equals the exact degree, sign included.
- Limit of the claim: the tests use double precision with relative slack 10⁻⁹ to 10⁻¹³, not directed rounding. This is strong evidence of completeness, not a formal interval proof.

**The record's "14 distinct roots, sign sum 6" is not reproducible.** Running the author's own `roots_trop.py system_deg8_N10_small.json 3` gives the same 10 roots and "sum of Jacobian signs: 8". I do not know which system the 14 / 6 figures came from (the unminimised `system_deg8_N10.json` written by `deg8.py` is not in the bundle; label 6 of the family has 12 roots with sum 6). §6 item 1 of the record should be withdrawn and replaced by the table above.

**1c. Own degree-8 systems (extra).** With my own hill-climber (`climb.cpp`) on the author's 22 realizations (shifted to nonnegative exponents) I obtained an exact degree ±8 system on **all 22** polytopal types, one more than the record's "21 of 22" (the 22nd is author index 11069). The 21 systems the record mentions are not in the bundle, so I could not check those themselves. Complete root search finished for 9 of my systems (types 11069, 79831, 80377, 82118, 85188, 93364, 95720, 97480, 194690): 8 to 14 roots, sign sum equal to the exact degree (±8) in every case, no unresolved box. The other 13 hit my 150-second limit (exponents up to 4,059) and are not counted. Two more matches from a first run (types 48162, 48831) are in the log only; those two systems were not saved.

So Bihan's formula is corroborated, with the sign, on 18 systems of this record.

## 2. Sphere 101272 is not polytopal

**2a.** 30 facets; every triangle in two facets; every vertex link a 2-sphere; f = (10, 40, 60, 30); every vertex in 12 facets; |Aut| = 20, one vertex orbit; the five non-edges (0,1), (2,9), (3,8), (4,6), (5,7) form a perfect matching, and the matching is not an automorphism. It is `manifold_3_10_249002` in Lutz's list, typed S³ there. Author index 128691 is `manifold_3_10_248808`, |Aut| = 240, matching is an automorphism, polytopal: consistent with the record's remark. My automorphism counts for all 247,882 spheres reproduce Lutz's Table 3 line by line (for example 240,683 with trivial group, 6,675 of order 2, three of order 20).

**2b. The reduction is right.**
- A simplicial polytope keeps its type under small perturbation of the vertices, so a polytopal K has a realization in general position.
- For a facet F, all other vertices are on one side, so s_F = sign det[(1,a_F),(1,a_w)] does not depend on w. For two facets sharing a triangle, with the extra vertices at positions i and j, a direct computation gives s_g = −(−1)^{i+j} s_f, which is the coherence rule. K is a connected closed pseudomanifold, so it has exactly two coherent orientations, ±ε. Hence s = ±ε. A reflection of the point set negates the chirotope and s, and preserves realizability and the Grassmann–Plücker relations, so s = ε may be assumed.
- Sufficiency of the 3-term relations is **not needed**. Every chirotope of a realization satisfies them, so the enumeration returns a superset of the chirotopes of realizations; non-realizability of every member of the superset proves the claim. (For uniform chirotopes the 3-term relations are in fact sufficient, so the four are genuine oriented matroids, but nothing rests on this.)

**2c.** Own encoder (`enc.Base(10, 5)`, clauses from a truth table of the identity): 252 variables, 67,200 clauses, 130 variables fixed by the 30 facets, 122 free. Full enumeration with blocking clauses gives **exactly 4** models with CaDiCaL 1.5.3, CaDiCaL 1.9.5, Glucose 4, MiniCard and Lingeling. The four are identical to the author's `chi27.json`.

**2d.** The inequality is right: three nonzero terms with sum 0 and signs not all equal; the minority-sign term is minus the sum of the other two, so its absolute value is strictly larger than each. In logarithms it is linear, and a nonnegative integer combination, not all zero, that sums to the zero vector is a contradiction. For each of the four chirotopes I regenerated the 8,400 inequalities from the chirotope (4,200 relations, each consistent with the signs), then checked for the author's certificate: every row is one of these inequalities; weights are nonnegative integers; the weighted sum is the zero vector in exact integer arithmetic. Rows 202, 205, 211, 210. All four hold.

**2e. The non-polytopality is already known.** Firsching (Math. Program. 166 (2017); arXiv:1508.02531) classifies the 10-vertex simplicial 3-spheres: 162,004 polytopal, with rational coordinates published (now at github.com/mo271/polytope-data, file `polytopalsimplicialspheres10_4.txt.xz`). His numbering differs from the file `3_manifolds_10_all.txt`, so I matched by convex hull and canonical form: the 162,004 realizations hit 162,004 distinct classes of Lutz's list (floating-point hulls; exact rational hulls for the 22 below). **Sphere 101272 is not among them.** So claim B's "not polytopal" is a re-proof of a known fact, by the same kind of certificate Firsching uses. I did not trace the sphere in Altshuler or Bokowski–Sturmfels.

**2f.** Own rank-5 test on Lutz's lists: V = 5, 6, 7 all SAT; V = 8: 37 SAT, 2 UNSAT (`manifold_3_8_19`, `_23`); V = 9: 1,143 SAT, 153 UNSAT. Against Firsching's 1,142 polytopal 9-vertex spheres: the 153 are all non-polytopal, and the one non-polytopal sphere that passes is `manifold_3_9_1121`. At V = 10: 0 UNSAT on 1,500 random polytopal spheres; 1,494 UNSAT on 1,500 random non-polytopal ones.

## 3. The degree-9 vertex map

`vertexmap_deg9_nonpolytopal.json`: entries at most 19 in absolute value; every 4 columns independent; the stored orientation is coherent; degree **9** at 200 random generic u and at u = −C_v for all ten v; ε_G · sign det C_G = +1 for all 30 facets; every tested u lies in exactly 9 facet cones. All as stated.

## 4. Upper bound sweep

**4a. Sphere list.** Regenerated with the author's code (247,882 in 259 s; entry 101272 equals `sphere_deg9.json` facet for facet). Lutz's `3_manifolds_10_all.txt.gz` and `_type.txt` (TU Berlin memorial mirror of the Manifold Page) contain 249,015 complexes, 247,882 typed S³. With my own canonical form (`canon.cpp`, no nauty): **identical sets of isomorphism classes, 0 only-Lutz, 0 only-author.**

**4b. Sweep.**
- The activation-literal construction in `run10.py` is sound. Each sphere's clauses carry ¬act and use fresh auxiliary variables (the offset arithmetic gives no collisions); solving under the assumption act decides GP ∧ sphere clauses; the later unit ¬act satisfies them for good; learned clauses are implied by the clause set, so they stay valid. "TRIV" (a vertex with fewer than d facets avoiding it) is a correct shortcut.
- Own encoder (`enc.py`, rewritten): circuit form of the cone test with one shared pair of variables per 5-subset, own sequential counter, spheres from Lutz's file, own orientation. It first reproduces the V ≤ 9 results (SAT counts 1; 2; 597, 24, 0).
- **All 247,882 spheres, "≥ 8"** (Glucose 4, incremental): UNSAT 247,712, trivial 131, **SAT 39**. Mapped through the canonical forms, this equals the author's results file on every sphere, including the 131.
- Fresh solver per sphere, no incremental state: 32,316 spheres (Glucose 4) and a random sample of 6,000 plus the 39 (CaDiCaL 1.5.3 and Lingeling): no disagreement.
- "≥ 9" on the 39: SAT only for `manifold_3_10_249002` (= 101272). "≥ 10" on the 39: all UNSAT. Five solvers each.
- Not done: the full sweep with a fresh solver per sphere; proof certificates.

**4c. Positive controls** (chirotope fixed to the true signs; SAT at the exact degree, UNSAT one above): 300 random integer matrices on random spheres (degrees 0–2); 200 hill-climbed on random spheres (degrees 2–7, mostly 5 and 6); 78 hill-climbed on the 39 exceptional spheres (degree 6: 4, 7: 22, 8: 51, 9: 1); the nine family systems and the degree-9 map, each with both orientations (−9 to 9). **598 controls, 0 mismatches.**

**4d.** The 22 realizations: exact integer hull has 10 vertices, (U_A) holds, and the type equals the stated one in all 22. All 22 are in Firsching's list. The other 17: 16 are UNSAT in my rank-5 test (three solvers) and absent from Firsching's list; the 17th is sphere 101272.

## 5. Logic of claim A

Sound, given 1–4.
- N = 10, (U_A), (U_C): conv(A) is a simplicial 4-polytope with V ≤ 10 vertices and deg F = ± degree of the vertex map on its boundary (Theorem 1.3).
- V = 10: the boundary is a polytopal sphere. 247,843 spheres have no vertex map of degree ≥ 8; of the 39 others, 17 are not polytopal and the 22 polytopal ones have none of degree ≥ 9.
- V ≤ 9: |deg| ≤ 6 by the reviewed N = 9 runs, which cover V = 5..9 with any (U_C) columns.
- 0..8 are attained with exactly 10 monomials (1a, corroborated by 1b).

What it rests on: Bihan's Theorem 1.3 (and Theorem 3 for D(5,4)), both in unrefereed preprints; SAT answers without certificates, now from two independent encoders and four solver families. Firsching's classification is used only as a cross-check, not in the proof.

## Novelty

- D(4,5), D(5,4): nothing found (alphaXiv, web). Bihan's paper states the maximum is unknown for n, k ≥ 3 and names n = k = 3 as the first open case. The paper is one day old; check again before posting.
- Non-polytopality of sphere 101272: **not new** (Firsching 2017).
- Degree 9 on a non-polytopal sphere against 8 on polytopes: I found no earlier statement. Nearby literature studies simplicial maps of triangulated spheres onto the boundary of a simplex (Ryabichev arXiv:2512.01137; Basak–Gupta–Trivedi arXiv:2409.00907; arXiv:2511.10870). Those are a special, degenerate case of vertex maps and reach much smaller degrees (degree 4 needs 10 vertices in dimension 3). They should be cited as related work.

## Statements stronger than the evidence

1. **§6 item 1** ("14 roots, sign sum 6"): wrong for the bundled system; see 1b.
2. **§3, "it is not polytopal"** reads as a new result. It is a known fact with a new, independent certificate. §6 item 6 says "should already contain it"; it does.
3. **§4, "confirmed by convex hull"**: `realize10.py` uses a floating-point hull. The exact check is done here and holds.
4. **§4, "On 21 of the 22, an exact degree-8 system was built"**: those systems are not in the bundle. Degree 8 holds on all 22 (my systems are in the package).
5. **§7, "From N = 10 on … realizability enters"**: shown for N = 10 only.
6. **§3, validation**: "none of the 1,064 types realized by sampling is UNSAT" was not checkable from the bundle; the stronger statement in 2f holds.
7. The bundle cannot be re-run as shipped: `pool10_all.txt`, `spheres10.json`, `pool10_sat8.txt`, `s8_a.txt`, the `search10` binary and the 21 systems are missing.

## Known gaps (record §6), status after this review

| Gap | Status |
|---|---|
| 1. Root count | Closed: 10 roots, sum +8 |
| 2. No certificates at V = 10 | **Open.** Mitigated: second encoder, other solvers, fresh-solver sample |
| 3. Completeness of the 4 chirotopes | Mitigated: second encoder, five solvers. No certificate |
| 4. Sphere list by count only | Closed: equal to Lutz's list |
| 5. No positive controls at V = 10 | Closed: 598 controls up to degree 9 |
| 6. Novelty of the non-polytopality | Closed: known (Firsching) |

## Required corrections

1. Replace §6 item 1 by the root table of 1b (10 roots, sum +8) and say how completeness was argued.
2. State that sphere 101272 is Lutz's `manifold_3_10_249002`, that its non-polytopality follows from Firsching's classification, and that the four final polynomials are an independent check.
3. Replace "count equals Lutz's" by the class-by-class comparison; name the files.
4. Say "exact convex hull" only where an exact hull was computed; cite Firsching's list for the 22 and the 17.
5. Change "21 of the 22" to all 22, or supply the 21 systems.
6. In §7 write "at N = 10", not "from N = 10 on".
7. Keep both caveats next to the claim: conditional on Bihan's Theorems 1.3 and 3; computer-assisted without certificates.
8. Before publishing: DRAT/LRAT certificates for the 247,712 + 38 + 1 UNSAT answers and for the final UNSAT of the chirotope enumeration.
9. Ship the missing inputs (item 7 above) or a script that rebuilds them.

## Reviewer's files (`topdeg-D45-coldreview-code.tar.gz`)

- `code/exact.py` exact verifier; `code/bnb.cpp`, `mkinput.py`, `rootcheck.py` root search; `code/climb.cpp`, `deg8_own.py` own degree-8 systems.
- `code/sph.py`, `enc.py`, `sweep.py`, `run_small.py`, `controls.py`, `cmp_sweep.py` encoder, sweep, controls.
- `code/canon.cpp`, `cmp_lists.py`, `prep10.py` list comparison; `code/poly5.py`, `item23.py`, `item2cd.py`, `item2f.py`, `item4d.py`, `firs.py`, `firs2.py`, `v9check.py` sphere 101272, rank-5 test, Firsching.
- `out/` results: full sweep, samples, root lists, own systems, logs.
- External data (not packaged): Lutz `3_manifolds_10_all.txt.gz`, `3_manifolds_10_all_type.txt`, `3_manifolds_all.txt`, `3_manifolds_all_type.txt` (www3.math.tu-berlin.de/IfM/Nachrufe/Frank_Lutz/stellar/); Firsching `polytopalsimplicialspheres10_4.txt.xz`, `polytopalsimplicialspheres9_4.txt` (github.com/mo271/polytope-data).
