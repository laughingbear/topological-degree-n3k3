# REVIEW of `topdeg-review-2026-11-27` (cold referee report, 2026-11-27)

## 0. Verdicts

| claim | verdict | basis (all referee code, exact integer arithmetic unless stated) |
|---|---|---|
| A1 (V=12, deg 13, q=0) | **CONFIRMED** | shelling ⇒ PL 3-sphere; degree 13 at 64 generic u; q=0; Thm M 52 < 70.5 |
| A2 (V=13, deg 16, q=5) | **CONFIRMED** | shelling; degree −16 (so \|deg\|=16); q=5; Thm M 64 < 82.5 |
| A3 (V=14, deg 17, q=0) | **CONFIRMED** | shelling; degree 17; q=0; Thm M 68 < 87.5 |
| A4 (V=16, deg 25, q=7) | **CONFIRMED** | shelling; degree 25; q=7; Thm M 100 < 138.5 |
| A5 (V=18, deg 35, q=0) | **CONFIRMED** | shelling; degree 35; q=0; Thm M 140 < 166.5 |
| A6 (V=20, deg 41, q=8) | **CONFIRMED** | shelling; degree 41; q=8; Thm M 164 < 233.5 |
| B (none of the six is polytopal) | **CONFIRMED** | own encoder, UNSAT by unit propagation for all six; machine-checked proofs of 6–60 steps; hand-checkable 7-line contradiction for V=12 (§2.3) |
| C (12-vertex 4-polytope, deg 11, q=6) | **CONFIRMED** | exact hull = the 48 listed facets, f=(12,60,96,48); degree 11; q=6 |
| D, combinatorial statement ((U_A), (U_C), conv(A), formula value 11) | **CONFIRMED** | all 792+495 minors nonzero; hull = same 48 facets; formula = 11 at 305 generic u |
| D, "deg F = 11" as a statement about the map F | **CONFIRMED conditional on Bihan's Theorem A / 1.3** (unrefereed preprints, 25 Nov 2026); independent support: 13 positive zeros certified in exact rational interval arithmetic, Jacobian signs 12×(+1), 1×(−1), signed count **+11**; a float64 branch-and-bound over all of ℝ⁴_{>0} finds exactly these 13 and no others (§4.3 — numerical evidence, **not a proof**). |

I found **no error that changes a claimed number**. I found a number of small errors and
weaknesses in the write-up and scripts (§6), and one framing problem about what the results are
being compared with (§5, E3).

Independence statement. Every confirming check above was computed by code in `referee/` that I
wrote from `CLAIMS.md`. I read `ours/` and ran `run_ours.py`, `verify_exact.py`, `polysat2.py`
once, as a comparison only (they agree with me everywhere). No function from `ours/` is imported.
The stored `orientation` arrays are not used; I recompute the orientation and then *compare*: mine
equals yours on all seven complexes (same normalisation, +1 on the first facet).
Two things are not independent of outside sources: Claim D's passage from "formula = 11" to
"deg F = 11" uses Bihan's theorem, and the sphere certificate uses the textbook theorem that a
shellable closed pseudomanifold is a PL sphere.

---

## 1. Claims A1–A6

### 1.1 What I ran
- `sphere_cert.py` — (a) closed 3-manifold test: every triangle in exactly two facets, every
  vertex link a connected closed surface with χ = 2, every edge link a single cycle, dual graph
  connected, χ(K) = 0. (b) **3-sphere by shelling**, not flips: a randomised greedy search finds
  a shelling order, and a separate function re-verifies it from the definition (for each k ≥ 2,
  F_k ∩ (F_1 ∪ … ∪ F_{k−1}) is a non-empty union of triangles of F_k). A shellable simplicial
  complex in which every triangle lies in exactly two facets is a PL sphere (Danaraj–Klee;
  Björner; Ziegler, *Lectures on Polytopes*, Ch. 8). Shellings are stored in `out/shelling_*.json`
  and can be checked with ten lines of code. (c) Integer homology by Smith normal form:
  H = (ℤ, 0, 0, ℤ), no torsion, for all seven complexes — over ℤ, not only mod small primes.
- `check_A.py` — orientation recomputed from the definition and verified on every triangle;
  all facet determinants nonzero; degree by Cramer's rule in Python integers at 4 fixed + 60
  random u with entries up to 10⁹ (none was non-generic); q; Theorem M as 8|deg| < 4F − E − q.
- `controls_A.py` — addresses "right for the wrong reason": (1) five random integer vertex maps
  per sphere give a u-independent signed count (values in {−3,…,1}); (2) flipping the orientation
  sign of a single facet destroys u-independence of the claimed maps for 50/51, 57/61, 62/63,
  94/100, 114/117, 152/165 facets (40 test vectors u). So constancy in u is a real test of the orientation, and it passes.
  (3) In `check_CD.py` the radial map of the Claim-C polytope onto itself has all signs +1 and
  degree +1 with the same code.

### 1.2 Results

| | f (recomputed) | deg (o=+1 on first facet) | q | facets whose cone contains u (min–max over u) | Thm M |
|---|---|---|---|---|---|
| A1 | (12,63,102,51) | +13 | 0 | 13–13 | 52 < 70.5 |
| A2 | (13,74,122,61) | −16 | 5 | 16–20 | 64 < 82.5 |
| A3 | (14,77,126,63) | +17 | 0 | 17–17 | 68 < 87.5 |
| A4 | (16,116,200,100) | +25 | 7 | 25–29 | 100 < 138.5 |
| A5 | (18,135,234,117) | +35 | 0 | 35–35 | 140 < 166.5 |
| A6 | (20,185,330,165) | +41 | 8 | 41–43 | 164 < 233.5 |

All agree with the table in `CLAIMS.md`. For the q = 0 cases every generic u is covered exactly
|deg| times, as it must be.

### 1.3 Remarks
- Error class (a) "manifold but not a sphere" is excluded by the shellings.
- Error class (b): excluded as far as I can test it (§1.1 controls, plus Claim D where the
  orientation comes from geometry rather than from the combinatorial rule and gives the same 11).
- The degree is only defined because no facet is degenerate; I assert this for every facet.

---

## 2. Claim B (non-polytopality)

### 2.1 The reasoning
- *General position is legitimate.* For a simplicial polytope every facet is a simplex whose
  supporting hyperplane has all other vertices strictly on one side; this is an open condition, so
  a small perturbation of the vertices keeps the same boundary complex, and a generic perturbation
  makes all 5×5 minors of (1,p_i) nonzero. So K polytopal ⇒ K = ∂P with a uniform chirotope.
- *One global sign is enough.* For a facet t the sign χ(t,w) is the same for all w ∉ t (all on
  the inner side), and t ↦ χ(t,w) is a coherent orientation of ∂P in the sense of `CLAIMS.md`
  (standard; I also verified it numerically, see below). A strongly connected pseudomanifold has
  exactly two coherent orientations, ±o. Hence χ(t,w) = ε·o_t with one ε. Fixing ε = +1 loses
  nothing because −χ is the chirotope of the mirror image. I ran both ε anyway.
- *Soundness of the test.* Only necessity of the constraints is used: every real point
  configuration satisfies the three-term Grassmann–Plücker identities, so "the three signed
  products are not all equal" is necessary. UNSAT therefore proves non-polytopality. (It in fact
  proves a little more: no uniform rank-5 oriented matroid, realisable or not, has K as its
  boundary complex in this sense.) SAT would decide nothing, as `polysat2.py` correctly prints.
- *The sign convention is the right one.* This is where a wrong-for-the-wrong-reason UNSAT could
  hide, so I tested it on real polytopes (`chiro_controls.py`): the Claim-C polytope, conv(A) of
  Claim D, the cyclic polytope C(12,4) on the moment curve, and a random 11-vertex polytope. For
  each, the actual chirotope computed from coordinates satisfies every unit condition with a
  single ε, violates 0 of the 27,720 (resp. 11,550) encoded GP relations, and the propagator
  reports no conflict and derives only true signs (690, 690, 792, 171 derived; 0 wrong).

### 2.2 What I ran
- `chiro.py` — my encoder and a queue-based propagator (rule: if two terms of a GP relation are
  known and equal, the third is forced). Conflict for all six spheres, both ε, in under a second.
  It writes a *minimised* proof (the backward cone of the conflict) to `out/proof_*.json`.
- `check_proof.py` — a stand-alone checker that shares no code with `chiro.py` (own permutation
  sign, own orientation, own GP expansion). All 12 proofs: VALID.
- `sat_cross.py` — CaDiCaL on the CNF from my encoder: V=12, 13, 14 UNSAT for both ε (443,928 /
  961,509 / 1,922,550 clauses); the Claim-C polytope SAT. (V ≥ 16 not run through CaDiCaL; the
  checked propagation proofs make it unnecessary.)

| sphere | proof size (facet conditions + GP steps) | final violated relation (s; q) |
|---|---|---|
| V=12 | 6 + 0 | (0,1,6; 2,7,9,11) |
| V=13 | 41 + 19 | (1,2,7; 0,6,9,11) |
| V=14 | 12 + 4 | (0,2,5; 1,3,4,9) |
| V=16 | 6 + 0 | (0,1,9; 2,8,10,13) |
| V=18 | 17 + 5 | (0,2,6; 1,3,16,17) |
| V=20 | 6 + 0 | (0,1,16; 2,11,15,17) |

"Found by unit propagation alone" is confirmed. For V = 12, 16, 20 it is even simpler: a single
GP relation is violated by the facet conditions themselves.

### 2.3 Human-readable contradiction for V = 12
Notation: χ[abcde] with a<b<c<d<e. Normalise o(0,1,2,6) = +1, ε = +1.

Orientations of the four facets used (each arrow is one shared triangle; rule o_t′ = −(−1)^{i+j} o_t,
i, j the 0-based positions of the deleted vertices):

    [0,1,2,6](+) → [0,1,2,7](−) → [0,1,7,11](−)
    [0,1,2,6](+) → [0,2,6,8](−) → [0,6,8,11](+) → [0,6,9,11](−)
    [0,1,2,6](+) → [0,1,2,7](−) → [1,2,3,7](+) → [1,3,7,9](−) → [1,6,7,9](+)

Facet conditions χ(t,w) = o_t, then sorting the indices:

1. facet 0126, w = 7, 9, 11:  χ[01267] = χ[01269] = χ[0126 11] = +1.
2. facet 0 1 7 11, w = 6:  χ(0,1,7,11,6) = −1; moving 6 two places is even ⇒ χ[0 1 6 7 11] = −1.
3. facet 0 6 9 11, w = 1:  χ(0,6,9,11,1) = −1; moving 1 three places is odd ⇒ χ[0 1 6 9 11] = +1.
4. facet 1 6 7 9, w = 0:  χ(1,6,7,9,0) = +1; moving 0 four places is even ⇒ χ[0 1 6 7 9] = +1.

Grassmann–Plücker with s = (0,1,6), (a,b,c,d) = (2,7,9,11):
[s,2,7][s,9,11] − [s,2,9][s,7,11] + [s,2,11][s,7,9] = 0. Here [0,1,6,2,x] = −χ[0,1,2,6,x] = −1
for x = 7, 9, 11 (one transposition). So the three terms are

    (−1)(+1) = −1,    −(−1)(−1) = −1,    (−1)(+1) = −1,

three strictly negative numbers summing to zero. Contradiction.

Geometric meaning: {0,1,6} is a triangle of K. Projecting a hypothetical polytope along the plane
of that triangle gives a convex polygon in which the triangle is a vertex and the other vertices
lie in a wedge of angle < π, so "angular order seen from the triangle" is a linear order. The
facets 0·1·7·11, 0·6·9·11 and 1·6·7·9 force 7 ≺ 9, 9 ≺ 11 and 11 ≺ 7. The V = 16 and V = 20
proofs have the same shape (triangles {0,1,9} and {0,1,16}; `b_human.py` prints them).

---

## 3. Claim C

`check_CD.py`, exact integers. For each of the 495 four-subsets of `polytope_points` I compute the
signs of the 5×5 determinants with every other point. Exactly 48 subsets have all other points
strictly on one side; none has a point on its hyperplane; no 5×5 minor of (1,p) vanishes at all
(the points are in general position). The 48 facets are exactly the listed ones, every point is a
vertex, f = (12,60,96,48), the complex is a closed 3-manifold and shellable. The side signs equal
ε·o_t with one ε. Degree of `C_int`: 11 at 80 generic u (covering numbers 11–13), q = 6.
Theorem M: 44 < 63. **CONFIRMED.**

The baseline "published construction gives 10 for 12 vertices" I could only trace to your own
same-day announcement (§5, E3); I did not verify it.

---

## 4. Claim D

### 4.1 Combinatorial statement — CONFIRMED
- A ⊂ [0,95]⁴, 12 distinct points; C ∈ [−14,14] (two entries are 0, which is harmless).
- (U_A): all 792 minors nonzero, smallest |minor| = 7,722. (U_C): all 495 minors nonzero,
  smallest |minor| = 15.
- conv(A): 48 simplicial facets, no degenerate support, all 12 points are vertices, and the facet
  list is *identical as labelled sets* to the facet list of Claim C.
- C is not a rescaling of Claim C's `C_int` column by column, but sign det C_Γ agrees with sign
  det C_int,Γ on all 48 facets up to one global sign, which is all the formula sees.
- Formula of Theorem 1.3 with ε_Γ = sign det[(1,a_s)_{s∈Γ},(1,a_i)] (i off Γ, last): value **11**
  at all 305 generic u tried (5 fixed, 300 random; number of contributing facets 11, 13 or 15).
  Genericity of u is tested exactly (all 220 three-column determinants with u nonzero).

### 4.2 The sign convention (your worry (d))
- I read the statement in both preprints (arXiv:2611.34336 Thm 1.3; arXiv:2610.29808v2 Thm A).
  Your transcription is accurate: (−1)ⁿ, columns of Γ increasing, i_Γ last, Â with first row 1.
- ε_Γ is a coherent orientation: on conv(A) it equals ε·o_Γ for my recomputed combinatorial
  orientation. Hence the formula is ± the degree of the vertex map a_s ↦ C_s on ∂conv(A) in the
  sense of `CLAIMS.md` §0, whatever the conventions; |value| = 11 cannot be an artefact of
  (−1)ⁿ, of the row/column layout of Â, or of where i_Γ is placed. Each of those can only flip
  the global sign.
- n = 1 sanity check done by hand: the formula gives [c_N>0] − [c_1>0], the correct degree.
- So the only non-elementary input is the theorem itself: deg F equals the degree of that
  piecewise-linear boundary map. I read the proof sketch in 2610.29808v2 (Lemmas 3.1, 3.2 and
  Prop. 1.3): moment-map compactification, straight-line homotopy on ∂P kept away from 0 by the
  half-space condition. I see no gap, but both preprints are two days old and unrefereed, and the
  author states they were written with AI assistance. That is why the second D row in the table
  is "conditional".

### 4.3 Independent count of positive zeros
Feasible in part.

*Discovery (floating point / multiprecision, `D_seek.py`, `D_seek_fast.py`, `D_refine.py`).*
Newton in t = log x from 104 "Viro" seeds (5-subsets of columns with a positive dependence),
1,800 random seeds at 80 digits and 640,000 random seeds in float64. Three of the zeros were hit
once or not at all in 640,000 random starts, so random search alone is not a completeness
argument for this system.

*Certification (exact, `D_certify.py`).* For each candidate a Krawczyk test in exact rational
interval arithmetic on a box of relative radius 10⁻³⁰ around a 200-bit rational point: unique
zero in the box, DF nonsingular on the box, and the sign of det DF from an interval determinant
that excludes 0. Monomial ranges over a box are exact because monomials are monotone on the
positive orthant. All 13 pass:

| x₁ | x₂ | x₃ | x₄ | sign det DF |
|---|---|---|---|---|
| 0.10206977 | 5.10817424 | 1.49708126 | 2.36572654 | + |
| 0.181126471 | 1.1168416 | 0.471066645 | 0.130128669 | + |
| 0.54656299 | 12.2030684 | 3.63806416 | 0.827519061 | + |
| 0.972025443 | 1.04928521 | 0.973463805 | 0.905281111 | + |
| 1.01767726 | 0.976537296 | 0.981653475 | 1.0516978 | **−** |
| 1.01869509 | 0.972962608 | 0.982796622 | 1.04077908 | + |
| 1.03717257 | 0.814547441 | 0.994252518 | 1.15096644 | + |
| 1.14724899 | 1.15835236 | 1.31989595 | 0.792418191 | + |
| 1.19491786 | 1.42674794 | 0.832895498 | 0.58659195 | + |
| 1.34397875 | 0.787602949 | 0.700462656 | 2.93605725 | + |
| 2.19749439 | 12.8670673 | 0.119746632 | 0.00600101998 | + |
| 3.3048786 | 0.00671951098 | 1.92003169 | 11.2812947 | + |
| 7.68032411 | 0.297335904 | 0.0481395315 | 37.2869888 | + |

**Proved:** F has at least these 13 nondegenerate positive zeros; their signed count is +11.
(60-digit values of log x are in `out/D_roots_all.json`.) This matches the formula, including
the sign with the convention det DF in the variables x.

*Completeness (float64, `D_exclude.py`) — evidence, not proof.* The search covers all of
t-space ℝ⁴ by 73 jobs: the cube max|t_k| ≤ 1/16, 64 conic shells 1/16 ≤ max|t_k| ≤ 512, and 8
far-field cells max|t_k| ≥ 512. Every cell is either excluded, shown to contain exactly one zero
(Krawczyk), or subdivided, until nothing is left. Exclusion test: after dividing by the locally
largest monomial, each monomial ratio lies in an interval over the cell; the cell has no zero if
some w ∈ ℝ⁴ makes Σ_i min((w·C_i)l_i, (w·C_i)h_i) > 0; w comes from a linear program (HiGHS) and
is then re-evaluated directly with a relative margin 10⁻⁷. In the far field the bound uses only
that the non-dominant ratios are ≤ exp(512·gap) and that the dominant set is a face.
Result: all 73 jobs completed, **0 unresolved cells**, 2,054,859 LPs, about 308,000 subdivisions;
the cells certified to contain a zero give **exactly the 13 zeros above and no others**
(2 in the central cube, of opposite sign; 11 in the shells; none beyond max|t_k| = 5.12).
A built-in self-test evaluated the excluding functional at 5,108,997 random points of excluded
cells: 0 failures. Logs: `out/excl_central.log`, `out/excl_even.log`, `out/excl_odd_a.log`,
`out/excl_odd_b.log`; summary by `D_union.py`.

Limits of this evidence: float64 throughout, an LP solver in the loop, and about 300 lines of my
own bounding code that nobody has reviewed. It is a careful numerical computation, not a proof,
and by your rule 2 it does not count as a confirming check. Making it rigorous would need the
same subdivision with outward-rounded interval arithmetic and exact verification of each w
(about 2·10⁶ cells; feasible in C, not attempted).

What this does and does not establish. Proved: 13 nondegenerate positive zeros with signed count
+11. Numerically supported but not proved: there are no others. If the list is complete then
deg F = +11 independently of Bihan's theorem (the zeros are nondegenerate and the zero set is
compact), so the direct count and the formula agree. I regard "deg F = 11" as correct; the
formal status is: theorem-dependent, plus a non-rigorous but exhaustive-looking direct count.

### 4.4 The consequence
"The maximal topological degree for n = 4, N = 12 (k = 7) is at least 11" follows from 4.1 +
Bihan's theorem. Bihan's own bound (Thm 2 of 2611.34336) gives |deg F| < ½·M(4,12) = 27, and
your Theorem M with a neighbourly f-vector gives < 18.75, so 11 is consistent with both.

---

## 5. Comments on E1–E3

### E1 (enumeration)
- The symmetry statements in the certificate notes are true (`check_sym.py`): V=12 and V=14 are
  rotation-invariant (and in fact also reflection-invariant, c = 5 resp. 4, which the notes do not
  say); V=16 is D₈-invariant with c = 1; V=18 is D₉-invariant with c = 0; V=20 is D₁₀-invariant
  with c = 1; the Claim-C complex is Z₆- but not reflection-invariant.
- Logic of `enum.c`: sound for what it does. Root loop + permanent exclusion of earlier roots
  gives each solution once (rooted at its smallest orbit); branching on a triangle orbit with
  count 1 gives a unique path to each solution, because exactly one of the candidate orbits can
  belong to a given solution; the sibling exclusion `state = −1` is harmless for the same reason.
- I reproduced the n = 5 count with an independent enumerator using a different search
  (include/exclude over orbits, `e1_enum.py`): there are 3,012 non-empty invariant unions with
  every triangle in 0 or 2 facets; 2,276 are strongly connected; **2,306** have components
  forming a single group orbit. Your "closed 2306" is the third number. For the dihedral problem
  (n = 5, c = 0) the three numbers are 292 / 164 / 170 and your binary prints 170. So the code is
  right, but the label is loose: 30 (resp. 6) of the counted objects are not pseudomanifolds in the
  usual (strongly connected) sense. No sphere is lost by this. I did not reproduce 1,503,130
  (n = 6); my Python enumerator is too slow for it.
- The link filter in `report()` (connected, 2·nv − nf = 4) is necessary but not sufficient for a
  2-sphere link: two spheres glued at two points pass it. As a pre-filter that is fine; it must
  not be described as a sphere test. Checking only a_0 and b_0 is correct by transitivity.
- Unchecked buffer: `tmp[64]` is written only for c < 64, but `nopt = c` and
  `memcpy(opts, tmp, sizeof(int)*c)` run with the full c; `optbuf[64][64]`, `chosen[64]`,
  `par[64]`, `used[64]` likewise assume ≤ 64. For V ≤ 20 a triangle orbit has at most 17
  candidate orbits and a solution at most ~33 orbits, so nothing overflows here; add asserts.
- Coverage question I cannot answer from the tarball: for even n the reflections with c and c+1
  are not conjugate by a shift (shifting block b by k changes c by 2k), so both c = 0 and c = 1
  must be run for n = 10 (V = 20), while one c suffices for n = 9. Please state which were run.
  Also the class fixes the rotation to act with the same step on both blocks; a_i ↦ a_{i+1},
  b_i ↦ b_{i+m} with m ≠ ±1 is a different class and is not covered. "Enumerated completely"
  should be read with those two qualifications.

### E2
Nothing to add; no maximality is claimed and none follows. Theorem M leaves room (e.g. it allows
|deg| ≤ 17 for the f-vector of A1).

### E3 (prior art) — what I found
Searches: web search, arXiv listing pages, alphaXiv index, the four references you gave.
Not exhaustive (no MathSciNet/zbMATH access).

1. A. Ryabichev, *On vertex-minimal simplicial maps to the sphere*, arXiv:2512.01137v3
   (24 Sep 2026), https://arxiv.org/abs/2512.01137 — λ(n,d)^h is linear in d; all constructions
   polytopal; the Appendix studies the periodic colouring φ_k of ∂C(k, n+1) and asks (Question 3)
   when it is degree-maximal; Prop. 4: a degree-d simplicial map K → ∂Δ^{n+1} needs F(K) ≥ d(n+2).
2. K. Apolonskaya, O. Musin, *Minimal simplicial spherical mappings with a given degree*,
   arXiv:2511.10870 (Math. Notes 119 (2026)), https://arxiv.org/abs/2511.10870 — λ(n,d) for small d.
3. B. Basak, R. K. Gupta, A. Trivedi, *Simplicial degree d self-maps on n-spheres*,
   arXiv:2409.00907v2, https://arxiv.org/abs/2409.00907 — constructions, facet-minimal maps.
4. Q. Wang, *Exact Periodic Degrees, Degree Spectra, and Sharp Thresholds for Cyclic-Polytope
   Colorings*, https://zenodo.org/records/22249508, and *Vertex–Facet Extremality of Simplicial
   Sphere Maps*, https://zenodo.org/records/23430458 — deg φ_k = q^{h−j}(q+1)^j; nine-vertex
   3-sphere census: maximum |degree| 3.
5. F. Bihan, *Topological degree of sparse polynomial maps*, arXiv:2610.29808v2,
   https://arxiv.org/abs/2610.29808, and *Descartes' rule of signs for the topological degree*,
   arXiv:2611.34336, https://arxiv.org/abs/2611.34336 — the latter says explicitly that the
   maximal topological degree is not known for n, k ≥ 3 and gives only the upper bound
   ½·min(M(n,N), M(k,N)); it contains no lower-bound construction for n = 4.
6. Your own public repository, https://github.com/laughingbear/topological-degree-n3k3 — banner
   dated 2026-11-27: (N−1)(N−3)/8 ≤ max|deg F| < (3N²−11N)/16 for odd N, n = 4; banner dated
   2026-11-26: D(4,5) = {0,…,8}, and maximum 9 over 10-vertex 3-spheres, attained only on a
   non-polytopal sphere.

Assessment.
- I found **no prior occurrence** of vertex maps of degree 13/16/17/25/35/41 on 12–20-vertex
  3-spheres, nor of a 12-monomial system in 4 variables with topological degree ≥ 11.
- **Framing problem.** Items 1–4 are about *simplicial maps to the 5-vertex sphere*
  (colourings). Your objects are not of that kind and cannot be compared with them: for a
  colouring, d·5 ≤ F, so with 12 vertices d ≤ ⌊54/5⌋ = 10, and A1 (13·5 = 65 > 51 = F) is
  impossible as a colouring. For the same reason Wang's formula gives only 4 for the periodic
  colouring of C(12,4) (k = 12 = 2·5 + 2: q = 2, j = 0, deg = q² = 4), not your "cyclic value" 10.
  So A1–A6 say nothing about λ(3,d), and "above the cyclic value" must not be read as "beats
  Ryabichev/Wang". The only literature your numbers can be compared with is Bihan's D(n,k), where
  the sole prior lower bounds for n = 4 are your own, published the same day. Say so explicitly;
  a reader who knows the λ(n,d) papers will otherwise misread the claim.
- The "cyclic value" for even V, (V−2)(V−4)/8, is not in the public statement I could read (the
  banner gives odd N only). It is implied by the odd case and monotonicity in N (stack a vertex on
  a facet and put its vector inside that facet's cone), e.g. N = 12 ≥ N = 11 value = 10. State
  that derivation, or the source.
- Your own N = 10 result (8, exact) already exceeds the "cyclic value" (10−2)(10−4)/8 = 6. So
  "the best construction published so far" is not the best published value at N = 10, and a
  reader should not take the cyclic value as the record in general. For N = 12 the best previously
  published lower bound I can derive from your public statements is max(10, 8) = 10, so Claims C/D
  do improve it by one.

---

## 6. Errors and weaknesses found (all minor; none changes a verdict)

1. `CLAIMS.md` §1 header says "above the cyclic value (except A3)", but A3's 17 is above its
   cyclic value 15. The exception refers to the bracketed odd-V formula evaluated at even V
   (17.875), which is not the value of any construction (it is not even an integer for V = 12,
   14, 16, 18, 20). Either drop the bracketed numbers or explain what they are a bound for.
2. The field `q` is absent from `V12_deg13.json`, `V14_deg17.json`, `V12_polytope_deg11.json`
   (the notes say "all facets same sign" for the first two, nothing for the third, whose q = 6 is
   claimed in the text). `V12_deg13.json` has an undocumented field `"h": 60`.
3. `degree` in the certificates is |deg|; with the stored orientation A2 has degree −16. Harmless,
   but say "absolute value".
4. `ours/README.txt` is garbled: line 8 ("store "C_int" and "degree" instead.") is a fragment and
   the `run_ours.py` line is out of order.
5. `polysat2.py` reports `clauses=443520`, which omits the unit clauses (my count including them
   is 443,928 for A1). Cosmetic. More substantively, your UNSAT evidence was a bare solver answer
   with no proof log; the proofs in `out/proof_*.json` close that gap.
6. `enum.c`: see E1 — mislabelled count, unchecked fixed-size buffers, link filter weaker than its
   comment suggests.
7. Claim B is stated for "the six spheres"; the argument actually shows more (no uniform oriented
   matroid polytope), and for V = 12, 16, 20 it reduces to one sentence. Worth saying, since a
   6-step contradiction is far more convincing than "CaDiCaL says UNSAT".
8. Claim D's consequence rests on a two-day-old preprint. The package should include an
   independent check of the degree; the 13 certified zeros in `out/D_roots_all.json` and the
   search in `D_exclude.py` are a start (your repository apparently has a `roots_krawczyk.py` for
   other systems; nothing of the kind is in this package).
9. `CLAIMS.md` §4 says conv(A) is "the same complex as in Claim C". True, and with the identical
   labelling — but the coefficient matrix is *not* `C_int` (not even up to column scaling); it has
   the same facet signs. If the intent was "same oriented-matroid data on the facets", say that.
10. Novelty framing: see E3.

## 7. What I did not do
- No CaDiCaL/DRAT run for V ≥ 16 (superseded by the checked propagation proofs).
- No *proof* of completeness of the zero set in Claim D; the completeness search of §4.3 is float64.
- No reproduction of the n = 6 enumeration count, and no audit of `post.py`, `climb.c`, `sph.py`,
  `drive.py`, `realize.py` beyond reading; they do not enter any verdict.
- No independent verification of Theorem M or of the "cyclic value" construction.
