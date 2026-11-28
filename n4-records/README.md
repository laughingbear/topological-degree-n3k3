# n = 4: a better lower bound at N = 12, and sphere records that are not polytopes

Setting: F. Bihan, *Descartes' rule of signs for the topological degree* (arXiv:2611.34336) and
*Topological degree of sparse polynomial maps* (arXiv:2610.29808). For a system of n = 4 real
polynomials in 4 variables with N monomials, the topological degree of the map on the positive
orthant equals, by Bihan's Theorem 1.3, the degree of a piecewise-linear "vertex map" on the
boundary of the Newton polytope. This folder adds two things to the bounds in `n4-bounds/`.

## 1. N = 12: the maximal degree is at least 11

`certs/system_N12_deg11.json` is an explicit system with 12 monomials (exponents in [0,95],
integer coefficients in [-14,14]). Conditions (U_A) and (U_C) hold, the Newton polytope is
simplicial with 48 facets, f = (12,60,96,48), and Bihan's formula gives degree 11.

- Previous lower bound at N = 12 (from `n4-bounds/`): 10, the even-N value (N-2)(N-4)/8 of the
  cyclic construction (explicit system in `n4-bounds/data/cyclic_systems.json`).
- Upper bound at N = 12 (Theorem M in `n4-bounds/`): 18.
- Independent of Bihan's formula: 13 positive zeros of the system are certified in exact rational
  interval arithmetic (Krawczyk), 12 with Jacobian sign +1 and one with -1, signed count +11
  (`referee/D_certify.py`, `referee/out/D_roots_all.json`). A floating-point branch-and-bound over
  the whole positive orthant finds exactly these 13 zeros. That search is numerical evidence, not
  a proof, so "degree 11" formally rests on Bihan's theorem.

`certs/V12_polytope_deg11.json` has the polytope (integer points) and the vertex map.

## 2. Vertex maps on 3-spheres with larger degree; none of these spheres is a polytope

A vertex map assigns a vector in R^4 to each vertex of a triangulated 3-sphere; its degree is
defined as in `n4-bounds/`. Polynomial systems correspond to spheres that are boundaries of
polytopes, so these records bound the sphere version of the question, not Bihan's D(4,k).

| vertices | degree | q (facets of the minority sign) | f-vector | best polytopal value known |
|---|---|---|---|---|
| 12 | 13 | 0 | (12,63,102,51) | 11 (section 1) |
| 13 | 16 | 5 | (13,74,122,61) | 15 |
| 14 | 17 | 0 | (14,77,126,63) | 16 (section 3; was 15) |
| 16 | 25 | 7 | (16,116,200,100) | 21 |
| 18 | 35 | 0 | (18,135,234,117) | 28 |
| 20 | 41 | 8 | (20,185,330,165) | 36 |

Last column: section 1 for 12 vertices; otherwise the cyclic construction of `n4-bounds/`,
(V-1)(V-3)/8 for odd V and (V-2)(V-4)/8 for even V. No maximality is claimed for any entry.

Each certificate in `certs/` lists the facets, the orientation, and an integer vertex map. Each
complex is a 3-sphere (an explicit shelling is in `referee/out/shelling_*.json`), and the degree
is checked in exact integer arithmetic.

**None of the six spheres is the boundary of a 4-polytope.** `referee/out/proof_*.json` are
short proofs (6 to 60 steps) from the three-term Grassmann-Pluecker relations, checked by the
stand-alone `referee/check_proof.py`. For 12 vertices the proof is seven lines by hand: around
the triangle {0,1,6}, the facets {0,1,7,11}, {0,6,9,11} and {1,6,7,9} force the cyclic order
7 before 9, 9 before 11 and 11 before 7 (see `referee/REVIEW.md`, section 2.3). The same already
happens at 10 vertices (`n4k5/`: 9 on a sphere, 8 on polytopes).

How they were found: the spheres at 12, 14, 16, 18, 20 vertices are invariant under a dihedral
group acting on two rings of V/2 vertices, with both rings rotating at the same speed. That class
was enumerated completely for V <= 20 and the vertex map optimised numerically. The 13-vertex
example came from a random search and has no symmetry.

## 3. What this is not

The papers on simplicial maps of spheres with few vertices (Ryabichev arXiv:2512.01137,
Apolonskaya-Musin arXiv:2511.10870, Basak-Gupta-Trivedi arXiv:2409.00907, and Q. Wang's Zenodo
preprints) study colourings, that is, maps to the 5-vertex sphere. A colouring of degree d needs
at least 5d facets, so with 12 vertices d <= 10. The vertex maps here are a wider class and the
numbers are not comparable with those results. The only comparable earlier bounds are the ones
in this repository.

## 4. Verify

```
cd referee && ./run_all.sh ../certs     # about 10 minutes, Python 3.10+; the last step (zeros) needs mpmath
python3 check_proof.py out/proof_*.json
```

`referee/` is the code of an independent cold review (written from the claims, without our
code); `referee/REVIEW.md` is its report: all claims confirmed.

Found with Claude (Anthropic). Laughing Bear, 2026-11-27.

## 3. N = 14: the maximal degree is at least 16 (added 2026-11-27, evening)

`n14/certs/system_N14_deg16.json` is an explicit system with 14 monomials (exponents in [0,117],
integer coefficients in [-20,20]). Conditions (U_A) and (U_C) hold, the Newton polytope is
simplicial with 14 vertices and 56 facets, f = (14,70,112,56), all 56 facets have the same sign,
and Bihan's formula gives degree -16.

- Previous lower bound at N = 14: 15. That is our own bound from `n4-bounds/` (the even-N value
  (N-2)(N-4)/8 of the cyclic construction, explicit system in `n4-bounds/cyclic_systems.json`).
  We know of no lower bound for D(4,9) by anyone else. For comparison,
  `n14/certs/C14_cyclic_deg15.json` is a vertex map of degree 15 on the cyclic polytope C(14,4)
  itself, built the same way.
- Upper bound at N = 14 (Theorem M in `n4-bounds/`): 27. On 14-vertex 3-spheres the best value
  known is 17 (section 2), on a sphere that is not a polytope.
- Independent of Bihan's theorem: the referee counted the positive solutions. The system has
  exactly 30, all nondegenerate, 7 with positive and 23 with negative Jacobian determinant, signed
  count -16. Existence, uniqueness and the signs are proved in exact rational interval arithmetic,
  and so is the absence of solutions outside the cube |log x| <= 8; the absence of a 31st solution
  inside the cube rests on floating-point interval arithmetic with outward rounding
  (`n14/referee/S_outer.py`, `S_inner.cpp`, `S_certify.py`).
- How it was found: the 14 exponent vectors are two orbits of the cyclic group Z_7 acting on R^4 by
  two rotations (angles 2 pi/7 and 4 pi/7), the coefficient vectors are two orbits under rotations
  by 2 pi/7 and 6 pi/7, and everything is then rounded to integers (`n14/search/`). No maximality
  is claimed: 16 is the best value found in this one symmetric family.
- `n14/certs/system_N14_deg16_small.json` (from the referee) has exponents in [0,51] and
  coefficients in [-17,16], the same 56 facets and formula value -16; its solutions were not counted.
- Three coefficients of the main system are zero, so two of the four polynomials have 13 and 12
  monomials; the system has 14 and (U_C) holds.

Review: `n14/referee/REVIEW.md` (independent exact code, all claims confirmed; the referee could
not locate the benchmark 15 because the review package did not cite `n4-bounds/`; it is cited above).
