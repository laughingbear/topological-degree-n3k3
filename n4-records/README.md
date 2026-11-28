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
| 14 | 17 | 0 | (14,77,126,63) | 15 |
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
