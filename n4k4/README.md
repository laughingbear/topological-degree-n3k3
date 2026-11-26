# n = k = 4: the maximal topological degree is 6

**Result (2026-11-26).** D(4,4) = {0, 1, 2, 3, 4, 5, 6}.

That is: for fewnomial maps F : ℝ⁴₊ → ℝ⁴ with N = 9 monomials whose exponent and coefficient matrices have nonzero maximal minors, the possible values of |deg F| are exactly 0, …, 6. Theorem 2 of F. Bihan, *Descartes' rule of signs for the topological degree* ([arXiv:2611.34336](https://arxiv.org/abs/2611.34336)), gives |deg F| ≤ 13 here.

Two things to know before relying on it:
- It is **conditional on Bihan's Theorem 1.3** (the facet-counting formula for the degree), which is in recent preprints.
- The upper bound is **computer-assisted**: 1,343 small SAT instances, all unsatisfiable, with checked DRUP proofs (see below).

## Lower bound: a system of degree 6

Exponents (9 points in ℝ⁴):

```
a = (3,2,6,6) (4,5,4,2) (6,4,4,1) (3,0,3,3) (2,0,5,5) (6,6,2,6) (0,3,1,4) (3,1,5,2) (2,1,1,0)
```

Coefficients:

```
C = [ -2 -2  3 -2 -1 -1  3 -3  1 ]
    [  5 -7  1 -6 -1  0 -2  3  2 ]
    [  8 -5  0  2 -3 -4  1 -1  2 ]
    [  0  3  0  2  7 -2 -3 -8  1 ]
```

- Exact rational arithmetic: all maximal minors of both matrices are nonzero, the Newton polytope is simplicial with 9 vertices and 22 facets, and Bihan's formula gives **deg F = −6** for every generic u tried. Run `python3 verify_exact.py data/system_deg6.json`.
- Independently of that formula: F has (at least) 8 positive zeros, with Jacobian signs −,−,−,−,−,+,−,−, which sum to −6. Run `python3 roots.py data/system_deg6.json x 3000`. Completeness of the list of zeros is not claimed.
- [data/D44_family.json](data/D44_family.json): the same exponents with seven coefficient matrices (entries in [−8, 8]) of degree 0, 1, …, 6. Run `python3 verify_family.py`.
- [data/deg6_systems.json](data/deg6_systems.json): degree-6 systems on nine different combinatorial types of 4-polytope, including a neighborly one where all 27 facets have the same sign.

## Upper bound: no degree 7

Under the hypotheses, the boundary of the Newton polytope is a simplicial 3-sphere K with V ≤ 9 vertices, and deg F is, up to sign, the degree of the map K → S³ that sends each facet to the spherical tetrahedron spanned by its four coefficient columns. For generic u,

    deg = Σ over facets G with u ∈ cone(C_G) of  ε_G · sign det(C_G).

1. For every vertex v, u = −C_v is generic, and no facet through v has −C_v in its cone. So the sum can be taken over the facets avoiding v, for each v separately.
2. Whether −C_v lies in cone(C_G) is a condition on the signs of 4×4 minors of C (Cramer's rule). So the degree is a function of the *chirotope* of C: the signs of its C(V,4) maximal minors.
3. A chirotope satisfies the three-term Grassmann–Plücker relations: [σab][σcd] − [σac][σbd] + [σad][σbc] = 0, so the three terms, all nonzero, do not have the same sign.
4. For a fixed sphere K and bound d, this gives a SAT formula: variables are the minor signs (126 for V = 9), clauses are the relations in (3) plus "the sum in (1) is at least d, for every vertex v". Every real C with degree ≥ d satisfies it. A reflection of C negates the degree, so unsatisfiable means no C has |degree| ≥ d on K.

Results ([omsat.py](omsat.py); at most 2,655 variables and 27,099 clauses per instance):

| V | 3-spheres | best found by search | SAT relaxation |
|---|---|---|---|
| 5 | 1 | 1 | ≥ 2 unsatisfiable |
| 6 | 2 | 1 | ≥ 2 unsatisfiable |
| 7 | 5 | 3 | ≥ 4 unsatisfiable for all |
| 8 | 39 | 4 | ≥ 5 unsatisfiable for all |
| 9 | 1,296 | 6 | **≥ 7 unsatisfiable for all 1,296** |

- V ≤ 8 reproduces D(4,1) = {0,1} and D(4,2) = {0,…,3} (Bihan) and D(4,3) = {0,…,4} (the n = 3 result in this repository plus Bihan's duality), by a different method.
- **Certificates.** [certify.py](certify.py) solves each instance with Glucose 4 with proof logging and checks the DRUP proof with [rupcheck.c](rupcheck.c), a 100-line forward RUP checker. For "≥ 7": 1,291 proofs verified, and 52 instances contain the empty clause outright (some vertex has fewer than 7 facets avoiding it). The longest proof has 168 lemmas. The proofs have not been run through drat-trim.
- **Other solvers.** CaDiCaL, Glucose and MiniSat agree on all 1,343 instances ([crosscheck.py](crosscheck.py)).
- **Controls.** With the minor signs fixed to those of a real matrix, the formula is satisfiable at the true |degree| and unsatisfiable one above: 450 random matrices ([control.py](control.py)) and the 16 exact systems above, of degree 0 to 6 ([control_exact.py](control_exact.py)).
- **The list of spheres.** [spheres.py](spheres.py) enumerates simplicial 3-spheres by bistellar flips. The result agrees, class by class, with F. Lutz's published list of triangulated 3-manifolds with up to 9 vertices (1, 2, 5, 39, 1,296 spheres).

## Independent review

A separate session reviewed the claim with its own code, written from the definitions: its own exact verifier, its own root search, and its own SAT encoder run directly on Lutz's file. It reproduced every number above, sphere by sphere. The report is [review/topdeg-D44-coldreview.md](review/topdeg-D44-coldreview.md) and its code is in [review/](review/).

## Reproduce

```
pip install pynauty python-sat numpy scipy mpmath
cd n4k4
python3 spheres.py 9          # enumerate 3-spheres with 5..9 vertices  (about 1 s)
python3 orient.py             # orient them
python3 omsat.py 9 7          # "degree >= 7" for all 1,296: every line UNSAT  (about 2 min)
gcc -O2 -o rupcheck rupcheck.c
for V in 5 6 7 8 9; do python3 certify.py $V 7; done      # checked proofs  (about 5 min)
python3 verify_exact.py data/system_deg6.json
python3 verify_family.py
python3 control_exact.py
```

[search.c](search.c) is the hill-climbing search that found the degree-6 maps; [glue4.py](glue4.py) checks that degrees add when two such maps are glued along a facet (degree 12, 18, 24 on 14, 19, 24 vertices, as maps of spheres; the corresponding polynomial systems have not been built).

Found with the AI assistant Claude (Anthropic).
