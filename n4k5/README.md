# n = 4, k = 5 (and n = 5, k = 4): the maximal topological degree is 8

**Result (2026-11-26).** D(4,5) = {0, 1, …, 8}, and by Bihan's duality D(5,4) = {0, 1, …, 8}.

That is: for fewnomial maps F : ℝ⁴₊ → ℝ⁴ with N = 10 monomials whose exponent and coefficient matrices have nonzero maximal minors, the possible values of |deg F| are exactly 0, …, 8. Notation and background: F. Bihan, *Descartes' rule of signs for the topological degree* ([arXiv:2611.34336](https://arxiv.org/abs/2611.34336)); method as in [../n4k4](../n4k4/README.md).

**A side result.** Over *all* triangulated 3-spheres with 10 vertices, the largest degree of such a map is 9, not 8. It is attained on exactly one sphere, and that sphere is not the boundary of any 4-polytope, so it does not come from a polynomial system. For N = 10 the answer therefore depends on which spheres are polytopes, not only on the combinatorics of spheres.

Two things to know before relying on this:
- It is **conditional on Bihan's Theorem 1.3** (the facet-counting formula for the degree) and, for D(5,4), on his Theorem 3 (duality). Both are in recent preprints.
- The upper bound is **computer-assisted**: 247,882 SAT instances, with checked DRUP proofs for every unsatisfiable one (see below). The checker is our own small program, not a formally verified one.

## Lower bound: a system of degree 8

Exponents (10 points in ℝ⁴):

```
a = (16,0,20,22) (35,17,8,4) (32,6,16,9) (22,24,11,30) (18,17,0,0)
    (20,14,32,14) (18,17,3,27) (34,30,25,27) (25,35,21,8) (0,11,3,19)
```

Coefficients:

```
C = [ -4  0  3  2 -2  2 -7 -2 -2  6 ]
    [  1 -3  2  6  5 -3 -3 -3  2 -8 ]
    [ -3 -3  0  6  4  0  0  3 -1 -3 ]
    [  4  0  2 -1  6 -5 -7  0  1 -2 ]
```

- Exact rational arithmetic: all maximal minors of both matrices are nonzero, the Newton polytope is simplicial with 10 vertices and 29 facets, and Bihan's formula gives **deg F = 8** for every generic u tried. File: [system_deg8.json](system_deg8.json).
- Independently of that formula, F has exactly 10 positive zeros, nine with positive Jacobian determinant and one with negative, so the signs sum to 8:

```
sign   x1          x2          x3           x4
 +     61.568781   1.2987929   44.069953    0.066049984
 +     4.8623292   19.354202   0.037323859  1.4019541
 +     1.048859    0.88457213  1.9991066    0.20463232
 +     1.3190655   0.86007838  0.88234585   1.1697061
 +     1.8152427   0.13761143  0.14442304   1.2706193
 +     0.45613304  1.128726    2.040299     2.6227826
 +     0.42321851  2.1397338   1.649363     1.5245866
 -     0.94701853  1.012609    1.0753688    0.88909773
 +     0.88654779  1.0447774   1.092965     0.86213811
 +     0.67029418  0.52682446  1.4433182    0.42625089
```

  The list is from the independent review (branch and bound in logarithmic coordinates with a Krawczyk test and an explicit bound on where zeros can lie; no unresolved box; double precision with slack, so strong evidence of completeness and not a formal interval proof). `python3 roots_trop.py system_deg8.json 3` finds the same ten zeros by Newton's method.
- [D45_family.json](D45_family.json): the same exponents with nine coefficient matrices of degree 0, 1, …, 8.
- Degree 8 is attained on all 22 polytope types that allow it; realizations of the 22 are in [realizations_22.json](realizations_22.json), and systems for all of them are in the review package.

## Upper bound

The boundary of the Newton polytope is a simplicial 3-sphere with at most 10 vertices. Fewer than 10 vertices gives |deg F| ≤ 6 by the N = 9 result. For 10 vertices:

1. **All spheres.** [enum10.py](enum10.py) enumerates the simplicial 3-spheres with 10 vertices by bistellar flips: 247,882, the same list, class by class, as F. Lutz's.
2. **Degree ≥ 8** ([run10.py](run10.py), same SAT formula as in n4k4): unsatisfiable for 247,843 of them. It is satisfiable for 39 ([sat8_idx.json](sat8_idx.json)).
3. **Degree ≥ 9** is satisfiable for exactly one of the 39, index 101272; **degree ≥ 10** is unsatisfiable there.
4. **The 39.** One is sphere 101272. Of the other 38, 16 are not polytopal (a rank-5 chirotope test, [polysat.py](polysat.py), is already unsatisfiable) and 22 are polytopal, with degree at most 8.
5. **Sphere 101272** ([sphere_deg9.json](sphere_deg9.json); Lutz's `manifold_3_10_249002`; 30 facets, every vertex in 12 of them, automorphism group of order 20) is not polytopal. This is known: it is absent from M. Firsching's classification of the 162,004 polytopal simplicial 3-spheres with 10 vertices ([arXiv:1508.02531](https://arxiv.org/abs/1508.02531)). We give an independent certificate of the same kind: exactly 4 uniform rank-5 chirotopes are compatible with its face lattice ([count27.py](count27.py), [chi27.json](chi27.json)), and each is non-realizable by a biquadratic final polynomial, verified in exact integer arithmetic ([bfp.py](bfp.py), [bfp27.json](bfp27.json)).
6. **Degree 9 on that sphere.** [vertexmap_deg9_nonpolytopal.json](vertexmap_deg9_nonpolytopal.json): integer vectors, every four independent, exact degree 9, all 30 facets of the same sign (a 9-sheeted branched cover of S³).

**Certificates.** [certify10.py](certify10.py) solves every instance again with a fresh Glucose 4 with proof logging and streams each DRUP proof to [rupserve.c](rupserve.c), a forward RUP checker.

| question | proofs verified | empty clause outright | satisfiable | failed |
|---|---|---|---|---|
| degree ≥ 8, all 247,882 spheres | 247,712 | 131 | 39 | 0 |
| degree ≥ 9, the 39 | 38 | 0 | 1 | 0 |
| degree ≥ 10, sphere 101272 | 1 | 0 | 0 | 0 |

Per-sphere outcomes: [certificates_ge8_results.txt.gz](certificates_ge8_results.txt.gz). The step "exactly 4 chirotopes" also has a checked proof. As negative controls, 240 truncated proofs or weakened formulas were all rejected.

## Independent review

A separate session reviewed both claims with its own code: its own exact verifier, a complete root search, its own SAT encoder run on Lutz's file for all 247,882 spheres (identical outcome on every sphere), 598 positive controls up to degree 9, its own enumeration of the 4 chirotopes with five solvers, an exact re-check of the four final polynomials, exact convex hulls of the 22 realizations, and a comparison with Firsching's data. Report: [review/topdeg-D45-coldreview.md](review/topdeg-D45-coldreview.md); code and results: [review/topdeg-D45-coldreview-code.tar.gz](review/topdeg-D45-coldreview-code.tar.gz). The corrections it asked for are incorporated here.

## Reproduce

```
pip install pynauty python-sat numpy scipy mpmath sympy
cd n4k5
python3 verify_all.py                    # exact checks without a SAT solver: systems, degree-9 map, final polynomials  (1 min)
python3 enum10.py                        # the 247,882 spheres -> spheres10.json  (6 min, 100 MB)
python3 orient10.py                      # -> pool10_all.txt
python3 make_pool39.py                   # -> pool10_sat8.txt
gcc -O2 -o rupserve rupserve.c
python3 certify10.py 0 39 8 c39_8.txt pool10_sat8.txt      # all 39 SAT
python3 certify10.py 0 39 9 c39_9.txt pool10_sat8.txt      # 38 verified, 1 SAT
python3 certify10.py 27 28 10 c_10.txt pool10_sat8.txt     # verified
python3 certify10.py 0 247882 8 c_all.txt                  # the full sweep: about 3.5 CPU-hours; resumable, ranges can run in parallel
python3 count27.py                       # the 4 chirotopes of sphere 101272
```

## Related work

- M. Firsching, *Realizability and inscribability for simplicial polytopes via nonlinear optimization*, Math. Program. 166 (2017); data at github.com/mo271/polytope-data.
- F. H. Lutz, *Combinatorial 3-manifolds with 10 vertices*, [arXiv:math/0604018](https://arxiv.org/abs/math/0604018).
- Simplicial maps of triangulated spheres onto the boundary of a simplex, a special degenerate case of the maps used here: Ryabichev [arXiv:2512.01137](https://arxiv.org/abs/2512.01137); Basak, Gupta, Trivedi [arXiv:2409.00907](https://arxiv.org/abs/2409.00907); [arXiv:2511.10870](https://arxiv.org/abs/2511.10870).

Maximal degree for n = 4 so far: 1, 1, 3, 4, 6, 8 for N = 5, …, 10.

Found with the AI assistant Claude (Anthropic).
