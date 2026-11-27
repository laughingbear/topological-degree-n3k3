# Referee code for the topdeg n = 4 claims

Written from the definitions in CLAIMS.md §0 alone. Nothing is imported from the attached scripts.
Needs Python 3, numpy, scipy and mpmath.

| file | what it does |
|---|---|
| `core.py` | My own coherent orientation (BFS over triangles), exact integer determinants (Bareiss) and Cramer signs, degree by cone counting, the **link degree w_v computed as a genuine 3-dimensional problem in C_v^⊥** (induced orientation, projected vectors, cone counting), the **edge winding w_e as an angle sum** of the oriented projected link polygon, closed-form 3-d solid angles, and 4-d orthant probabilities. |
| `cyclic.py` | Item 5(a). Builds the facets of C(V,4) from **Gale's evenness condition** and checks them against the recipe {i,i+1,j,j+1} and against a Qhull hull of points on the trigonometric moment curve. Orientation is my own. Computes C_i = (−1)^i γ(πi/V) in 60-digit arithmetic and the degree by cone counting at 5 random u and 3 values u = γ(τ), for odd V = 5..41. |
| `randmaps.py` | Item 5(b). Random simplicial 4-polytopes with 6..14 vertices (integer coordinates, simpliciality checked exactly) and three families of integer vertex maps. Also runs the example JSON files, read directly. For every map it checks: F = E − V; the exact degree; exact w_v; w_e; the per-vertex and global Gauss–Bonnet relations; w_e ≤ (p_e−1)/2; Σ_{t∋v} s_t μ(L_{t,v}) = w_v (closed-form solid angles); Σ s_t μ(cone t) = deg; the per-facet relation Σ_v μ(L_{t,v}) = 4μ(cone t) + D_t; the identity 4·deg = Σw_v − X + Y; 0 < D_t < 1/2; and 4·deg < 2F − E/2 − q/2. |
| `search.py` | Item 5(c). Hill-climbing and annealing on 4·deg − (2F − E/2 − q/2) over vertex maps on cyclic and random polytopes, seeded partly from the Theorem C maps. Tracks the same-sign refinement. Every record map is re-verified exactly and run through the full `randmaps.analyse`. Also probes D_t < 1/2 with random bases, Nelder–Mead, and an explicit family where D_t → 1/2. |
| `logs/` | Output of the runs reported in review.md: `cyclic.log`, `examples.log`, `rand_A.log`, `rand_B2.log`, `search.log`, `search_Dprobe.log`, and `att_*.log` for the attached scripts. `search.log` ends with a traceback from the abandoned Nelder–Mead stage, which crashed after all its results had been printed. |

Run:

```
python3 -I cyclic.py
python3 -I randmaps.py path/to/topdeg4 105 111   # examples + 105 random maps (seed 111)
python3 -I randmaps.py - 105 333                 # 105 more random maps (seed 333)
python3 -I search.py 800                         # 800-second search budget (log: logs/search.log)
python3 -I dprobe_family.py                      # D_t -> 1/2 family (log: logs/search_Dprobe.log)
```

Precision: signs, degrees and w_v are exact (integer arithmetic). w_e is checked to be within 1e-7 of an integer. The 4-d orthant probabilities come from Genz's algorithm (scipy) at 1e-6 tolerance, so the identity is verified to about 1e-5, not exactly.
