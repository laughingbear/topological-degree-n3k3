# Referee code for `topdeg-review-2026-11-27`

All files here were written by the referee from the definitions in `CLAIMS.md`.
Nothing is imported from `ours/`; the stored `orientation` arrays are only *compared* with the
recomputed ones, never used.

Requirements: Python 3.10+. Confirming checks need only the standard library.
Optional: `python-sat` (`sat_cross.py`), `mpmath` (`D_seek.py`, `D_refine.py`, `D_certify.py`),
`numpy` + `scipy` (`D_seek_fast.py`, `D_exclude.py`).

```
CERTS=/path/to/topdeg-review-2026-11-27/certs
./run_all.sh $CERTS            # everything except the long Claim-D searches (about 6 minutes)
```

| file | claim | arithmetic | what it does |
|---|---|---|---|
| `common.py` | all | exact int | determinant, closed-manifold test, coherent orientation, facet signs, degree |
| `sphere_cert.py` | A (i), C | exact int | finds a **shelling** and re-verifies it from the definition; integer homology by Smith normal form |
| `check_A.py` | A (ii),(iii) | exact int | recomputed orientation, degree at 64 generic u, q, Theorem M |
| `controls_A.py` | A | exact int | random maps give u-independent degree; flipping one facet sign breaks it |
| `hull.py`, `check_CD.py` | C, D | exact int | convex hull by brute force over 4-subsets, (U_A), (U_C), Bihan's formula at 305 u |
| `chiro.py` | B | combinatorial | own encoder + unit propagator, writes minimal proofs to `out/proof_*.json` |
| `check_proof.py` | B | combinatorial | stand-alone proof checker (shares no code with `chiro.py`); `--print` gives a readable proof |
| `b_human.py` | B | combinatorial | hand-checkable form of the 6-step proofs (V = 12, 16, 20) |
| `chiro_controls.py` | B | exact int | real polytopes satisfy every encoded constraint; propagator derives only true signs |
| `sat_cross.py` | B | SAT | CaDiCaL on the CNF of our encoder (V = 12, 13, 14 + SAT control) |
| `check_sym.py`, `e1_enum.py` | E1 | combinatorial | symmetry of the certificates; independent count for n = 5 |
| `D_seek.py`, `D_seek_fast.py`, `D_refine.py` | D | floating / multiprecision | *discovery only* of positive zeros |
| `D_certify.py` | D | **exact rational intervals** | Krawczyk existence+uniqueness and sign of det DF for each zero |
| `D_exclude.py`, `D_union.py`, `D_assign.py` | D | float64 | branch-and-bound completeness search over all of t-space (numerical evidence, not a proof); 73 jobs, about 75 CPU-minutes: `python3 D_exclude.py system.json 512 0 2 0` (central cube), `... 512 0 2` (even jobs), `... 512 1 2` (odd jobs) |

`out/` contains the artefacts of the referee's run: shellings, proofs, certified roots, logs.
