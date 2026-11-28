#!/bin/sh
# usage: ./run_all.sh /path/to/certs      (run from the directory containing this file)
set -e
D="$1"
mkdir -p out
echo "== A (i): shelling + integer homology"
python3 sphere_cert.py "$D"/V*.json
echo "== A (ii),(iii) and degree part of C"
python3 check_A.py "$D"/V*.json
echo "== A controls"
python3 controls_A.py "$D"/V*.json
echo "== C and D (combinatorial part)"
python3 check_CD.py "$D"/V12_polytope_deg11.json "$D"/system_N12_deg11.json
echo "== B: encoder controls on real polytopes"
python3 chiro_controls.py "$D"/V12_polytope_deg11.json "$D"/system_N12_deg11.json
echo "== B: unit propagation with proof logging"
python3 chiro.py "$D"/V12_deg13.json "$D"/V13_deg16.json "$D"/V14_deg17.json "$D"/V16_deg25.json "$D"/V18_deg35.json "$D"/V20_deg41.json
echo "== B: independent proof checker"
python3 check_proof.py out/proof_*.json
python3 check_proof.py --print out/proof_eps+1_V12_deg13.json
python3 b_human.py out/proof_eps+1_V12_deg13.json
echo "== E1"
python3 check_sym.py "$D"/V*.json
python3 e1_enum.py 5 0
python3 e1_enum.py 5 1 0
echo "== D: exact certification of the zeros listed in out/D_roots_all.json (slow, ~4 min)"
python3 D_certify.py "$D"/system_N12_deg11.json out/D_roots_all.json
echo "Optional (not run here): sat_cross.py, D_seek.py, D_seek_fast.py, D_refine.py, D_exclude.py"
