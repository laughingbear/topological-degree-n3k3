"""Certificates for the upper bound: for every simplicial 3-sphere with V = 5..9 vertices, write the CNF "degree >= 7"
(omsat.build), solve it with Glucose 4 with proof logging (DRUP), and check the proof with our own forward RUP checker
(rupcheck.c, ~100 lines, independent of the solver).  usage: certify.py V [d] [keep]"""
import itertools, sys, os, subprocess, hashlib
from pysat.solvers import Glucose4
from omsat import build
V = int(sys.argv[1]); d = int(sys.argv[2]) if len(sys.argv) > 2 else 7; keep = len(sys.argv) > 3
sub = list(itertools.combinations(range(9), 4)); os.makedirs('cert', exist_ok=True)
ok = triv = fail = 0; h = hashlib.sha256(); maxlem = 0
for i, l in enumerate(open(f'pool{V}.txt').read().split('\n')[1:]):
    if not l.strip(): continue
    z = list(map(int, l.split())); ofac = [(sub[z[1 + 2 * k]], z[2 + 2 * k]) for k in range(z[0])]
    cnf, X = build(V, ofac, d)
    if any(len(c) == 0 for c in cnf.clauses): triv += 1; continue     # bound exceeds the number of available facets
    base = f'cert/V{V}_s{i}_d{d}'
    cnf.to_file(base + '.cnf')
    with Glucose4(bootstrap_with=cnf.clauses, with_proof=True) as s:
        r = s.solve(); proof = s.get_proof() if not r else None
    if r: print('SAT', V, i); fail += 1; continue
    open(base + '.drup', 'w').write('\n'.join(proof) + '\n'); maxlem = max(maxlem, len(proof))
    out = subprocess.run(['./rupcheck', base + '.cnf', base + '.drup'], capture_output=True, text=True).stdout
    if out.startswith('s VERIFIED'): ok += 1
    else: fail += 1; print('FAIL', V, i, out.strip())
    h.update(open(base + '.cnf', 'rb').read())
    if not keep: os.remove(base + '.cnf'); os.remove(base + '.drup')
print(f'V={V} d>={d}: verified {ok}, trivially UNSAT (empty clause in CNF) {triv}, failed {fail}, longest proof {maxlem} lines, sha256 of all CNFs {h.hexdigest()[:16]}')
