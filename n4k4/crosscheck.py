import itertools, sys, json
from pysat.solvers import Glucose4, Minisat22
from omsat import build
sub = list(itertools.combinations(range(9), 4)); tot = 0; sat = []; ncl = 0; nv = 0
for V in (5, 6, 7, 8, 9):
    for i, l in enumerate(open(f'pool{V}.txt').read().split('\n')[1:]):
        if not l.strip(): continue
        z = list(map(int, l.split())); ofac = [(sub[z[1 + 2 * k]], z[2 + 2 * k]) for k in range(z[0])]
        cnf, X = build(V, ofac, 7); tot += 1
        if any(len(c) == 0 for c in cnf.clauses): continue
        ncl = max(ncl, len(cnf.clauses)); nv = max(nv, cnf.nv)
        with (Glucose4 if sys.argv[1] == 'g' else Minisat22)(bootstrap_with=cnf.clauses) as s:
            if s.solve(): sat.append((V, i))
print(sys.argv[1], 'instances', tot, 'SAT:', sat, 'max vars', nv, 'max clauses', ncl)
