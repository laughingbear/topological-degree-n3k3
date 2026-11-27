"""V = 10 sweep (towards D(4,5) = D(5,4)): for spheres start..end-1 of pool10_all.txt decide the relaxation "degree >= d".
Incremental: one CaDiCaL instance holds the Grassmann-Pluecker clauses; each sphere's clauses (same as omsat.build, see
omsat_fast.check_same) get fresh auxiliary variables and an activation literal, are solved under that assumption and then
switched off.  The solver is rebuilt every BATCH spheres.  Output lines: "index SAT|UNSAT|TRIV".
usage: run10.py start end d out.txt"""
import itertools, sys, time
from pysat.solvers import Cadical153
from omsat_fast import build_fast, _gp
V = 10; BATCH = 400
start, end, d, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
sub = list(itertools.combinations(range(V), 4)); NX = len(sub)
lines = open(sys.argv[5] if len(sys.argv) > 5 else 'pool10_all.txt').read().split('\n')[1:]
done = set()
try: done = {int(l.split()[0]) for l in open(outp) if l.strip()}
except FileNotFoundError: pass
out = open(outp, 'a'); t0 = time.time(); solver = None; nin = 0; top = 0; cnt = {'SAT': 0, 'UNSAT': 0, 'TRIV': 0}
for i in range(start, min(end, len(lines))):
    if i in done or not lines[i].strip(): continue
    z = list(map(int, lines[i].split())); ofac = [(sub[z[1 + 2 * k]], z[2 + 2 * k]) for k in range(z[0])]
    r = build_fast(V, ofac, d)
    if r[0] is None: res = 'TRIV'
    else:
        gp, cl, X = r
        if solver is None or nin >= BATCH:
            if solver is not None: solver.delete()
            solver = Cadical153(bootstrap_with=gp); nin = 0; top = NX
        mx = max(abs(l) for c in cl for l in c); off = top - NX; act = mx + off + 1; top = act
        sh = lambda l: l if abs(l) <= NX else (l + off if l > 0 else l - off)
        for c in cl: solver.add_clause([sh(l) for l in c] + [-act])
        res = 'SAT' if solver.solve(assumptions=[act]) else 'UNSAT'
        solver.add_clause([-act]); nin += 1
    cnt[res] += 1; out.write(f'{i} {res}\n')
    if (i - start) % 300 == 0: out.flush()
out.close(); print('range', start, end, 'd', d, cnt, round(time.time() - t0), 's', flush=True)
