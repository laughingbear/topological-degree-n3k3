"""usage: run.py V d [solver]   -> SAT/UNSAT of 'deg >= d' for every Lutz 3-sphere with V vertices (own encoder)."""
import sys, json, time
from multiprocessing import Pool
from pysat.solvers import Solver
from sph import parse_lutz, orient
from enc import Enc
V = int(sys.argv[1]); d = int(sys.argv[2]); sv = sys.argv[3] if len(sys.argv) > 3 else 'cadical153'
L = parse_lutz('/home/claude/ext/3_manifolds_all.txt', '/home/claude/ext/3_manifolds_all_type.txt')
names = sorted((n for n, (v, F, t) in L.items() if v == V and t == 'S^3'), key=lambda n: int(n.split('_')[3]))
E = None
def work(n):
    global E
    if E is None: E = Enc(V)
    cl, nv = E.formula(orient(L[n][1]), d)
    if any(len(c) == 0 for c in cl): return n, False, nv, len(cl)
    with Solver(name=sv, bootstrap_with=cl) as s: return n, s.solve(), nv, len(cl)
if __name__ == '__main__':
    t = time.time()
    with Pool(2) as p: res = p.map(work, names, chunksize=8)
    sat = [n for n, r, _, _ in res if r]
    json.dump({n: r for n, r, _, _ in res}, open(f'res_{V}_{d}_{sv}.json', 'w'))
    print(f'V={V} d>={d} {sv}: instances {len(res)}, SAT {len(sat)}, UNSAT {len(res) - len(sat)}, max vars {max(x[2] for x in res)}, max clauses {max(x[3] for x in res)}, {time.time() - t:.0f}s', ('SAT: ' + ' '.join(sat)) if len(sat) <= 30 else '')
