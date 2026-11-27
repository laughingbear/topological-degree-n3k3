"""Checked proofs for the V = 10 sweep.  For spheres start..end-1 of the pool: build "degree >= d" (omsat_fast.build_fast, identical
to omsat.build), solve with a FRESH Glucose 4 with DRUP proof logging, and have the proof checked by rupserve (our forward RUP checker,
independent of the solver).  Output lines: "index VERIFIED n" | "index TRIV" | "index SAT" | "index FAIL ...".
usage: certify10.py start end d out.txt [pool]"""
import itertools, sys, subprocess, time, os
from pysat.solvers import Glucose4
from omsat_fast import build_fast, _gp
V = 10
start, end, d, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
pool = sys.argv[5] if len(sys.argv) > 5 else 'pool10_all.txt'
sub = list(itertools.combinations(range(V), 4))
lines = open(pool).read().split('\n')[1:]
X, gp = _gp(V)
base = f'gp{V}.cnf'
if not os.path.exists(base):
    with open(base + '.tmp%d' % os.getpid(), 'w') as f:
        f.write(f'p cnf {len(X)} {len(gp)}\n' + ''.join(' '.join(map(str, c)) + ' 0\n' for c in gp))
    os.replace(base + '.tmp%d' % os.getpid(), base)
srv = subprocess.Popen(['./rupserve', base], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
done = set()
try: done = {int(l.split()[0]) for l in open(outp) if l.strip()}
except FileNotFoundError: pass
out = open(outp, 'a'); t0 = time.time(); cnt = {}
for i in range(start, min(end, len(lines))):
    if i in done or not lines[i].strip(): continue
    z = list(map(int, lines[i].split())); ofac = [(sub[z[1 + 2 * k]], z[2 + 2 * k]) for k in range(z[0])]
    r = build_fast(V, ofac, d)
    if r[0] is None: res = 'TRIV'
    else:
        _, cl, _ = r
        with Glucose4(bootstrap_with=gp, with_proof=True) as s:
            s.append_formula(cl); sat = s.solve(); proof = None if sat else s.get_proof()
        if sat: res = 'SAT'
        else:
            srv.stdin.write(f'J {i}\n' + ''.join(' '.join(map(str, c)) + ' 0\n' for c in cl) + 'P\n' + '\n'.join(proof) + '\nE\n'); srv.stdin.flush()
            ans = srv.stdout.readline().split(); assert ans and int(ans[0]) == i, ans
            res = ' '.join(ans[1:])
    k = res.split()[0]; cnt[k] = cnt.get(k, 0) + 1; out.write(f'{i} {res}\n')
    if (i - start) % 200 == 0: out.flush()
out.close(); srv.stdin.close(); print('range', start, end, 'd', d, cnt, round(time.time() - t0), 's', flush=True)
