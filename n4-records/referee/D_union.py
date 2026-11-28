r"""Collect the zeros reported by all D_exclude.py logs given on the command line, compare with the
certified list out/D_roots_all.json, and summarise job coverage."""
import json, re, sys, numpy as np
cert = np.array([[float(x) for x in r] for r in json.load(open("out/D_roots_all.json"))])
found = []; jobs = set(); st = [0, 0]; unres = 0
for f in sys.argv[1:]:
    for line in open(f):
        if line.startswith("ZERO"):
            found.append([float(x) for x in re.search(r"\[(.*?)\]", line).group(1).split()])
        m = re.match(r"# job (\d+) ", line)
        if m: jobs.add(int(m.group(1)))
        if line.startswith("DONE"):
            m = re.search(r"selftest points/failures \[(\d+), (\d+)\]", line)
            if m: st[0] += int(m.group(1)); st[1] += int(m.group(2))
        unres += "UNRESOLVED" in line
found = np.array(found)
uniq = []
for z in found:
    if not any(np.abs(z - u).max() < 1e-7 for u in uniq): uniq.append(z)
match = [int(np.abs(cert - z).max(axis=1).argmin()) for z in uniq]
ok = all(np.abs(cert[m] - z).max() < 1e-8 for m, z in zip(match, uniq))
print(f"jobs completed: {len(jobs)} of 73 (missing: {sorted(set(range(73)) - jobs)}); unresolved cells: {unres}")
print(f"distinct zeros reported: {len(uniq)}; each equals a certified zero to 1e-8: {ok}; "
      f"certified zeros not reported: {sorted(set(range(len(cert))) - set(match))}")
print(f"self-test points / failures: {st}")
