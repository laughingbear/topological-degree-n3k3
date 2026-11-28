r"""Which cell job of D_exclude.py contains each certified zero (bookkeeping for REVIEW 4.3)."""
import json, itertools, numpy as np
roots = [[float(x) for x in r] for r in json.load(open("out/D_roots_all.json"))]
R0 = 0.0625
jobs = [('in',)]
for face in range(4):
    for sign in (-1, 1):
        for sg in itertools.product((-0.5, 0.5), repeat=3): jobs.append(('shell', face, sign, sg))
        jobs.append(('far', face, sign))
for t in roots:
    t = np.array(t); m = np.abs(t).max()
    if m <= R0: idx = 0
    else:
        face = int(np.abs(t).argmax()); sign = int(np.sign(t[face])); nu = t / m
        others = [k for k in range(4) if k != face]
        sg = tuple(0.5 if nu[o] > 0 else -0.5 for o in others)
        idx = jobs.index(('shell', face, sign, sg))
    print(idx, "odd" if idx % 2 else "even", np.round(t, 4))
