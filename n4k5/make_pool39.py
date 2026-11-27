"""pool10_sat8.txt: the 39 spheres on which "degree >= 8" is satisfiable (indices in sat8_idx.json), from pool10_all.txt."""
import json
lines = open('pool10_all.txt').read().split('\n')[1:]; S = json.load(open('sat8_idx.json'))
open('pool10_sat8.txt', 'w').write(f'{len(S)}\n' + '\n'.join(lines[i] for i in S) + '\n'); print(len(S), 'spheres')
