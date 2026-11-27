import itertools, json
from pysat.solvers import Cadical153
from polysat import gp5
V=10; S=json.load(open('sphere_deg9.json')); sub=list(itertools.combinations(range(10),4))
lines=open('pool10_all.txt').read().split('\n')[1:]; z=list(map(int,lines[S['index']].split())); ofac=[(sub[z[1+2*k]],z[2+2*k]) for k in range(z[0])]
X,cl,chi=gp5(V); units=[]
for G,eps in ofac:
    for w in range(V):
        if w in G: continue
        l=chi(tuple(G)+(w,)); units.append([l if eps>0 else -l])
fixed={abs(u[0]) for u in units}; free=[X[s] for s in X if X[s] not in fixed]; print('fixed',len(fixed),'free',len(free))
models=[]
with Cadical153(bootstrap_with=cl+units) as s:
    while len(models)<3000 and s.solve():
        m=set(s.get_model()); a=[v if v in m else -v for v in free]; models.append(a); s.add_clause([-l for l in a])
print('chirotopes compatible with the face lattice:',len(models),'(capped at 3000)' if len(models)>=3000 else '(complete)')
json.dump({'free':free,'models':models,'units':[u[0] for u in units]},open('chi27.json','w'))
