"""Search for a 'gadget': vertex map on a 7-vertex triangulated sphere with |deg| = 3 and exactly ONE triangle
whose orientation sign is opposite to deg. Enumerates all labelled triangulations of S^2 on V vertices."""
import numpy as np, itertools, sys, pickle, os
from n8_exhaustive import flips, orient
V=int(sys.argv[1]); want=int(sys.argv[2]); trials=int(sys.argv[3]); seed=int(sys.argv[4])
fn=f'tri{V}.pkl'
if os.path.exists(fn): Ts=pickle.load(open(fn,'rb'))
else:
    start=frozenset(frozenset((i,(i+1)%(V-2),a)) for i in range(V-2) for a in (V-2,V-1))
    seen={start}; fr=[start]
    while fr:
        nx=[]
        for T in fr:
            for U in flips(T):
                if U not in seen: seen.add(U); nx.append(U)
        fr=nx
    Ts=list(seen); pickle.dump(Ts,open(fn,'wb'))
print('labelled triangulations',len(Ts),flush=True)
F=2*V-4
O=np.array([orient(T) for T in Ts])
trip=list(itertools.combinations(range(V),3)); tid={t:i for i,t in enumerate(trip)}
idx=np.zeros(O.shape[:2],int); sgn=np.zeros(O.shape[:2])
for m in range(O.shape[0]):
    for f in range(F):
        a,b,c=O[m,f]; s=tuple(sorted((a,b,c))); idx[m,f]=tid[s]
        perm=[s.index(x) for x in (a,b,c)]; inv=sum(1 for i in range(3) for j in range(i+1,3) if perm[i]>perm[j])
        sgn[m,f]=-1 if inv%2 else 1
def solid(a,b,c):
    a,b,c=a/np.linalg.norm(a),b/np.linalg.norm(b),c/np.linalg.norm(c)
    return 2*np.arctan2(a@np.cross(b,c),1+a@b+b@c+c@a)
rng=np.random.default_rng(seed); found=0
for t in range(trials):
    C=rng.normal(size=(V,3))
    w=np.array([solid(C[i],C[j],C[k]) for (i,j,k) in trip])
    W=w[idx]*sgn                      # signed solid angle per triangle
    deg=np.round(W.sum(1)/(4*np.pi)).astype(int)
    s=np.sign(deg)[:,None]
    nopp=((np.sign(W)*s)<0).sum(1)   # triangles opposite to the degree's sign
    ok=np.where((np.abs(deg)==want)&(nopp==1))[0]
    if len(ok):
        j=int(ok[0]); np.save(f'gadget_V{V}_d{want}.npy',{'C':C,'T':O[j],'deg':int(deg[j])})
        print(f'trial {t}: FOUND gadget, deg={deg[j]}, triangles={len(O[j])}',flush=True); found=1; break
print('found' if found else 'not found')
