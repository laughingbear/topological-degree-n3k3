"""Discovery search (how the degree-3 example was found).
Triangulation: vertex 0 coned over hexagon 1..6; other side = triangle (1,3,5) + three ears.
Searches vertex images C on S^2 for a degree-3 vertex map (all triangles positively oriented).
Usage: python3 search_hexagon.py <seed> <tries>   -> writes hex_hit.npy on success.
"""
import numpy as np, sys
from deg import solid
cone=[(0,i,i%6+1) for i in range(1,7)]
bottom=[(1,5,3),(1,3,2),(3,5,4),(5,1,6)]
F=cone+bottom
from collections import Counter
E=Counter(); [E.update([(i,j),(j,k),(k,i)]) for i,j,k in F]
assert all(E[(b,a)]==1 for (a,b) in E) and max(E.values())==1, 'orientation'
rng=np.random.default_rng(int(sys.argv[1]))
def sph(lon,lat): return np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
best=0
def deg_and_neg(C):
    w=np.array([solid(C[i],C[j],C[k]) for i,j,k in F]); return w.sum()/4/np.pi, w
for t in range(int(sys.argv[2])):
    C=np.zeros((7,3)); C[0]=[0,0,1]
    for i in range(1,7):
        C[i]=sph(np.deg2rad(120*i)+rng.normal()*0.6, rng.uniform(-1.3,1.3))
    d,w=deg_and_neg(C)
    # local improve: minimise negative area
    for it in range(400):
        C2=C.copy(); v=rng.integers(0,7); C2[v]+=rng.normal(size=3)*0.2; C2[v]/=np.linalg.norm(C2[v])
        d2,w2=deg_and_neg(C2)
        key=lambda d,w:(round(abs(d)), -np.clip(-np.sign(d)*w,0,None).sum())
        if key(d2,w2)>=key(d,w): C,d,w=C2,d2,w2
    if round(abs(d))>best or round(abs(d))>=3:
        best=max(best,round(abs(d))); print(t,'deg',round(d,6),'tri angles/pi',np.round(w/np.pi,3),flush=True)
        if round(abs(d))>=3: np.save('hex_hit.npy',C); break
print('best',best)
