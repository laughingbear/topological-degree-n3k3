"""Exact certificate for the gluing gadget: an oriented 7-vertex triangulated sphere T and integer vectors C_v in Z^3
with all triples independent; degree of the vertex map computed exactly by counting, for several generic u,
the oriented triangles whose cone contains u (signed by det).  Also lists the triangle orientation signs."""
import json, sys, itertools
from fractions import Fraction as Fr
def det3(a,b,c): return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))
def deg_exact(T, C, u):
    tot=0
    for a,b,c in T:
        d=det3(C[a],C[b],C[c])
        l=[Fr(det3(u,C[b],C[c]),d), Fr(det3(C[a],u,C[c]),d), Fr(det3(C[a],C[b],u),d)]
        if all(x>0 for x in l): tot += 1 if d>0 else -1
    return tot
if __name__=='__main__':
    G=json.load(open(sys.argv[1])); T=[tuple(t) for t in G['T']]; C=G['C']; V=len(C)
    # sphere check: every edge in exactly two triangles with opposite directions
    from collections import Counter
    E=Counter((a,b) for t in T for a,b in [(t[0],t[1]),(t[1],t[2]),(t[2],t[0])])
    assert all(E[(b,a)]==1 for (a,b) in E) and max(E.values())==1 and V-len(E)//2+len(T)==2, 'not an oriented sphere'
    assert all(det3(C[i],C[j],C[k])!=0 for i,j,k in itertools.combinations(range(V),3)), 'dependent triple'
    for u in ([1,2,3],[-3,1,7],[5,-2,-1],[11,13,-17],[2,-7,5]):
        print('u',u,'deg',deg_exact(T,C,u))
    print('triangle signs', [ (1 if det3(C[a],C[b],C[c])>0 else -1) for a,b,c in T])
