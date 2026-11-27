import numpy as np, itertools, sys
from mixed_check import orient
rng=np.random.default_rng(0)
def cyc_facets(V):
    F=set()
    for i in range(V):
        for j in range(V):
            f={i,(i+1)%V,j,(j+1)%V}
            if len(f)==4: F.add(tuple(sorted(f)))
    return sorted(F)
def degree(facets,o,C,nu=5):
    M=np.array([C[:,list(f)] for f in facets]); d=np.linalg.det(M)
    if np.min(np.abs(d))<1e-9: return None
    s=np.sign(d)*np.array(o); Mi=np.linalg.inv(M); res=set()
    for _ in range(nu):
        u=rng.standard_normal(4); c=Mi@u
        res.add(int(round((s*(c>0).all(1)).sum())))
    return res.pop() if len(res)==1 else None
if __name__ == "__main__":
  for V in [9,10,11,12,13,15,17,19,21,25,29]:
      facets,o=orient(cyc_facets(V)); best=(0,None); i=np.arange(V)
      for r in range(1,V):
          for s in range(r+1,V):
              C=np.array([np.cos(2*np.pi*r*i/V),np.sin(2*np.pi*r*i/V),np.cos(2*np.pi*s*i/V),np.sin(2*np.pi*s*i/V)])
              C=C+1e-7*rng.standard_normal(C.shape)
              d=degree(facets,o,C)
              if d is not None and abs(d)>best[0]: best=(abs(d),(r,s))
      print(V,"F=",len(facets),"best",best,"V^2/16=",V*V/16,"V^2/8=",V*V/8, "ratio",best[0]/V**2, flush=True)
