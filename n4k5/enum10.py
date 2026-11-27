import sys, json, time
from spheres import enumerate_all
t=time.time(); L=enumerate_all(10); print(len(L), 'spheres, max facets', max(len(F) for F in L), round(time.time()-t), 's', flush=True)
json.dump(L, open('spheres10.json','w'))
