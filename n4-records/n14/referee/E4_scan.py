#!/usr/bin/env python3
"""Joint scan of exponent scale 4..26 and coefficient scale 2..20 (exact verification of each candidate).
Usage: python3 -I E4_scan.py <V14_polytope_deg16.json>"""
import sys,json,math,os
sys.argv=[sys.argv[0],sys.argv[1]]
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import io,contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import E_checks as E
from itertools import combinations
hits=[]
for sa in range(4,27):
    Ai=E.rounded(E.Xf,E.mx,sa)
    if len(set(map(tuple,Ai)))<14: continue
    r0=E.analyse(Ai,E.Cfull,nu=3)
    if not(r0['UA'] and r0['nvert']==14): continue
    for sc in range(2,21):
        Ci=E.rounded(E.Yf,E.my,sc)
        if E.good(Ai,Ci): hits.append((sa,sc)); 
print('A-scales with (U_A) and 14 vertices and some working C scale, (sa,sc) pairs:',hits)
