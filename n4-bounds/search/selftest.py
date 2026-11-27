import sys, itertools
sys.argv = ["walker.py", "9", "6", "99", "0", "0", "selftest"]
src = open("walker.py").read().split("t0 = time.time(); S = starts()")[0]
exec(src)
C9 = cyclic(9)
assert sat(C9, 6) is True and sat(C9, 7) is False, "SAT control failed"
print("SELFTEST OK", sys.version.split()[0])
