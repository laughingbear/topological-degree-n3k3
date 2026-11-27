#!/bin/bash
# 003-sphsearch: overnight search for simplicial 3-spheres whose maximal degree exceeds the cyclic value
#   V=11: is there a sphere with "deg >= 11"?   V=12: how high does the maximum go?   V=13: is there one with "deg >= 16"?
# Uses the chirotope SAT relaxation (same encoder as D(4,4), D(4,5)).  python-sat is bundled (no install, no network).
set -e
PV=$(python3 -c 'import sys; print("cp%d%d" % sys.version_info[:2])')
echo "python tag $PV"
[ -f "pydeps-$PV.tar.gz" ] || { echo "NO BUNDLED python-sat FOR $PV"; exit 3; }
mkdir -p pydeps && tar xzf "pydeps-$PV.tar.gz" -C pydeps
export PYTHONPATH="$PWD/pydeps/$PV"
python3 selftest.py || { echo "SELFTEST FAILED"; exit 2; }
SECS=${SECS:-25200}
python3 walker.py 11 10 11 $SECS 101 w11a &
python3 walker.py 11 10 11 $SECS 102 w11b &
python3 walker.py 11 10 11 $SECS 103 w11c &
python3 walker.py 12 10 11 $SECS 201 w12a &
python3 walker.py 12 10 11 $SECS 202 w12b &
python3 walker.py 13 15 16 $SECS 301 w13a &
python3 walker.py 13 15 16 $SECS 302 w13b &
python3 walker.py 13 15 16 $SECS 303 w13c &
wait
echo "ALL WALKERS DONE"
for f in w*.log; do echo "== $f"; tail -n 1 "$f"; done
echo "hits:"; cat w*.hits.jsonl 2>/dev/null | python3 -c "
import sys, json
L = [json.loads(l) for l in sys.stdin]
for V in (11, 12, 13):
    gs = sorted(x['g'] for x in L if x['V'] == V)
    print('V', V, 'records/hits', len(gs), 'max g', gs[-1] if gs else None)"
