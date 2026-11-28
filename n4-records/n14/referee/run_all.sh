#!/bin/sh
# Re-run every referee check.  Usage: sh run_all.sh <path to unpacked topdeg-review-N14-2026-11-27>
# Needs python3 (numpy, mpmath) and g++.  Runtime about 3 minutes.  Logs go to ./logs/.
set -e
D="$1"; H="$(cd "$(dirname "$0")" && pwd)"; mkdir -p "$H/logs" "$H/build"
P="$D/certs/V14_polytope_deg16.json"; S="$D/certs/system_N14_deg16.json"
python3 -I "$H/check_P.py" "$P" "$S" 2000            | tee "$H/logs/check_P.log"
python3 -I "$H/S_outer.py" "$S" 8                     | tee "$H/logs/S_outer.log"
python3 -I "$H/mk_system_txt.py" "$S" > "$H/build/system.txt"
g++ -O2 -std=c++17 -o "$H/build/S_inner" "$H/S_inner.cpp"
"$H/build/S_inner" 8 < "$H/build/system.txt"          | tee "$H/logs/S_inner.log"
python3 -I "$H/S_certify.py" "$S" "$H/logs/S_inner.log" | tee "$H/logs/S_certify.log"
python3 -I "$H/E_checks.py" "$P" "$H/system_N14_deg16_small.json" | tee "$H/logs/E_checks.log"
python3 -I "$H/E4_scan.py" "$P" | tee "$H/logs/E4_scan.log"
