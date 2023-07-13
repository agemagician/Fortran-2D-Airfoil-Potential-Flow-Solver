#!/bin/sh
# Runs the debug executable on the default input and on a bad input, in a temporary directory.
set -e
exe=build/debug/airfoil
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

sed -e "s|^ *vtk_path.*|  vtk_path = '$tmp/flow.vtk'|" -e "s|^ *cp_path.*|  cp_path = '$tmp/cp.csv'|" input.nml > "$tmp/good.nml"
"$exe" "$tmp/good.nml" > "$tmp/good.log"
head -1 "$tmp/flow.vtk" | grep -q '^# vtk DataFile Version 3.0' || { echo "FAIL: vtk header"; exit 1; }
pts=$(grep '^POINTS' "$tmp/flow.vtk" | awk '{print $2}')
[ "$pts" = "$((257 * 96))" ] || { echo "FAIL: POINTS $pts"; exit 1; }
rows=$(($(wc -l < "$tmp/cp.csv") - 1))
[ "$rows" = "256" ] || { echo "FAIL: csv rows $rows"; exit 1; }

sed -e 's/^ *ni .*/  ni = 7/' "$tmp/good.nml" > "$tmp/bad.nml"
rm -f "$tmp/flow.vtk" "$tmp/cp.csv"
if "$exe" "$tmp/bad.nml" > "$tmp/bad.log" 2>&1; then echo "FAIL: bad input exited 0"; exit 1; fi
[ ! -e "$tmp/flow.vtk" ] && [ ! -e "$tmp/cp.csv" ] || { echo "FAIL: bad input wrote files"; exit 1; }
echo "check_main.sh: ok"
