#!/bin/bash
# SEM3D PML-heterogeneous-inherit test: a solid cube whose material is a heterogeneous FILE
# field (Kappa_Mu_Rho spacedef=file, Kappa linear in x) with PML on the boundary. The PML
# materials are declared `copy=0; domain=solidpml`, so they must inherit the base FILE
# definition and read the SAME per-GLL field (plan 2026-07-10_pml-heterogeneous-inherit).
# Reuses the NON-REGR/TEST_0006 geometry (26 PML materials) but swaps the random field for a
# controlled gradient so the result is checkable. Runs a short sim with the material snapshot
# (item 8) on and validates via check.py.
#
#   ./run_all.sh [BUILD_DIR]
set -u
here=$(cd "$(dirname "$0")" && pwd); cd "$here"
BUILD=$(cd "${1:-../../../build}" && pwd)
MESHER="$BUILD/MESH/mesher"; SEM="$BUILD/SEM3D/sem3d.exe"
PY=${PYTHON:-python3}
SRC="$here/../NON-REGR/TEST_0006_rand_cube_pulse_ricker"
work="$here/_work"; rm -rf "$work"; mkdir -p "$work/mat/h5"

cp "$SRC"/{mat.dat,mesh.input,material.spec,mater.in,stations.txt,input.spec} "$work/"
$PY gen_field.py "$work/mat/h5" >/dev/null          # controlled gradient (Kappa linear in x)
# short sim, keep the material snapshot (geometry*.h5); np=1
sed -i.bak -E 's/sim_time *= *[0-9.]+;/sim_time = 0.02;/; s/snap_interval *= *[0-9.]+;/snap_interval = 0.02;/' "$work/input.spec"
printf "1\n1\n" > "$work/mesh.input"

( cd "$work" && "$MESHER" < mesh.input > mesher.log 2>&1 ) || { echo "mesher-fail"; exit 1; }
( cd "$work" && mkdir -p sem && mv mesh4spec* sem/ 2>/dev/null && mpirun -np 1 "$SEM" > sem.log 2>&1 ) \
    || { echo "sem-fail (see $work/sem.log)"; exit 1; }

echo "domain: $(grep -iE 'solid .*elem.*ngll|pml .*elem.*ngll' "$work/sem.log" | tr '\n' ' ')"
$PY check.py "$work"
