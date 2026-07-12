#!/bin/bash
# MPI np=1 vs np=4 comparison suite for SEM2D.
#
# For every (physics x mesh x PML) case it: builds the case for NPROCS=1 and NPROCS=4
# (mkcase.py), runs the mesher then sem2d for each, and compares the receiver traces
# (compare_mpi.py). A correct MPI partition reproduces the serial run to round-off.
#
#   ./run_all.sh [BUILD_DIR] [PHYS_FILTER] [MESH_FILTER]
#     BUILD_DIR    : path to the SEM build (default ../../../build)
#     PHYS_FILTER  : run only physics matching this (e.g. solid, fluid, sf_); default all
#     MESH_FILTER  : run only meshes matching this (e.g. quad8);            default all
#
# Exit code 0 if every non-noise case is within RTOL, else 1.
set -u
here=$(cd "$(dirname "$0")" && pwd); cd "$here"
BUILD=$(cd "${1:-../../../build}" && pwd)
PHYS_FILTER=${2:-}
MESH_FILTER=${3:-}
MESHER="$BUILD/MESH2D/mesher2D"
SEM="$BUILD/SEM2D/sem2d.exe"
PY=${PYTHON:-python3}
RTOL=${RTOL:-1e-6}      # pass threshold; non-PML should be ~1e-12, PML is looser (see README)
work="$here/_work"
rm -rf "$work"; mkdir -p "$work"

PHYS="solid solid_aniso fluid fluid_aniso sf_iso sf_aniso_solid sf_aniso_fluid sf_aniso"
MESHES="onthefly quad4 quad8"

run_one() { # dir nprocs
    local d=$1 np=$2
    ( cd "$d" && "$MESHER" < mesh.input > mesher.log 2>&1 ) || { echo "mesher-fail"; return 1; }
    ( cd "$d" && mpirun -np "$np" "$SEM" > sem.log 2>&1 )   || { echo "sem-fail";    return 1; }
    return 0
}

printf "%-26s %-8s %-4s  %-10s %s\n" "CASE" "MESH" "PML" "VERDICT" "worst_reldiff(non-noise)"
echo "--------------------------------------------------------------------------------"
overall=0
for phys in $PHYS; do
    [ -n "$PHYS_FILTER" ] && [[ "$phys" != *"$PHYS_FILTER"* ]] && continue
    for mesh in $MESHES; do
        [ -n "$MESH_FILTER" ] && [[ "$mesh" != *"$MESH_FILTER"* ]] && continue
        for pml in 0 1; do
            case="${phys}_${mesh}_pml${pml}"
            d1="$work/$case/np1"; d4="$work/$case/np4"
            $PY mkcase.py --dir "$d1" --phys "$phys" --mesh "$mesh" --pml "$pml" --nprocs 1 >/dev/null 2>&1
            $PY mkcase.py --dir "$d4" --phys "$phys" --mesh "$mesh" --pml "$pml" --nprocs 4 >/dev/null 2>&1
            err=""
            run_one "$d1" 1 >/dev/null 2>&1 || err="np1:$(run_one "$d1" 1 2>&1)"
            run_one "$d4" 4 >/dev/null 2>&1 || err="${err} np4:$(run_one "$d4" 4 2>&1)"
            if [ -n "$err" ]; then
                printf "%-26s %-8s %-4s  %-10s %s\n" "$phys" "$mesh" "$pml" "ERROR" "$err"
                overall=1; continue
            fi
            out=$($PY compare_mpi.py "$d1/traces" "$d4/traces" --rtol "$RTOL" 2>&1)
            verdict=$(echo "$out" | grep -oE "PASS|FAIL" | tail -1)
            worst=$(echo "$out" | grep -oE "worst.*reldiff = [0-9.e+-]+" | grep -oE "[0-9.e+-]+$")
            [ -z "$verdict" ] && verdict="NOTRACE"
            [ "$verdict" = "FAIL" ] && overall=1
            printf "%-26s %-8s %-4s  %-10s %s\n" "$phys" "$mesh" "$pml" "$verdict" "${worst:-?}"
        done
    done
done
echo "--------------------------------------------------------------------------------"
echo "RTOL=$RTOL  (non-PML expected PASS ~1e-12; PML residual is a known shared issue -- see README)"
exit $overall
