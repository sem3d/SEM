#!/usr/bin/env bash
# Run every mesher2D test in this directory and report the outputs.
#
# Usage:
#   ./run_all.sh [path/to/mesher2D]
# Default binary: ../../build/MESH2D/mesher2D  (override with $1 or $MESHER2D).
#
# Each test feeds mesh.input to mesher2D on stdin and checks that the per-proc
# mesh4spec.NNNN.h5 files (and, for on-the-fly, material.input) are produced.

set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MESHER2D="${1:-${MESHER2D:-$HERE/../../build/MESH2D/mesher2D}}"

if [ ! -x "$MESHER2D" ]; then
    echo "ERROR: mesher2D not found/executable at: $MESHER2D"
    echo "       build SEM first, or pass the path:  ./run_all.sh /path/to/mesher2D"
    exit 1
fi
echo "Using mesher2D: $MESHER2D"
PY="${PYTHON:-python3}"

run_case () {
    local dir="$1"; local nproc="$2"; local prep="$3"
    echo
    echo "================  $dir  ================"
    cd "$HERE/$dir" || return 1
    rm -f mesh4spec*.h5
    [ -n "$prep" ] && eval "$prep"
    "$MESHER2D" < mesh.input > outputmesh.log 2>&1
    local rc=$?
    echo "  exit=$rc   $(grep -c . outputmesh.log) log lines"
    grep -E "Creating grid|Wrote material|Nodes, .*Quads|ERR" outputmesh.log | sed 's/^/    /'
    local n=$(ls mesh4spec.*.h5 2>/dev/null | wc -l | tr -d ' ')
    echo "  -> $n / $nproc mesh4spec.NNNN.h5 written"
    if [ -f material.input ]; then echo "  -> material.input ($(head -1 material.input) materials)"; fi
    [ "$n" = "$nproc" ] && echo "  PASS" || echo "  CHECK (expected $nproc proc file(s))"
    cd "$HERE"
}

# on-the-fly variations
run_case onthefly_1mat   1 ""
run_case onthefly_2layer 1 ""
run_case onthefly_pml    1 ""
run_case onthefly_mpi    4 ""

# external mesh inputs (generate the input mesh first)
run_case unv_input  1 "$PY $HERE/gen_unv_mesh.py mesh.unv"
run_case hdf5_input 1 "$PY $HERE/gen_hdf5_mesh.py mesh_input.h5"

# imported mesh + PML added by extrusion (pml.input); base grid has no PML
run_case hdf5_pml   1 "$PY $HERE/gen_hdf5_mesh.py mesh_input.h5"

echo
echo "Done. Inspect each dir's outputmesh.log + material.input, and h5dump mesh4spec.0000.h5."
