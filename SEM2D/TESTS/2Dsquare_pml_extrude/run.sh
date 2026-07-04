#!/bin/bash
# End-to-end 2D example: build a bare quad mesh, let the mesher add PML by
# extrusion (pml.input + mater.in), then run sem2d.
#
#   MESHER2D=/path/to/mesher2D SEM2D=/path/to/sem2d ./run.sh
#
# Only the mesher step is required to exercise the feature; the sem2d run is
# optional (skipped if $SEM2D is unset / not found).
set -e
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"

PY=${PYTHON:-python3}
MESHER2D=${MESHER2D:-mesher2D}

echo ">> generating bare quad mesh (no PML): [0,500]x[0,300], 10x6"
$PY gen_mesh2d.py mesh_input.h5 500 300 10 6

echo ">> running mesher2D (adds PML from pml.input + mater.in, writes material.input)"
"$MESHER2D" < mesh.input

if [ -n "$SEM2D" ] && command -v "$SEM2D" >/dev/null 2>&1; then
    echo ">> running sem2d"
    "$SEM2D"
else
    echo ">> SEM2D not set; stopping after mesh generation."
    echo "   Inspect mesh4spec.0000.xmf in ParaView; material.input lists the PML materials."
fi
