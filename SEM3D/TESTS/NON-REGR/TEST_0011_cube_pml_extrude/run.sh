#!/bin/bash
# End-to-end 3D example: build a bare cube mesh, let the mesher add PML by
# extrusion (pml.input), then run sem3d.
#
#   MESHER=/path/to/mesher SEM3D=/path/to/sem3d ./run.sh
#
# Only the mesher step is required to exercise the PML-extrusion feature; the
# sem3d run is optional (skipped if $SEM3D is not set / not found).
set -e
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"

PY=${PYTHON:-python3}
MESHER=${MESHER:-mesher}

echo ">> generating bare cube mesh (no PML)"
$PY gen_cube_h5.py cube.h5 300 6

# The mesher rewrites material.input and *appends* to material.spec, so start
# each run from the pristine base files (keeps re-runs idempotent).
printf '1\nS 6300. 2500. 2800. 0. 0.\n' > material.input
cat > material.spec <<'EOF'
material 0 {
    domain = solid;
    deftype = Vp_Vs_Rho;
    rho = 2800.;
    vp = 6300.;
    vs = 2500.;
};
EOF

echo ">> running mesher (adds PML layers from pml.input, rewrites material.input/.spec)"
"$MESHER" < mesh.input

if [ -n "$SEM3D" ] && command -v "$SEM3D" >/dev/null 2>&1; then
    echo ">> running sem3d"
    "$SEM3D"
else
    echo ">> SEM3D not set; stopping after mesh generation."
    echo "   Inspect mesh4spec.0000.xmf in ParaView; material.input now lists the PML materials."
fi
