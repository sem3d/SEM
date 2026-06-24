#!/usr/bin/env python3
"""Generate the 2D quad mesh (HDF5) that mesh2dc reads, for the aniso2d tests.

Writes the format expected by MESH2D/main_mesh2d.cpp::read_sem_mesh:
    /Nodes        (Nnode, 2) float64   -> column 0 = x, column 1 = z
    /Sem2D/Quad4  (Nquad, 4) int32     -> 0-based node ids, CCW order [BL,BR,TR,TL]
    /Sem2D/Mat    (Nquad,)   int32     -> material index per quad (0-based)

A regular nx*nz grid over [0,Lx] x [0,Lz], ONE material (index 0). The same mesh
serves every test (T1 solid, T2/T3 fluid): the solid/fluid choice is in
material.input (material 0 = "S" or "F"), not in the mesh.

Then partition it for SEM2D with the built mesher:
    mesh2dc  Nproc  mesh_input.h5  mesh4spec      ->  mesh4spec.0000.h5 ...

Usage:
    python3 gen_test_mesh2d.py                       # default [0,0.05]x[0,0.03], 100x60
    python3 gen_test_mesh2d.py Lx Lz nx nz [out.h5]
"""
import sys
import numpy as np
import h5py

# ---- parameters (keep the box consistent with gen_test_cstar2d.py) ---------
Lx, Lz = 0.05, 0.03
nx, nz = 100, 60
out = "mesh_input.h5"
if len(sys.argv) >= 5:
    Lx, Lz = float(sys.argv[1]), float(sys.argv[2])
    nx, nz = int(sys.argv[3]), int(sys.argv[4])
if len(sys.argv) >= 6:
    out = sys.argv[5]

dx, dz = Lx / nx, Lz / nz
nnode = (nx + 1) * (nz + 1)
nquad = nx * nz


def nid(ix, iz):
    """0-based node id on the (nx+1) x (nz+1) grid."""
    return ix + iz * (nx + 1)


# nodes
nodes = np.empty((nnode, 2), dtype=np.float64)
for iz in range(nz + 1):
    for ix in range(nx + 1):
        nodes[nid(ix, iz), 0] = ix * dx
        nodes[nid(ix, iz), 1] = iz * dz

# quads (CCW: bottom-left, bottom-right, top-right, top-left)
quads = np.empty((nquad, 4), dtype=np.int32)
q = 0
for ez in range(nz):
    for ex in range(nx):
        quads[q, 0] = nid(ex,     ez)      # BL
        quads[q, 1] = nid(ex + 1, ez)      # BR
        quads[q, 2] = nid(ex + 1, ez + 1)  # TR
        quads[q, 3] = nid(ex,     ez + 1)  # TL
        q += 1

# single material -> index 0 everywhere
mat = np.zeros(nquad, dtype=np.int32)

with h5py.File(out, "w") as f:
    f.create_dataset("Nodes", data=nodes)
    g = f.create_group("Sem2D")
    g.create_dataset("Quad4", data=quads)
    g.create_dataset("Mat", data=mat)

print(f"written {out}: {nnode} nodes, {nquad} quads, domain [0,{Lx}]x[0,{Lz}], grid {nx}x{nz}")
print(f"partition it:  mesh2dc <Nproc> {out} mesh4spec")
