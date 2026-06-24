#!/usr/bin/env python3
"""Generate an HDF5 quad mesh for mesher2D menu option 4 (HDF5 Quad files).

Writes the format read by Mesh2D::read_sem_mesh (same as mesh2dc input):
    /Nodes        (Nnode, 2) float64   -> col 0 = x, col 1 = z
    /Sem2D/Quad4  (Nquad, 4) int32     -> 0-based node ids, CCW [BL,BR,TR,TL]
    /Sem2D/Mat    (Nquad,)   int32     -> material index per quad (0-based)

Usage:
    python3 gen_hdf5_mesh.py [out.h5] [Lx Lz nx nz]
    (default mesh_input.h5, [0,0.05]x[0,0.03], 10x6, one material 0)
"""
import sys
import numpy as np
import h5py

out = sys.argv[1] if len(sys.argv) >= 2 else "mesh_input.h5"
Lx, Lz, nx, nz = 0.05, 0.03, 10, 6
if len(sys.argv) >= 6:
    Lx, Lz = float(sys.argv[2]), float(sys.argv[3])
    nx, nz = int(sys.argv[4]), int(sys.argv[5])

dx, dz = Lx / nx, Lz / nz
nnode = (nx + 1) * (nz + 1)
nquad = nx * nz


def nid(ix, iz):
    return ix + iz * (nx + 1)


nodes = np.empty((nnode, 2), dtype=np.float64)
for iz in range(nz + 1):
    for ix in range(nx + 1):
        nodes[nid(ix, iz), 0] = ix * dx
        nodes[nid(ix, iz), 1] = iz * dz

quads = np.empty((nquad, 4), dtype=np.int32)
q = 0
for ez in range(nz):
    for ex in range(nx):
        quads[q] = [nid(ex, ez), nid(ex + 1, ez), nid(ex + 1, ez + 1), nid(ex, ez + 1)]
        q += 1

mat = np.zeros(nquad, dtype=np.int32)

with h5py.File(out, "w") as f:
    f.create_dataset("Nodes", data=nodes)
    g = f.create_group("Sem2D")
    g.create_dataset("Quad4", data=quads)
    g.create_dataset("Mat", data=mat)

print(f"written {out}: {nnode} nodes, {nquad} quads, [0,{Lx}]x[0,{Lz}], grid {nx}x{nz}")
