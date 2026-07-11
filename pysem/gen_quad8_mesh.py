#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_quad8_mesh.py - Write a structured 8-node (serendipity) quad mesh in the HDF5
format read by mesher2D menu option 4 (Mesh2D::read_sem_mesh).

    /Nodes        (Nnode, 2) float64   x, z
    /Sem2D/Quad8  (Nquad, 8) int32     0-based: 4 corners CCW [BL,BR,TR,TL] then
                                       4 edge-mids [bottom, right, top, left]
                                       (matches SEM2D shape8.F90 node convention)
    /Sem2D/Mat    (Nquad,)   int32     material index per quad

    Ex: python3 gen_quad8_mesh.py mesh_input.h5 500 300 10 6
"""
import sys
import numpy as np
import h5py

out = sys.argv[1] if len(sys.argv) >= 2 else "mesh_input.h5"
Lx, Lz, nx, nz = 0.05, 0.03, 10, 6
if len(sys.argv) >= 6:
    Lx, Lz, nx, nz = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
dx, dz = Lx / nx, Lz / nz

# Unique node dedup by coordinate (corners + edge midpoints shared between elements).
coords, cmap = [], {}
def node(x, z):
    k = (round(x / (min(dx, dz) * 1e-6)), round(z / (min(dx, dz) * 1e-6)))
    if k in cmap:
        return cmap[k]
    cmap[k] = len(coords); coords.append((x, z)); return cmap[k]

quads = []
for ez in range(nz):
    for ex in range(nx):
        x0, z0 = ex * dx, ez * dz
        x1, z1 = x0 + dx, z0 + dz
        xm, zm = 0.5 * (x0 + x1), 0.5 * (z0 + z1)
        c0, c1, c2, c3 = node(x0, z0), node(x1, z0), node(x1, z1), node(x0, z1)
        mb, mr, mt, ml = node(xm, z0), node(x1, zm), node(xm, z1), node(x0, zm)
        quads.append([c0, c1, c2, c3, mb, mr, mt, ml])

nodes = np.array(coords, dtype=np.float64)
quads = np.array(quads, dtype=np.int32)
mat = np.zeros(len(quads), dtype=np.int32)
with h5py.File(out, "w") as f:
    f.create_dataset("Nodes", data=nodes)
    g = f.create_group("Sem2D")
    g.create_dataset("Quad8", data=quads)
    g.create_dataset("Mat", data=mat)
print(f"written {out}: {len(nodes)} nodes, {len(quads)} Quad8, [0,{Lx}]x[0,{Lz}], grid {nx}x{nz}")
