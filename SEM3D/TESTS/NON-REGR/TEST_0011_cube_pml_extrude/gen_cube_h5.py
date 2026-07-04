#!/usr/bin/env python3
"""Generate an HDF5 Hexa8 cube for the 3D mesher (menu option 4), WITHOUT PML.

The PML layers are added afterwards by the mesher itself, from pml.input
(feature: PML by extrusion of boundary faces). Node ordering follows the SEM
convention used by RectMesh::create_linear_element (v0-3 low-z face, v4-7
high-z face), which the mesher's orientation check expects.

Writes the format read by Mesh3D::read_mesh_file:
    /Nodes        (Nnode, 3) float64
    /Sem3D/Hexa8  (Nhex, 8)  int32   0-based node ids
    /Sem3D/Mat    (Nhex,)    int32   material index per hex (0-based)

Usage: python3 gen_cube_h5.py [out.h5] [L n]   (default cube.h5, 300 m, 6^3)
"""
import sys
import numpy as np
import h5py

out = sys.argv[1] if len(sys.argv) > 1 else "cube.h5"
L, n = 300.0, 6
if len(sys.argv) >= 4:
    L, n = float(sys.argv[2]), int(sys.argv[3])
d = L / n


def nid(i, j, k):
    return i + j * (n + 1) + k * (n + 1) * (n + 1)


nodes = np.empty(((n + 1) ** 3, 3), dtype=np.float64)
for k in range(n + 1):
    for j in range(n + 1):
        for i in range(n + 1):
            nodes[nid(i, j, k)] = [i * d, j * d, k * d]

hexa = np.empty((n ** 3, 8), dtype=np.int32)
e = 0
for k in range(n):
    for j in range(n):
        for i in range(n):
            hexa[e] = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                       nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
            e += 1

mat = np.zeros(n ** 3, dtype=np.int32)

with h5py.File(out, "w") as f:
    f.create_dataset("Nodes", data=nodes)
    g = f.create_group("Sem3D")
    g.create_dataset("Hexa8", data=hexa)
    g.create_dataset("Mat", data=mat)

print(f"written {out}: {len(nodes)} nodes, {len(hexa)} hexa, cube [0,{L}]^3 grid {n}^3")
