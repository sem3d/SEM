#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_hexa27_cube.py - Write a structured 27-node (triquadratic) hex cube in the HDF5
format read by the 3D mesher (menu 4), WITHOUT PML (added later by extrusion).

    /Nodes         (Nnode, 3) float64
    /Sem3D/Hexa27  (Nhex, 27) int32   SEM node order (decoded from shape27.F90):
                   8 corners, 12 edge-mids, 6 face-centers, 1 body-center.
    /Sem3D/Mat     (Nhex,)    int32

    Ex: python3 gen_hexa27_cube.py cube.h5 300 6
"""
import sys
import numpy as np
import h5py

out = sys.argv[1] if len(sys.argv) > 1 else "cube27.h5"
L, n = (300.0, 6)
if len(sys.argv) >= 4:
    L, n = float(sys.argv[2]), int(sys.argv[3])
d = L / n

# SEM Hexa27 node signatures (sx,sy,sz) in {0=min,1=mid,2=max}, order matching shape27_func.
SIG27 = [
    (0,0,0),(2,0,0),(2,2,0),(0,2,0),(0,0,2),(2,0,2),(2,2,2),(0,2,2),      # corners
    (1,0,0),(2,1,0),(1,2,0),(0,1,0),(0,0,1),(2,0,1),(2,2,1),(0,2,1),      # edges
    (1,0,2),(2,1,2),(1,2,2),(0,1,2),                                      # edges
    (1,1,0),(1,0,1),(2,1,1),(1,2,1),(0,1,1),(1,1,2),                      # faces
    (1,1,1),                                                             # center
]

# Fine grid at half-spacing (each element spans 2 fine cells per axis).
nf = 2 * n + 1
def fid(a, b, c):
    return a + b * nf + c * nf * nf

nodes = np.empty((nf ** 3, 3), dtype=np.float64)
for c in range(nf):
    for b in range(nf):
        for a in range(nf):
            nodes[fid(a, b, c)] = [a * d / 2, b * d / 2, c * d / 2]

hexa = np.empty((n ** 3, 27), dtype=np.int32)
e = 0
for k in range(n):
    for j in range(n):
        for i in range(n):
            hexa[e] = [fid(2 * i + sx, 2 * j + sy, 2 * k + sz) for (sx, sy, sz) in SIG27]
            e += 1

mat = np.zeros(n ** 3, dtype=np.int32)
with h5py.File(out, "w") as f:
    f.create_dataset("Nodes", data=nodes)
    g = f.create_group("Sem3D")
    g.create_dataset("Hexa27", data=hexa)
    g.create_dataset("Mat", data=mat)
print(f"written {out}: {len(nodes)} nodes, {len(hexa)} Hexa27, cube [0,{L}]^3 grid {n}^3")
