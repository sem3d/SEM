#!/usr/bin/env python3
"""Generate an Ideas .unv quad mesh for mesher2D menu option 3 (Ideas .unv files).

Writes the three blocks the SEM .unv reader (COMMON/read_unv.cpp) consumes:
  2411 Nodes    : per node  ->  "label 0 0 0"  then  "x y z"   (z=0; x,z carry the 2D coords)
  2412 Elements : per quad  ->  "label 44 0 0 0 4"  then 4 node ids (1-based, CCW)
  2467 Groups   : one element group (-> material index 0). NOTE this reader's 2467 header is
                  6 dummies + n_entities + name, and element entities are "8 elemtag 0 0".

All quads go to a single group, so every element gets material index 0 (the reader sets
m_mat = group index). Without a group, m_mat would be -1 and Mesh2D::read_mesh would abort.

Usage:
    python3 gen_unv_mesh.py [out.unv] [Lx Lz nx nz]
    (default mesh.unv, [0,0.05]x[0,0.03], 10x6)
"""
import sys

out = sys.argv[1] if len(sys.argv) >= 2 else "mesh.unv"
Lx, Lz, nx, nz = 0.05, 0.03, 10, 6
if len(sys.argv) >= 6:
    Lx, Lz = float(sys.argv[2]), float(sys.argv[3])
    nx, nz = int(sys.argv[4]), int(sys.argv[5])

dx, dz = Lx / nx, Lz / nz


def nid1(ix, iz):
    """1-based node id on the (nx+1)x(nz+1) grid."""
    return ix + iz * (nx + 1) + 1


lines = []

# ---- 2411 : nodes ----------------------------------------------------------
lines += ["-1", "2411"]
for iz in range(nz + 1):
    for ix in range(nx + 1):
        lab = nid1(ix, iz)
        lines.append(f"{lab:>10d}         0         0         0")
        lines.append(f"{ix*dx:25.16E}{iz*dz:25.16E}{0.0:25.16E}")
lines.append("-1")

# ---- 2412 : elements (Quad4 = type 44) -------------------------------------
lines += ["-1", "2412"]
elab = 0
for ez in range(nz):
    for ex in range(nx):
        elab += 1
        conn = [nid1(ex, ez), nid1(ex + 1, ez), nid1(ex + 1, ez + 1), nid1(ex, ez + 1)]
        lines.append(f"{elab:>10d}        44         0         0         0         4")
        lines.append("".join(f"{c:>10d}" for c in conn))
lines.append("-1")

# ---- 2467 : one element group -> material 0 --------------------------------
nelem = nx * nz
lines += ["-1", "2467"]
lines.append(f"{1:>10d}         0         0         0         0         0         0{nelem:>10d}")
lines.append("mat0")
ent = []
for e in range(1, nelem + 1):
    ent.append(f"         8{e:>10d}         0         0")
lines += ent
lines.append("-1")

with open(out, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"written {out}: {(nx+1)*(nz+1)} nodes, {nelem} quads (type 44), 1 group -> mat 0")
