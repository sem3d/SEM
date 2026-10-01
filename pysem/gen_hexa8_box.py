#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_hexa8_box.py - Write a structured 8-node (trilinear) hex box/cube in the HDF5
format read by the SEM3D mesher (menu option 4), WITHOUT PML (added later by extrusion).

Datasets created:
    /Nodes         (Nnode, 3) float64   Coordinates [x, y, z]
    /Sem3D/Hexa8   (Nhex, 8)  int32     0-based SEM/VTK node connectivity:
                                        4 bottom corners (z min, counter-clockwise)
                                        4 top corners (z max, counter-clockwise)
    /Sem3D/Mat     (Nhex,)    int32     Material ID per hex (0-based)

Usage:
    # 1. Individual dimensions and element counts (positional):
    python3 gen_hexa8_box.py box.h5 300 200 100 6 4 2

    # 2. Or using flags:
    python3 gen_hexa8_box.py -o box.h5 --size 300 200 100 --num 6 4 2

    # 3. Shorthand uniform cube (compatible with gen_hexa27_cube.py):
    python3 gen_hexa8_box.py cube.h5 300 6

    # 4. Generate companion XDMF file for ParaView inspection:
    python3 gen_hexa8_box.py box.h5 300 200 100 6 4 2 --xmf

    # 5. Default run (cube [0, 300]^3 with 6^3 elements -> box8.h5):
    python3 gen_hexa8_box.py
"""

import os
import sys
import argparse
import numpy as np
import h5py

# SEM Hexa8 node signatures (sx, sy, sz) in {0=min, 1=max}.
# Order matches SEM / VTK convention:
# Bottom face (z=min): (0,0,0), (1,0,0), (1,1,0), (0,1,0)
# Top face    (z=max): (0,0,1), (1,0,1), (1,1,1), (0,1,1)
SIG8 = [
    (0, 0, 0),  # 0: x-min, y-min, z-min
    (1, 0, 0),  # 1: x-max, y-min, z-min
    (1, 1, 0),  # 2: x-max, y-max, z-min
    (0, 1, 0),  # 3: x-min, y-max, z-min
    (0, 0, 1),  # 4: x-min, y-min, z-max
    (1, 0, 1),  # 5: x-max, y-min, z-max
    (1, 1, 1),  # 6: x-max, y-max, z-max
    (0, 1, 1),  # 7: x-min, y-max, z-max
]


def write_xmf(xmf_path, h5_filename, num_nodes, num_elements):
    """Write an XDMF XML wrapper to open the mesh directly in ParaView."""
    base_h5 = os.path.basename(h5_filename)
    xmf_content = f"""<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
<Domain>
<Grid Name="mesh" GridType="Uniform">
<Geometry Type="XYZ">
<DataItem Format="HDF" Dimensions="{num_nodes} 3" NumberType="Float" Precision="8">{base_h5}:/Nodes</DataItem>
</Geometry>
<Topology Type="Hexahedron" NumberOfElements="{num_elements}">
<DataItem Format="HDF" Dimensions="{num_elements} 8" NumberType="Int" Precision="4">{base_h5}:/Sem3D/Hexa8</DataItem>
</Topology>
<Attribute Name="Mat" Center="Cell" AttributeType="Scalar" Dimensions="{num_elements}">
<DataItem Format="HDF" Dimensions="{num_elements}" NumberType="Int" Precision="4">{base_h5}:/Sem3D/Mat</DataItem>
</Attribute>
</Grid>
</Domain>
</Xdmf>
"""
    with open(xmf_path, "w") as f:
        f.write(xmf_content)


def generate_hexa8_box(Lx, Ly, Lz, nx, ny, nz, out=None, mat_id=0, create_xmf=False):
    """
    Generate a structured grid of Hexa8 (8-node linear) elements.

    Parameters:
        Lx, Ly, Lz : float
            Box dimensions along X, Y, Z (meters).
        nx, ny, nz : int
            Number of elements along X, Y, Z.
        out : str, optional
            Path to output HDF5 file. If provided, saves the mesh to disk.
        mat_id : int, optional
            Material ID for the elements (default: 0).
        create_xmf : bool, optional
            Whether to write a companion .xmf file for ParaView (default: False).

    Returns:
        nodes : ndarray of shape ((nx+1)*(ny+1)*(nz+1), 3), float64
        hexa  : ndarray of shape (nx*ny*nz, 8), int32
        mat   : ndarray of shape (nx*ny*nz,), int32
    """
    nx, ny, nz = int(nx), int(ny), int(nz)
    Lx, Ly, Lz = float(Lx), float(Ly), float(Lz)

    if nx <= 0 or ny <= 0 or nz <= 0:
        raise ValueError(f"nx, ny, nz must be positive integers (got {nx}, {ny}, {nz})")
    if Lx <= 0 or Ly <= 0 or Lz <= 0:
        raise ValueError(f"Lx, Ly, Lz must be positive numbers (got {Lx}, {Ly}, {Lz})")

    dx = Lx / nx
    dy = Ly / ny
    dz = Lz / nz

    # 1D node indexing helper
    def nid(i, j, k):
        return i + j * (nx + 1) + k * (nx + 1) * (ny + 1)

    num_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    nodes = np.empty((num_nodes, 3), dtype=np.float64)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nodes[nid(i, j, k)] = [i * dx, j * dy, k * dz]

    num_elements = nx * ny * nz
    hexa = np.empty((num_elements, 8), dtype=np.int32)
    e = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                hexa[e] = [nid(i + sx, j + sy, k + sz) for (sx, sy, sz) in SIG8]
                e += 1

    mat = np.full(num_elements, mat_id, dtype=np.int32)

    if out:
        with h5py.File(out, "w") as f:
            f.create_dataset("Nodes", data=nodes)
            g = f.create_group("Sem3D")
            g.create_dataset("Hexa8", data=hexa)
            g.create_dataset("Mat", data=mat)

        msg = (
            f"written {out}: {len(nodes)} nodes, {len(hexa)} Hexa8, "
            f"box [0,{Lx}]x[0,{Ly}]x[0,{Lz}], grid {nx}x{ny}x{nz} "
            f"(dx={dx:.4g}, dy={dy:.4g}, dz={dz:.4g})"
        )

        if create_xmf:
            xmf_path = f"{out}.xmf"
            write_xmf(xmf_path, out, len(nodes), len(hexa))
            msg += f" [+ {xmf_path}]"

        print(msg)

    return nodes, hexa, mat


def parse_arguments():
    """Flexible command line parser supporting both positional arguments and flags."""
    raw_args = sys.argv[1:]

    # Separate any trailing flags like --xmf if passed with positional arguments
    create_xmf = "--xmf" in raw_args
    filtered_args = [a for a in raw_args if a != "--xmf"]

    has_other_flags = any(arg.startswith("-") for arg in filtered_args)

    if not has_other_flags and len(filtered_args) > 0:
        out = filtered_args[0]
        # Defaults
        Lx, Ly, Lz = 300.0, 300.0, 300.0
        nx, ny, nz = 6, 6, 6
        mat_id = 0

        # Uniform cube shorthand: out L n [mat]
        if len(filtered_args) in (3, 4):
            L = float(filtered_args[1])
            n = int(filtered_args[2])
            Lx, Ly, Lz = L, L, L
            nx, ny, nz = n, n, n
            if len(filtered_args) == 4:
                mat_id = int(filtered_args[3])
        # Full box specification: out Lx Ly Lz nx ny nz [mat]
        elif len(filtered_args) >= 7:
            Lx = float(filtered_args[1])
            Ly = float(filtered_args[2])
            Lz = float(filtered_args[3])
            nx = int(filtered_args[4])
            ny = int(filtered_args[5])
            nz = int(filtered_args[6])
            if len(filtered_args) >= 8:
                mat_id = int(filtered_args[7])

        return out, Lx, Ly, Lz, nx, ny, nz, mat_id, create_xmf

    # Argparse mode for standard flags and --help
    parser = argparse.ArgumentParser(
        description="Write a structured 8-node linear hexahedron (Hexa8) box in HDF5 format for SEM3D mesher.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("out", nargs="?", default="box8.h5", help="Output HDF5 file name")
    parser.add_argument(
        "--size",
        nargs=3,
        type=float,
        metavar=("LX", "LY", "LZ"),
        default=None,
        help="Box dimensions in meters (Lx Ly Lz)",
    )
    parser.add_argument(
        "--num",
        nargs=3,
        type=int,
        metavar=("NX", "NY", "NZ"),
        default=None,
        help="Number of elements along (nx ny nz)",
    )
    parser.add_argument("--lx", type=float, default=300.0, help="Dimension along X (m)")
    parser.add_argument("--ly", type=float, default=300.0, help="Dimension along Y (m)")
    parser.add_argument("--lz", type=float, default=300.0, help="Dimension along Z (m)")
    parser.add_argument("--nx", type=int, default=6, help="Number of elements along X")
    parser.add_argument("--ny", type=int, default=6, help="Number of elements along Y")
    parser.add_argument("--nz", type=int, default=6, help="Number of elements along Z")
    parser.add_argument("--mat", type=int, default=0, help="Material index for all elements")
    parser.add_argument("--xmf", action="store_true", help="Generate companion .xmf file for ParaView")

    args = parser.parse_args()

    Lx = args.size[0] if args.size is not None else args.lx
    Ly = args.size[1] if args.size is not None else args.ly
    Lz = args.size[2] if args.size is not None else args.lz

    nx = args.num[0] if args.num is not None else args.nx
    ny = args.num[1] if args.num is not None else args.ny
    nz = args.num[2] if args.num is not None else args.nz

    return args.out, Lx, Ly, Lz, nx, ny, nz, args.mat, args.xmf


def main():
    out, Lx, Ly, Lz, nx, ny, nz, mat_id, create_xmf = parse_arguments()
    generate_hexa8_box(Lx, Ly, Lz, nx, ny, nz, out=out, mat_id=mat_id, create_xmf=create_xmf)


if __name__ == "__main__":
    main()
