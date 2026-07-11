#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mkcase.py - Generate a self-contained SEM2D case directory for the MPI np=1 vs np=4
comparison suite. Writes the mesh (external HDF5 Quad4/Quad8 or on-the-fly mat.dat), the
materials (mater.in [+ material.spec + Cstar.h5 for anisotropic-from-file], the input.spec
(with/without PML), capteurs.dat and the mesher stdin (mesh.input, parametrized by NPROCS).

Domain: [0,500] x [0,300], base grid 10 x 6 (dx=dz=50). Coupling cases: solid top half
(z>150), fluid bottom half (z<150). Sources/receivers sit at element MIDPOINTS (an
on-boundary source/receiver is injected/read inconsistently across MPI ranks -- see
plan 2026-07-09_sf-coupling-2d-solution-plan.md -- which is a harness pitfall, not a bug).

    python3 mkcase.py --dir DIR --phys solid|fluid|fluid_aniso|sf_iso|sf_aniso \\
                      --mesh onthefly|quad4|quad8 --pml 0|1 --nprocs N
"""
import argparse
import os

import numpy as np
import h5py

LX, LZ, NX, NZ = 500.0, 300.0, 10, 6
DX, DZ = LX / NX, LZ / NZ
ZIFACE = 150.0  # solid/fluid interface for coupling cases

# base isotropic materials
SOLID = ("S", 2000.0, 1150.0, 2200.0)   # type, Vp, Vs, Rho
FLUID = ("F", 1500.0, 0.0, 1000.0)


def is_coupling(phys):
    return phys in ("sf_iso", "sf_aniso")


def is_aniso(phys):
    return phys in ("solid_aniso", "fluid_aniso", "sf_aniso")


# ---------------------------------------------------------------- external meshes
def _elem_mat(cx, cz, phys):
    """Material index of an element with centroid (cx,cz)."""
    if is_coupling(phys):
        return 0 if cz > ZIFACE else 1   # 0=solid(top), 1=fluid(bottom)
    return 0


def write_quad4(path, phys):
    def nid(ix, iz):
        return ix + iz * (NX + 1)
    nodes = np.array([[ix * DX, iz * DZ] for iz in range(NZ + 1) for ix in range(NX + 1)],
                     dtype=np.float64)
    quads, mats = [], []
    for ez in range(NZ):
        for ex in range(NX):
            quads.append([nid(ex, ez), nid(ex + 1, ez), nid(ex + 1, ez + 1), nid(ex, ez + 1)])
            mats.append(_elem_mat((ex + 0.5) * DX, (ez + 0.5) * DZ, phys))
    with h5py.File(path, "w") as f:
        f.create_dataset("Nodes", data=nodes)
        g = f.create_group("Sem2D")
        g.create_dataset("Quad4", data=np.array(quads, dtype=np.int32))
        g.create_dataset("Mat", data=np.array(mats, dtype=np.int32))


def write_quad8(path, phys):
    coords, cmap = [], {}
    def node(x, z):
        k = (round(x / (min(DX, DZ) * 1e-6)), round(z / (min(DX, DZ) * 1e-6)))
        if k not in cmap:
            cmap[k] = len(coords); coords.append((x, z))
        return cmap[k]
    quads, mats = [], []
    for ez in range(NZ):
        for ex in range(NX):
            x0, z0, x1, z1 = ex * DX, ez * DZ, (ex + 1) * DX, (ez + 1) * DZ
            xm, zm = 0.5 * (x0 + x1), 0.5 * (z0 + z1)
            c = [node(x0, z0), node(x1, z0), node(x1, z1), node(x0, z1),
                 node(xm, z0), node(x1, zm), node(xm, z1), node(x0, zm)]
            quads.append(c)
            mats.append(_elem_mat(xm, zm, phys))
    with h5py.File(path, "w") as f:
        f.create_dataset("Nodes", data=np.array(coords, dtype=np.float64))
        g = f.create_group("Sem2D")
        g.create_dataset("Quad8", data=np.array(quads, dtype=np.int32))
        g.create_dataset("Mat", data=np.array(mats, dtype=np.int32))


# ---------------------------------------------------------------- on-the-fly mat.dat
def write_matdat(path, phys, pml):
    npml = 4 if pml else 0
    lines = ["0.0", "%g" % LX, "%g" % DX, "%g" % LZ]
    if is_coupling(phys):
        lines += ["2", "%g %d" % (LZ - ZIFACE, NZ // 2), "%g %d" % (ZIFACE, NZ // 2)]
    else:
        lines += ["1", "%g %d" % (LZ, NZ)]
    lines += ["%d" % npml, "0 1", "5 2 0.001 0. 0."]  # Rc=0.001 -> positive Apow
    lines += ["4"]
    open(path, "w").write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- materials
def write_materin(path, phys):
    mats = []
    if is_coupling(phys):
        mats = [SOLID, FLUID]
    elif phys in ("fluid", "fluid_aniso"):
        mats = [FLUID]
    else:
        mats = [SOLID]
    with open(path, "w") as f:
        f.write("%d\n" % len(mats))
        for t, vp, vs, rho in mats:
            f.write("%s %g %g %g 0. 0.\n" % (t, vp, vs, rho))


def _cstar_group(f, name, val):
    g = f.create_group(name)
    g.create_dataset("samples", data=np.full((2, 2), val, dtype=np.float64))
    g.attrs["xMinGlob"] = np.array([-LX, -LZ], dtype=np.float64)   # cover PML region too
    g.attrs["xMaxGlob"] = np.array([2 * LX, 2 * LZ], dtype=np.float64)


def write_cstar_fluid(path):
    vp, rho = FLUID[1], FLUID[3]
    kappa = rho * vp * vp
    with h5py.File(path, "w") as f:
        _cstar_group(f, "iRho11", 1.0 / rho)
        _cstar_group(f, "iRho12", 0.0)
        _cstar_group(f, "iRho22", 1.0 / rho)
        _cstar_group(f, "iKappa", 1.0 / kappa)


def write_cstar_solid(path):
    vp, vs, rho = SOLID[1], SOLID[2], SOLID[3]
    mu = rho * vs * vs
    lam = rho * (vp * vp - 2 * vs * vs)
    C11 = lam + 2 * mu
    with h5py.File(path, "w") as f:
        for k, v in [("C11", C11), ("C22", C11), ("C33", mu),
                     ("C12", lam), ("C13", 0.0), ("C23", 0.0), ("Rho", rho)]:
            _cstar_group(f, k, v)


def write_materialspec(d, phys):
    """material.spec + Cstar files for anisotropic-from-file cases."""
    if phys == "solid_aniso":
        write_cstar_solid(os.path.join(d, "Cstar_solid.h5"))
        open(os.path.join(d, "material.spec"), "w").write(
            'material 0 {\n  domain = solid;\n  deftype = Hooke_Aniso;\n'
            '  spacedef = file;\n  filename0 = "Cstar_solid.h5";\n};\n')
    elif phys == "fluid_aniso":
        write_cstar_fluid(os.path.join(d, "Cstar_fluid.h5"))
        open(os.path.join(d, "material.spec"), "w").write(
            'material 0 {\n  domain = fluid;\n  deftype = Fluid_Aniso;\n'
            '  spacedef = file;\n  filename0 = "Cstar_fluid.h5";\n};\n')
    elif phys == "sf_aniso":
        write_cstar_solid(os.path.join(d, "Cstar_solid.h5"))
        write_cstar_fluid(os.path.join(d, "Cstar_fluid.h5"))
        open(os.path.join(d, "material.spec"), "w").write(
            'material 0 {\n  domain = solid;\n  deftype = Hooke_Aniso;\n'
            '  spacedef = file;\n  filename0 = "Cstar_solid.h5";\n};\n'
            'material 1 {\n  domain = fluid;\n  deftype = Fluid_Aniso;\n'
            '  spacedef = file;\n  filename0 = "Cstar_fluid.h5";\n};\n')


# ---------------------------------------------------------------- input.spec, capteurs
def write_inputspec(path, phys, pml):
    # source in the solid region for coupling, else domain centre; both at midpoints
    sx, sz = 275.0, (225.0 if is_coupling(phys) else 175.0)
    stype = "fluidpulse" if phys in ("fluid", "fluid_aniso") else "impulse"
    pml_block = "pml_infos {\n    pml_type = PML;\n};\n\n" if pml else ""
    open(path, "w").write(
        '# -*- mode: perl -*-\nrun_name = "mpi_compare";\nsim_time = 1.0;\n'
        'mesh_file = "mesh4spec";\nmat_file = "material.input";\ndim=2;\nngll=5;\n\n'
        + pml_block +
        'snapshots {\n    save_snap = false;\n    snap_interval = 0.05;\n};\n\n'
        'save_traces = true;\nstation_file = "capteurs.dat";\ntraces_format=hdf5;\n\n'
        'source {\n    coords = %g %g;\n    type = %s;\n    dir = 1. 0.;\n'
        '    func = ricker;\n    tau = .3;\n    freq = 4.;\n};\n\n' % (sx, sz, stype) +
        'time_scheme {\n    accel_scheme = false;\n    veloc_scheme = true;\n'
        '    alpha = 0.5;\n    beta = 0.5;\n    gamma = 1;\n    courant = 0.2;\n};\n')


def write_capteurs(path, phys):
    if is_coupling(phys):
        recs = [(275.0, 225.0), (275.0, 75.0), (175.0, 225.0)]   # solid, fluid, solid
    else:
        recs = [(275.0, 175.0), (175.0, 125.0)]
    open(path, "w").write("".join("%g %g\n" % r for r in recs))


def write_meshinput(path, mesh, nprocs):
    if mesh == "onthefly":
        open(path, "w").write("%d\n1\n" % nprocs)              # NPROCS, choice 1
    else:
        open(path, "w").write("%d\n4\nmesh_input.h5\n" % nprocs)  # NPROCS, choice 4, file


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", required=True)
    p.add_argument("--phys", required=True,
                   choices=["solid", "solid_aniso", "fluid", "fluid_aniso", "sf_iso", "sf_aniso"])
    p.add_argument("--mesh", required=True, choices=["onthefly", "quad4", "quad8"])
    p.add_argument("--pml", type=int, default=0)
    p.add_argument("--nprocs", type=int, default=1)
    a = p.parse_args()
    d = a.dir
    os.makedirs(d, exist_ok=True)

    write_materin(os.path.join(d, "mater.in"), a.phys)
    if is_aniso(a.phys):
        write_materialspec(d, a.phys)
    write_inputspec(os.path.join(d, "input.spec"), a.phys, a.pml)
    write_capteurs(os.path.join(d, "capteurs.dat"), a.phys)
    write_meshinput(os.path.join(d, "mesh.input"), a.mesh, a.nprocs)

    if a.mesh == "onthefly":
        write_matdat(os.path.join(d, "mat.dat"), a.phys, a.pml)
    else:
        if a.mesh == "quad4":
            write_quad4(os.path.join(d, "mesh_input.h5"), a.phys)
        else:
            write_quad8(os.path.join(d, "mesh_input.h5"), a.phys)
        if a.pml:
            open(os.path.join(d, "pml.input"), "w").write("x- 4\nx+ 4\nz- 4\n")
    print("case ready: %s (phys=%s mesh=%s pml=%d nprocs=%d)"
          % (d, a.phys, a.mesh, a.pml, a.nprocs))


if __name__ == "__main__":
    main()
