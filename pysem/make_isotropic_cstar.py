# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
make_isotropic_cstar.py - Write a SEM3D anisotropic-material HDF5
(MATDEF_HOOKE_ANISO / CSTAR) for a spatially UNIFORM isotropic material,
given Vp, Vs, rho directly -- no homo3dfft run needed.

Used to test SEM3D's Cstar.h5 reading path (build_prop_files.F90) in
isolation, by comparing against a reference run of the same material read
via the simple mater.in/material.input path.

Same group/dataset layout as cstar2h5.py's write_h5 (one group per property,
"samples" dataset + xMinGlob/xMaxGlob attributes) -- reuses ELASTIC_NAMES
from cstar2h5.py so the two stay in sync.

Usage:
    python make_isotropic_cstar.py out.h5 --vp 6300 --vs 2500 --rho 2800 \
        --xmin -500 --xmax 500 --ymin -500 --ymax 500 --zmin -500 --zmax 500
"""
import argparse

import h5py
import numpy as np

from cstar2h5 import ELASTIC_NAMES


def isotropic_cij(vp, vs, rho):
    """Kelvin/Voigt Cij (Pa) + density (kg/m^3) for an isotropic material."""
    c11 = rho * vp ** 2
    c44 = 2.0 * rho * vs ** 2
    c12 = c11 - c44
    values = dict.fromkeys(ELASTIC_NAMES, 0.0)
    values["C11"] = values["C22"] = values["C33"] = c11
    values["C44"] = values["C55"] = values["C66"] = c44
    values["C12"] = values["C13"] = values["C23"] = c12
    values["Rho"] = rho
    return values


def write_h5(values, out_path, xmin, xmax, npts=2):
    """One group per property; each holds a constant (npts,npts,npts) field.

    npts=2 (the default) is enough since the field is uniform -- SEM3D's
    regular-grid reader (build_prop_files.F90::init_prop_file_field_Cstar)
    interpolates linearly, so any grid reproduces a constant exactly.
    """
    xmin = np.asarray(xmin, dtype=np.float64)
    xmax = np.asarray(xmax, dtype=np.float64)
    shape = (npts, npts, npts)  # numpy (Nz,Ny,Nx) <-> Fortran reader var(Nx,Ny,Nz)
    with h5py.File(out_path, "w") as f:
        for name in ELASTIC_NAMES:
            g = f.create_group(name)
            g.create_dataset("samples", data=np.full(shape, values[name]))
            g.attrs["xMinGlob"] = xmin
            g.attrs["xMaxGlob"] = xmax


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", help="output HDF5 path")
    p.add_argument("--vp", type=float, required=True)
    p.add_argument("--vs", type=float, required=True)
    p.add_argument("--rho", type=float, required=True)
    p.add_argument("--xmin", type=float, required=True)
    p.add_argument("--xmax", type=float, required=True)
    p.add_argument("--ymin", type=float, required=True)
    p.add_argument("--ymax", type=float, required=True)
    p.add_argument("--zmin", type=float, required=True)
    p.add_argument("--zmax", type=float, required=True)
    p.add_argument("--npts", type=int, default=2,
                    help="grid points per axis (default 2; field is constant)")
    args = p.parse_args()

    values = isotropic_cij(args.vp, args.vs, args.rho)
    write_h5(values, args.out,
              xmin=[args.xmin, args.ymin, args.zmin],
              xmax=[args.xmax, args.ymax, args.zmax],
              npts=args.npts)
    print("Wrote %s: C11=%.6g C44=%.6g C12=%.6g Rho=%.6g"
          % (args.out, values["C11"], values["C44"], values["C12"], values["Rho"]))


if __name__ == "__main__":
    main()
