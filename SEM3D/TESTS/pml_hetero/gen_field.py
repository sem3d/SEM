#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_field.py DIR  -- write a CONTROLLED heterogeneous Kappa/Mu/Density field (h5) for the 3D
PML-heterogeneous-inherit test. Kappa varies LINEARLY in x across the physical box; Mu and
Density constant. Each file has a group <PropName> with dataset 'samples' (3D grid) and attrs
xMinGlob/xMaxGlob. The field box = the PHYSICAL domain, so PML GLL points (outside it) sample
the field CLAMPED at the face -> the material snapshot must show the PML "frozen at the face".

Physical box below is the non-PML interior of the TEST_0006 geometry ([-100,500] mesh, 100 of
PML on each side -> physical [0,400]). Prints expected frozen Kappa at x=XMIN / x=XMAX.
"""
import sys
import numpy as np
import h5py

XMIN, XMAX = 0.0, 400.0
YMIN, YMAX = 0.0, 400.0
ZMIN, ZMAX = 0.0, 400.0
NN = 9
K0, K1 = 1.0e10, 4.0e10             # Kappa(x=XMIN), Kappa(x=XMAX)
MU  = 1.75e10
RHO = 2800.0


def write_prop(path, name, grid):
    with h5py.File(path, "w") as f:
        g = f.create_group(name)
        g.create_dataset("samples", data=grid.astype(np.float64))
        g.attrs["xMinGlob"] = np.array([XMIN, YMIN, ZMIN], dtype=np.float64)
        g.attrs["xMaxGlob"] = np.array([XMAX, YMAX, ZMAX], dtype=np.float64)


def main(d):
    x = np.linspace(XMIN, XMAX, NN)
    kap = K0 + (K1 - K0) * (x - XMIN) / (XMAX - XMIN)
    kgrid = np.repeat(kap[:, None, None], NN, axis=1).repeat(NN, axis=2)
    # group name = the deftype property name (Kappa_Mu_Rho): "Kappa","Mu","Rho" (NOT "Density")
    write_prop(d + "/Mat_0_Kappa.h5",   "Kappa", kgrid)
    write_prop(d + "/Mat_0_Mu.h5",      "Mu",    np.full((NN, NN, NN), MU))
    write_prop(d + "/Mat_0_Density.h5", "Rho",   np.full((NN, NN, NN), RHO))
    print("%.6e %.6e" % (K0, K1))


if __name__ == "__main__":
    main(sys.argv[1])
