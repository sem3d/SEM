#!/usr/bin/env python3
"""Generate test HDF5 material file for fluid anisotropic domain (constant Kij).

Values match TEST_0009: Vp=1500 m/s, Rho=1000 kg/m3
  K11 = K22 = K33 = rho * Vp^2 = 2.25e9 Pa
  K12 = K13 = K23 = 0 (isotropic)
"""
import numpy as np
import h5py

OUTFILE = "mat_fluid_aniso.h5"
XMIN = np.array([0.0, 0.0, 0.0])
XMAX = np.array([500.0, 500.0, 500.0])
NPTS = 2  # 2 points per axis: minimum for trilinear interpolation

K_diag = 1000.0 * 1500.0**2   # 2.25e9 Pa
K_off  = 0.0
RHO    = 1000.0                # kg/m3

components = {
    "K11": K_diag,
    "K22": K_diag,
    "K33": K_diag,
    "K12": K_off,
    "K13": K_off,
    "K23": K_off,
    "Rho": RHO,
}

shape = (NPTS, NPTS, NPTS)

with h5py.File(OUTFILE, "w") as f:
    for name, value in components.items():
        grp = f.create_group(name)
        grp.attrs["xMinGlob"] = XMIN
        grp.attrs["xMaxGlob"] = XMAX
        data = np.full(shape, value, dtype=np.float64)
        grp.create_dataset("samples", data=data)

print(f"Written: {OUTFILE}")
