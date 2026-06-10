#!/usr/bin/env python3
"""Generate test HDF5 material file for fluid anisotropic domain.

Anisotropic: different wave speeds along each axis.
  Vx = 1500 m/s  ->  K11 = rho * Vx^2 = 2.250e9 Pa
  Vy = 1200 m/s  ->  K22 = rho * Vy^2 = 1.440e9 Pa
  Vz =  900 m/s  ->  K33 = rho * Vz^2 = 8.100e8 Pa
  K12 = K13 = K23 = 0  (orthorhombic, no axis coupling)
  Rho = 1000 kg/m3
"""
import numpy as np
import h5py

OUTFILE = "mat_fluid_aniso.h5"
XMIN = np.array([0.0, 0.0, 0.0])
XMAX = np.array([500.0, 500.0, 500.0])
NPTS = 2  # 2 points per axis: minimum for trilinear interpolation

RHO = 1000.0
VX, VY, VZ = 1500.0, 1200.0, 900.0

components = {
    "K11": RHO * VX**2,   # 2.250e9 Pa
    "K22": RHO * VY**2,   # 1.440e9 Pa
    "K33": RHO * VZ**2,   # 8.100e8 Pa
    "K12": 0.0,
    "K13": 0.0,
    "K23": 0.0,
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
for name, value in components.items():
    print(f"  {name} = {value:.4e}")
