#!/usr/bin/env python3
"""Generate the HDF5 material file for the anisotropic-DENSITY acoustic domain.

Density formulation (Capdeville & Cance 2015): the anisotropy is carried by the
density, not by the bulk modulus. The seven datasets are the inverse-density
tensor rho^{-1}_ij and the inverse bulk modulus 1/kappa. For I/O compatibility
the legacy HDF5 group keys K11..K23 / Rho are kept, but they now hold:

  K11 = rho^{-1}_11   K22 = rho^{-1}_22   K33 = rho^{-1}_33
  K12 = rho^{-1}_12   K13 = rho^{-1}_13   K23 = rho^{-1}_23
  Rho = 1/kappa       (inverse bulk modulus, NOT a density)

We reproduce the same orthorhombic wave speeds as the legacy stiffness test
(Vx=1500, Vy=1200, Vz=900 m/s) but via density anisotropy, by choosing a
uniform bulk modulus kappa and setting rho^{-1}_ii = V_i^2 / kappa, so that
V_i = sqrt(kappa * rho^{-1}_ii). With kappa = 2.25e9 Pa this gives an
anisotropic density (rho_11=1000, rho_22=1562.5, rho_33=2777.8 kg/m^3).
"""
import numpy as np
import h5py

OUTFILE = "mat_fluid_aniso.h5"
XMIN = np.array([0.0, 0.0, 0.0])
XMAX = np.array([500.0, 500.0, 500.0])
NPTS = 2  # 2 points per axis: minimum for trilinear interpolation

KAPPA = 2.25e9                 # uniform scalar bulk modulus [Pa]
VX, VY, VZ = 1500.0, 1200.0, 900.0

inv_kappa = 1.0 / KAPPA
components = {
    "iRho11": VX**2 / KAPPA,   # rho^{-1}_11  -> rho_11 = 1000.0   kg/m^3
    "iRho22": VY**2 / KAPPA,   # rho^{-1}_22  -> rho_22 = 1562.5   kg/m^3
    "iRho33": VZ**2 / KAPPA,   # rho^{-1}_33  -> rho_33 = 2777.78  kg/m^3
    "iRho12": 0.0,
    "iRho13": 0.0,
    "iRho23": 0.0,
    "iKappa": inv_kappa,       # 1/kappa  (solver inverts: lambda = 1/iKappa)
}

shape = (NPTS, NPTS, NPTS)

with h5py.File(OUTFILE, "w") as f:
    for name, value in components.items():
        grp = f.create_group(name)
        grp.attrs["xMinGlob"] = XMIN
        grp.attrs["xMaxGlob"] = XMAX
        data = np.full(shape, value, dtype=np.float64)
        grp.create_dataset("samples", data=data)

print(f"Written: {OUTFILE}   (kappa = {KAPPA:.3e} Pa, density formulation)")
for name, value in components.items():
    print(f"  {name} = {value:.6e}")
print("Wave speeds: Vx=%.0f Vy=%.0f Vz=%.0f m/s  (V_i = sqrt(kappa * rho^-1_ii))"
      % (np.sqrt(KAPPA*components['iRho11']), np.sqrt(KAPPA*components['iRho22']),
         np.sqrt(KAPPA*components['iRho33'])))
