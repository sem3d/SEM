#!/usr/bin/env python3
"""Generate 2D Cstar.h5 test files for the SEM2D anisotropic-from-file path
(Phase 2). Produces the exact layout the reader build_prop_files_2d.F90 expects:

  - one HDF5 GROUP per property; each group has
      * dataset "samples"  (2D uniform grid, dims NNx x NNz)
      * attributes "xMinGlob","xMaxGlob"  (length 2: [xmin,zmin],[xmax,zmax])
  - elastic (Hooke_Aniso) groups : C11,C12,C13,C22,C23,C33,Rho   (Voigt 2D 1=xx,2=zz,3=xz)
  - fluid   (Fluid_Aniso)  groups : iRho11,iRho12,iRho22,iKappa   (inverse-density tensor + 1/kappa)

All test fields are CONSTANT per property (uniform medium) so the HDF5 storage
order (C vs Fortran) is irrelevant -- this keeps the tests robust.

Cases:
  elastic_iso  : Cij built from an isotropic medium (lambda,mu,rho) -> must match a
                 plain isotropic material.input run (Block C isotropic-limit gate).
  fluid_iso    : iRho = (1/rho) I, iKappa = 1/kappa -> must match isotropic acoustic
                 (Block D isotropic-limit gate).
  fluid_aniso  : different iRho11 != iRho22 -> directional acoustic speeds.

Usage:  python3 gen_test_cstar2d.py            # writes all three .h5 in ./
        python3 gen_test_cstar2d.py 0 0.05 0 0.03   # custom domain [xmin xmax zmin zmax]
"""
import sys
import numpy as np
import h5py

# ---- domain bounding box (must cover the SEM2D mesh) -----------------------
if len(sys.argv) == 5:
    XMIN, XMAX, ZMIN, ZMAX = (float(a) for a in sys.argv[1:5])
else:
    XMIN, XMAX, ZMIN, ZMAX = 0.0, 0.05, 0.0, 0.03   # default test domain
NN = 2  # 2x2 uniform grid is enough for a constant field

xmin = np.array([XMIN, ZMIN], dtype=np.float64)
xmax = np.array([XMAX, ZMAX], dtype=np.float64)


def write_cstar(fname, components):
    """components: dict name->scalar (constant field)."""
    shape = (NN, NN)
    with h5py.File(fname, "w") as f:
        for name, value in components.items():
            g = f.create_group(name)
            g.attrs["xMinGlob"] = xmin          # length-2 (matches reader)
            g.attrs["xMaxGlob"] = xmax
            g.create_dataset("samples", data=np.full(shape, value, dtype=np.float64))
    print(f"written {fname}:")
    for k, v in components.items():
        print(f"   {k:8s} = {v:.6e}")


# ---------------------------------------------------------------------------
# Case 1: ELASTIC isotropic limit (Hooke_Aniso). Cij from (lambda,mu).
#   Voigt 2D plane strain: C11=C22=lambda+2mu, C12=lambda, C33=mu, C13=C23=0.
def case_elastic_iso():
    rho = 2000.0
    vp, vs = 3000.0, 1700.0
    mu = rho * vs**2
    lam = rho * vp**2 - 2.0 * mu
    write_cstar("Cstar_elastic_iso.h5", {
        "C11": lam + 2 * mu, "C12": lam, "C13": 0.0,
        "C22": lam + 2 * mu, "C23": 0.0,
        "C33": mu,
        "Rho": rho,
    })


# Case 2: FLUID isotropic limit (Fluid_Aniso). iRho = (1/rho) I, iKappa = 1/kappa.
def case_fluid_iso():
    rho = 1000.0
    vp = 1500.0
    kappa = rho * vp**2          # 2.25e9
    write_cstar("Cstar_fluid_iso.h5", {
        "iRho11": 1.0 / rho, "iRho12": 0.0, "iRho22": 1.0 / rho,
        "iKappa": 1.0 / kappa,
    })


# Case 3: FLUID anisotropic density. Different inverse density along x and z.
def case_fluid_aniso():
    rho_x, rho_z = 1000.0, 1600.0   # heavier in z -> slower vertical speed
    vp = 1500.0
    kappa = 1000.0 * vp**2
    write_cstar("Cstar_fluid_aniso.h5", {
        "iRho11": 1.0 / rho_x, "iRho12": 0.0, "iRho22": 1.0 / rho_z,
        "iKappa": 1.0 / kappa,
    })


if __name__ == "__main__":
    print(f"domain [x] {XMIN}..{XMAX}  [z] {ZMIN}..{ZMAX}")
    case_elastic_iso()
    case_fluid_iso()
    case_fluid_aniso()
    print("\nRename the chosen file to the name in material.spec (e.g. Cstar.h5).")
