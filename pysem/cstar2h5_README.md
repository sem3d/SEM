# cstar2h5

Convert a homogenisation **Cstar** binary (`homo/src/`, `icode = -82`, GLL grid)
into the **regular-grid anisotropic HDF5** material file read by SEM3D
(`MATDEF_HOOKE_ANISO` / `MATDEF_FLUID_ANISO`).

The Cstar tensor is stored on GLL nodes; SEM's regular-grid reader expects a
uniform grid, so the script performs the GLL → uniform Lagrange interpolation
internally (the same operation SEM does in
`SEM3D/SRC/build_prop_files.F90 :: init_prop_file_field_Cstar`).

## Usage

```bash
python cstar2h5.py <cstar_file> <out.h5> [--kind auto|elastic|acoustic] [--split]
python cstar2h5.py --selftest      # internal interpolation / mapping checks
```

* `--split` — write one file per property (`<out_base>_<name>.h5`, root
  `samples` dataset) instead of a single file with one group per property.

Requires `numpy` and `h5py`.

## Output layout

One group per property; each group holds dataset `samples` plus attributes
`xMinGlob`, `xMaxGlob`. `samples` has numpy shape `(Nz, Ny, Nx)` so the Fortran
HDF5 reader sees `var(Nx, Ny, Nz)` (C ↔ Fortran dimension reversal). Grid:
`Nx = nelx·ndeg + 1`, etc.; bounds `[0,0,0] … [nelx·xel, nely·yel, nelz·zel]`.

### Elastic (`ncomp = 22`, Nd = 6)
Groups: `C11,C22,C33,C44,C55,C66,C12,C13,C14,C15,C16,C23,C24,C25,C26,C34,C35,C36,C45,C46,C56,Rho`
— full stiffness tensor + density.

### Acoustic (`ncomp = 7`, Nd = 3) — Capdeville & Cancès (2015) density formulation
Groups: `iRho11,iRho22,iRho33,iRho12,iRho13,iRho23,iKappa`.
The `iRho..` components are the effective **inverse-density tensor** `ρ*⁻¹_ij`
(anisotropy lives in the density); `iKappa` is the scalar effective **inverse
bulk modulus** `1/κ*`. A `physical_meaning` attribute documents each group.

> These group names match `prop_field` `propName`s in
> `SEM3D/SRC/build_prop_files.F90` (`MATDEF_FLUID_ANISO` / `CSTAR_FLUID`),
> renamed from the former misleading `K11..K23`/`Rho`. SEM reads the fields by
> position, so the rename is purely a naming/contract change.
