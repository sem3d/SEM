# SEM3D PML heterogeneous-inherit test (pml_hetero)

Validates `2026-07-10_pml-heterogeneous-inherit`: a PML whose base material is a heterogeneous
FILE field (`Kappa_Mu_Rho`, `spacedef=file`) must **inherit** that FILE definition and read the
same per-GLL field, instead of falling back to a homogeneous placeholder.

## Run
```sh
./run_all.sh [BUILD_DIR]      # default ../../../build ; needs mpirun, h5py, numpy
```
Reuses the `NON-REGR/TEST_0006` cube geometry (26 `copy=0; solidpml` PML materials) but swaps
the random field for a CONTROLLED gradient (`gen_field.py`: Kappa linear in x, Mu/Rho const),
runs a short sim with the material snapshot on, and checks with `check.py`. Output under
`_work/` (wiped each run).

## What it asserts (PASS)
1. The solver log shows the PML materials inheriting the base FILE definition
   (`"PML material N inherits FILE definition from base 0 deftype 3"`, x26).
2. The material snapshot's `Lambda` is genuinely HETEROGENEOUS across the model
   (spread (max-min)/|mean| ≈ 2.6) — i.e. the field is read per-GLL, not a placeholder.
   (Mu/Rho are constant by construction, a sanity anchor.)
3. The run completes without NaN.

## Status (2026-07-12, local build)
PASS — 26 PML materials inherit the FILE base; Lambda spread 2.574; completed.

## Note
The precise "frozen at the face" clamp (property constant along the PML depth axis) is a
secondary interpolation detail and is NOT isolated here: the geometry snapshot does not tag PML
elements distinctly from the base solid, so a clean per-axis frozen check would need extra
instrumentation. The **core** feature (PML inherits the heterogeneous FILE field, not a
placeholder) is what this test locks in.
