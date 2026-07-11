# SEM2D MPI comparison suite (np=1 vs np=4)

Validates that an MPI-partitioned run reproduces the serial run. For every
**physics × mesh × PML** case it builds the mesh for NPROCS=1 and NPROCS=4, runs the
mesher then `sem2d` for each, and compares the receiver traces. A correct MPI partition
matches serial to round-off (~1e-12); anything larger is a real MPI defect.

## Run
```sh
./run_all.sh [BUILD_DIR] [PHYS_FILTER] [MESH_FILTER]
#   BUILD_DIR   default ../../../build
#   PHYS_FILTER e.g. solid, fluid, sf_   (substring; default all)
#   MESH_FILTER e.g. quad8               (substring; default all)
# env: RTOL (pass threshold, default 1e-6), PYTHON
```
Everything is generated under `_work/` (wiped each run). Needs `mpirun`, `h5py`, numpy.

## Files
- `mkcase.py`     — generate one self-contained case (mesh + materials + input.spec + capteurs).
- `compare_mpi.py`— compare two trace dirs; per-receiver reldiff + PASS/FAIL (skips noise columns).
- `run_all.sh`    — driver over the whole matrix; prints a verdict table.

## Matrix
Physics: `solid`, `solid_aniso`, `fluid`, `fluid_aniso`, `sf_iso` (solid+fluid coupling),
`sf_aniso` (solid + anisotropic-density fluid coupling).
Mesh: `onthefly` (mat.dat), `quad4` (external HDF5), `quad8` (external HDF5, 2nd order).
PML: `0` / `1`. Domain [0,500]×[0,300], base grid 10×6, sources/receivers at element
midpoints (an on-boundary source/receiver is inconsistent across ranks — a harness
pitfall, not a solver bug).

## Status (2026-07-11, local Debug build, after the PML wall-flag fix)
`worst non-noise reldiff`, serial vs 4-rank. **quad8 == quad4 == onthefly in every case**
(the mesh order never changes the verdict → the 2nd-order/MPI partition is correct).

| physics       | PML=0        | PML=1        |
|---------------|--------------|--------------|
| solid         | ✅ ~3e-14    | ✅ ~8e-15    |
| solid_aniso   | ✅ ~3e-14    | ✅ ~1e-14    |
| fluid         | ✅ ~7e-14    | ✅ ~2e-14    |
| fluid_aniso   | ✅ ~1e-14    | ✅ ~1e-14    |
| sf_iso        | ✅ ~2e-14    | ✅ ~2e-13    |
| sf_aniso      | ❌ ~2e-2     | ❌ ~2e-2     |

### Reading the table
- **Everything except `sf_aniso` is bit-exact**, PML included.
- The original PML+MPI residual (solid/fluid ~1.6e-3, sf_iso+PML ~9e-2) was caused by
  wrong PML/Reflex flags on partition-wall faces (`Near_Element(1) = -1` locally made
  them look like domain boundaries) which cascaded to their vertices — fixed in
  `PML_def.F90` by cross-rank flag handshakes; see
  `log/plans/2026-07-11_pml-vertex-mpi-exchange.md` for the full diagnosis.
- **`sf_aniso` fails even without PML** (~2e-2): the split-potential solid↔aniso-fluid
  coupling is not MPI-consistent across ranks — the known residual from
  `2026-07-09_sf-coupling-2d-solution-plan.md` (its interface DOFs aren't in the sWall comm
  tables; needs an interface-preserving partition or cross-rank SF exchange). This is the
  only remaining ❌ and has its own plan.
