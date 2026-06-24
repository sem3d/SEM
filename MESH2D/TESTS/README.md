Exercises every `mesher2D` input path: on-the-fly grid generation (with variations)
and the external-mesh readers (`.unv`, HDF5). Run them all after building SEM:

    ./run_all.sh [path/to/mesher2D]      # default ../../build/MESH2D/mesher2D

`run_all.sh` feeds each test's `mesh.input` to `mesher2D` on stdin, then reports the
`mesh4spec.NNNN.h5` files (and `material.input` for on-the-fly) produced.

`mesh.input` is the stdin menu answer: line 1 = NPROCS, line 2 = mesh source
(`1` on-the-fly, `3` .unv, `4` HDF5), and for `3`/`4` line 3 = the input filename.

--------------------------------------------------------------------------------
## Test matrix

| Dir              | Source       | NPROCS | What it checks |
|------------------|--------------|:------:|----------------|
| `onthefly_1mat`  | on the fly   | 1      | basic grid + `material.input` (1 solid), 10x6 |
| `onthefly_2layer`| on the fly   | 1      | 2 layers -> material 0 (solid, top) / 1 (fluid, bottom) |
| `onthefly_pml`   | on the fly   | 1      | PML rings (npml=2, W/E/bottom) + derived `P` materials |
| `onthefly_mpi`   | on the fly   | 4      | metis partition into 4 proc files, 20x12 |
| `unv_input`      | `.unv` (3)   | 1      | Ideas `.unv` reader (menu 3); mesh from `gen_unv_mesh.py` |
| `hdf5_input`     | HDF5 (4)     | 1      | HDF5 Quad reader (menu 4); mesh from `gen_hdf5_mesh.py` |

All meshes are the `[0,0.05] x [0,0.03]` box. Geometry is in `mat.dat`; materials
in `mater.in`. The two external-mesh tests carry a static `material.input` (the
mesher only generates `material.input` on the on-the-fly path).

--------------------------------------------------------------------------------
## Expected results

- **onthefly_1mat** — `Creating grid mesh 10 x 6`; `mesh4spec.0000.h5`;
  `material.input` = 1 material (`S ...`).
- **onthefly_2layer** — same grid; `material.input` = 2 materials (`S` then `F`);
  top 3 rows reference material 0, bottom 3 rows material 1.
- **onthefly_pml** — `Creating grid mesh 14 x 8 ... (npml=2)`; `material.input` has
  **6** materials: 1 base `S` + 5 derived `P` (left edge, right edge, bottom edge,
  bottom-left corner, bottom-right corner) and a PML block with 5 descriptor lines
  `F npow Apow Px Left Pz Down omegac kc`:
      left   : T T F F   (Px,Left ; no z-PML)
      right  : T F F F
      bottom : F F T T   (Pz,Down ; no x-PML)
      BL cnr : T T T T
      BR cnr : T F T T
- **onthefly_mpi** — `mesh4spec.0000.h5 .. mesh4spec.0003.h5` (4 files); the metis
  `OPT:` lines + a `Proc` field in `mesh4spec.h5`.
- **unv_input / hdf5_input** — `10 Quads`/`60 Quads` read, `mesh4spec.0000.h5`
  written. Parity check: both should yield the same partitioned mesh as
  `onthefly_1mat` (same 10x6 geometry, single material).

