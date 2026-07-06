#!/usr/bin/env python3
"""Isotropic-equivalent sanity check for the anisotropic-density fluid material.

Compares two SEM3D trace HDF5 files element by element:
  - one produced with the fluid-aniso material (deftype=Fluid_Aniso in
    material.spec, mater.in), and
  - one produced with the plain isotropic fluid material (no material.spec,
    mater_fluid_iso.in),
run on the *same* mesh / source / stations with the same Vp and Rho. Both use
the same DM_FLUID_CG domain (dom%aniso just flips the kernel internally).

Because the fluid-aniso constant case reduces to rho^{-1}_ij = (1/rho) delta_ij
and kappa = rho*Vp^2, the two runs must agree to machine precision. This script
walks every dataset common to both files and reports the maximum absolute and
RMS differences; it exits non-zero if any difference exceeds --tol.

Usage:
    python compare_aniso_vs_iso.py traces_aniso.h5 traces_iso.h5 [--tol 1e-9]

Typical workflow on the cluster:
    # 1) aniso run (mater.in + material.spec with deftype=Fluid_Aniso)
    cp mater.in material.input ; mpirun -n 4 sem3d.exe ; mv res res_aniso
    # 2) iso run (regular fluid; move material.spec out of the way first, or
    #    its deftype=Fluid_Aniso would still apply to material 0)
    mv material.spec material.spec.bak
    cp mater_fluid_iso.in material.input ; mpirun -n 4 sem3d.exe ; mv res res_iso
    mv material.spec.bak material.spec
    # 3) compare the trace files written under each res*/ directory
    python compare_aniso_vs_iso.py res_aniso/traces*.h5 res_iso/traces*.h5
"""
import sys
import argparse
import numpy as np
import h5py


def collect_datasets(h5obj, prefix=""):
    """Return {path: shape} for every dataset under an HDF5 group, recursively."""
    out = {}
    for key, item in h5obj.items():
        path = f"{prefix}/{key}"
        if isinstance(item, h5py.Group):
            out.update(collect_datasets(item, path))
        else:
            out[path] = item.shape
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("aniso_h5", help="trace HDF5 from the FluidAniso ('A') run")
    ap.add_argument("iso_h5", help="trace HDF5 from the regular Fluid ('F') run")
    ap.add_argument("--tol", type=float, default=1e-9,
                    help="max allowed |diff| / scale before failing (default 1e-9)")
    args = ap.parse_args()

    with h5py.File(args.aniso_h5, "r") as fa, h5py.File(args.iso_h5, "r") as fi:
        da = collect_datasets(fa)
        di = collect_datasets(fi)
        common = sorted(set(da) & set(di))
        only_a = sorted(set(da) - set(di))
        only_i = sorted(set(di) - set(da))

        if only_a:
            print(f"WARNING: {len(only_a)} dataset(s) only in aniso file (e.g. {only_a[:3]})")
        if only_i:
            print(f"WARNING: {len(only_i)} dataset(s) only in iso file (e.g. {only_i[:3]})")
        if not common:
            print("ERROR: no common datasets to compare.")
            return 2

        worst = 0.0
        worst_name = None
        n_fail = 0
        print(f"{'dataset':<40} {'max|diff|':>12} {'rms|diff|':>12} {'rel':>10}")
        for name in common:
            if da[name] != di[name]:
                print(f"{name:<40} SHAPE MISMATCH {da[name]} vs {di[name]}")
                n_fail += 1
                continue
            a = np.asarray(fa[name][()], dtype=np.float64)
            b = np.asarray(fi[name][()], dtype=np.float64)
            if a.size == 0:
                continue
            diff = np.abs(a - b)
            scale = max(np.max(np.abs(a)), np.max(np.abs(b)), 1e-300)
            maxd = float(np.max(diff))
            rmsd = float(np.sqrt(np.mean(diff**2)))
            rel = maxd / scale
            if rel > worst:
                worst, worst_name = rel, name
            if rel > args.tol:
                n_fail += 1
            print(f"{name:<40} {maxd:12.4e} {rmsd:12.4e} {rel:10.2e}")

        print("-" * 78)
        print(f"worst relative difference: {worst:.3e}  ({worst_name})")
        if n_fail == 0 and worst <= args.tol:
            print(f"PASS: anisotropic-density domain matches regular fluid (tol={args.tol:g}).")
            return 0
        print(f"FAIL: {n_fail} dataset(s) exceed tol={args.tol:g}.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
