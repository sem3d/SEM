#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_mpi.py - Compare SEM2D receiver traces between two runs (typically np=1 vs np=4)
and report the relative difference per receiver/column, with a PASS/FAIL verdict.

A correct MPI partition reproduces the serial result to round-off (~1e-12). Columns whose
peak amplitude is below --floor are treated as noise and skipped (their reldiff is
meaningless). Exit code 0 if every non-noise column is within --rtol, else 1.

Traces are read from the HDF5 station files (`capteurs.NNNN.h5`, one per rank -- each
receiver written by whichever rank owns it, so we merge across ranks by dataset name).
Falls back to the legacy ASCII `rec_*.vel` files if no h5 is present.

    python3 compare_mpi.py runA/traces runB/traces [--rtol 1e-9] [--floor 1e-11]
"""
import argparse
import glob
import os
import sys

import numpy as np


def _load_h5(d):
    """dict receiver_name -> (ntime, ncol) array from every capteurs.*.h5 in dir d.
    A receiver dataset is 2D (col0 = time); *_pos, Variables and Energy* are skipped."""
    import h5py
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "capteurs.*.h5"))):
        with h5py.File(f, "r") as h:
            for name, dset in h.items():
                if getattr(dset, "ndim", 0) != 2:
                    continue
                if name.endswith("_pos") or name.startswith("Energy") or name == "Variables":
                    continue
                out[name] = dset[()]
    return out


def load(d):
    """Receiver traces from h5 station files, or legacy rec_*.vel if no h5 present."""
    if glob.glob(os.path.join(d, "capteurs.*.h5")):
        return _load_h5(d)
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "rec_*.vel"))):
        if "prot" in f:
            continue
        out[os.path.basename(f)] = np.loadtxt(f)
    return out


def compare(dirA, dirB, rtol=1e-9, floor=1e-11):
    A, B = load(dirA), load(dirB)
    common = sorted(set(A) & set(B))
    if not common:
        print("  no common receiver traces found in the two dirs")
        return False, []
    rows, worst, ok = [], 0.0, True
    for name in common:
        a, b = A[name], B[name]
        if a.ndim == 1:
            a, b = a[:, None], b[:, None]
        n = min(len(a), len(b))
        for c in range(1, min(a.shape[1], b.shape[1])):  # col 0 = time
            peak = float(np.max(np.abs(a[:n, c])))
            d = float(np.max(np.abs(a[:n, c] - b[:n, c])))
            reld = d / peak if peak > 0 else 0.0
            noise = peak < floor
            if not noise:
                worst = max(worst, reld)
                if reld > rtol:
                    ok = False
            rows.append((name, c, peak, reld, noise))
    return ok, rows, worst


def main():
    p = argparse.ArgumentParser(description="Compare SEM2D traces of two runs (np=1 vs np=4).")
    p.add_argument("dirA")
    p.add_argument("dirB")
    p.add_argument("--rtol", type=float, default=1e-9, help="pass threshold on reldiff (default 1e-9)")
    p.add_argument("--floor", type=float, default=1e-11, help="peak below this = noise, skipped")
    p.add_argument("--label", default="", help="optional case label for the header")
    a = p.parse_args()
    res = compare(a.dirA, a.dirB, a.rtol, a.floor)
    if len(res) == 2:
        sys.exit(2)
    ok, rows, worst = res
    hdr = "compare {} vs {}".format(a.dirA, a.dirB)
    if a.label:
        hdr = "[{}] ".format(a.label) + hdr
    print(hdr)
    for name, c, peak, reld, noise in rows:
        tag = " (noise, skipped)" if noise else ""
        print("  {} col{}: reldiff={:.3e}  peak={:.3e}{}".format(name, c, reld, peak, tag))
    print("  worst non-noise reldiff = {:.3e}  -> {} (rtol={:g})".format(
        worst, "PASS" if ok else "FAIL", a.rtol))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
