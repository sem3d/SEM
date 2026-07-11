# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
check_nan.py - Scan a SEM (2D/3D) run for non-finite values (NaN/Inf).

Checks, in one pass over a results directory:
  * geometry<rg>.h5      -- static per-node fields (Lambda/Mu/Density/Mass/Jac/Cij/...)
  * Rsem<it>/sem_field.<rg>.h5 -- dynamic fields (veloc/displ/accel/pressure/...),
                           reporting the FIRST iteration where any field goes non-finite
  * rec_*.vel            -- receiver traces (text columns)

Exits non-zero if any NaN/Inf is found (handy in scripts / CI).

    Ex.1 : quick health check of a run's res/ directory
        python3 check_nan.py res/

    Ex.2 : also map the non-finite snapshot nodes to coordinates and say whether
           they fall inside the physical domain or in the PML (2D: xmin xmax zmin zmax)
        python3 check_nan.py res/ --locate --domain 0 500 0 300

    Ex.3 : only look at receiver traces
        python3 check_nan.py res/ --only traces

    Ex.4 : self-test (no data needed)
        python3 check_nan.py --selftest
"""
import argparse
import glob
import os
import sys

import numpy as np

osj = os.path.join


def _float_datasets(h5grp):
    """Yield (name, ndarray) for every floating-point dataset in an open h5 file/group."""
    import h5py as hf
    for name, item in h5grp.items():
        if isinstance(item, hf.Dataset) and item.dtype.kind == 'f':
            yield name, item[...]


def nonfinite_report(arr):
    """Return (n_bad, n_tot, absmax_finite) for an array."""
    finite = np.isfinite(arr)
    n_bad = int((~finite).sum())
    absmax = float(np.max(np.abs(arr[finite]))) if finite.any() else float('nan')
    return n_bad, int(arr.size), absmax


def check_geometry(wkd):
    """Scan geometry<rg>.h5 static fields. Returns list of (file, field, n_bad, absmax)."""
    import h5py as hf
    out = []
    for geo in sorted(glob.glob(osj(wkd, 'geometry*.h5'))):
        with hf.File(geo, 'r') as h5f:
            for name, arr in _float_datasets(h5f):
                n_bad, _, absmax = nonfinite_report(arr)
                if n_bad:
                    out.append((os.path.basename(geo), name, n_bad, absmax))
    return out


def _snap_iters(wkd):
    """Sorted list of snapshot iteration indices from Rsem<it>/ directories."""
    its = []
    for d in glob.glob(osj(wkd, 'Rsem*')):
        b = os.path.basename(d)
        if b.startswith('Rsem') and b[4:].isdigit():
            its.append(int(b[4:]))
    return sorted(its)


def check_snapshots(wkd):
    """Scan Rsem<it>/sem_field.<rg>.h5. Returns (first_bad_it, details, n_files) where
    details is a list of (it, rg, field, n_bad, absmax) for the FIRST bad iteration
    (or [] if clean); n_files is how many snapshot files were scanned."""
    import h5py as hf
    first_bad = None
    details = []
    n_files = 0
    for it in _snap_iters(wkd):
        bad_here = []
        for f in sorted(glob.glob(osj(wkd, 'Rsem{:04d}'.format(it), 'sem_field.*.h5'))):
            n_files += 1
            rg = int(os.path.basename(f).split('.')[1])
            with hf.File(f, 'r') as h5f:
                for name, arr in _float_datasets(h5f):
                    n_bad, _, absmax = nonfinite_report(arr)
                    if n_bad:
                        bad_here.append((it, rg, name, n_bad, absmax))
        if bad_here:
            first_bad = it
            details = bad_here
            break
    return first_bad, details, n_files


def check_traces(wkd):
    """Scan rec_*.vel text traces. Returns list of (file, col, n_bad, first_bad_row, absmax)."""
    out = []
    for f in sorted(glob.glob(osj(wkd, 'rec_*.vel'))) + sorted(glob.glob(osj(wkd, '**', 'rec_*.vel'), recursive=True)):
        if 'prot' in f:  # skip protection/checkpoint copies
            continue
        try:
            d = np.loadtxt(f)
        except Exception:
            continue
        if d.ndim == 1:
            d = d[:, None]
        for c in range(d.shape[1]):
            finite = np.isfinite(d[:, c])
            n_bad = int((~finite).sum())
            if n_bad:
                first_row = int(np.where(~finite)[0][0])
                absmax = float(np.max(np.abs(d[finite, c]))) if finite.any() else float('nan')
                out.append((os.path.relpath(f, wkd), c, n_bad, first_row, absmax))
    # de-dup (the two globs can overlap)
    seen, uniq = set(), []
    for r in out:
        if r[0] not in seen:
            seen.add(r[0]); uniq.append(r)
    return uniq


def locate_nan(wkd, domain=None):
    """For the first bad snapshot, map non-finite nodes to coordinates and classify
    inside/outside the physical domain. domain = (xmin,xmax,zmin,zmax) [+ ymin,ymax for 3D]."""
    import h5py as hf
    first_bad, details, _ = check_snapshots(wkd)
    if first_bad is None:
        print("  locate: no non-finite snapshot field found.")
        return
    # pick the rank/field of the first bad entry
    it, rg, field = details[0][0], details[0][1], details[0][2]
    geo = osj(wkd, 'geometry{:04d}.h5'.format(rg))
    snap = osj(wkd, 'Rsem{:04d}'.format(it), 'sem_field.{:04d}.h5'.format(rg))
    with hf.File(geo, 'r') as h5f:
        nodes = h5f['Nodes'][...]
    with hf.File(snap, 'r') as h5f:
        v = h5f[field][...]
    if nodes.shape[0] == 3 and nodes.shape[1] != 3:
        nodes = nodes.T
    if v.ndim == 2 and v.shape[0] != nodes.shape[0]:
        v = v.T
    bad = ~np.isfinite(v).reshape(v.shape[0], -1).all(axis=1)
    xz = nodes[:, :2]
    print("  locate: first bad snapshot it={} rg={} field='{}' -> {}/{} nodes non-finite"
          .format(it, rg, field, int(bad.sum()), len(bad)))
    if bad.any():
        b = xz[bad]
        print("    non-finite node bbox: x[{:.1f},{:.1f}]  z[{:.1f},{:.1f}]"
              .format(b[:, 0].min(), b[:, 0].max(), b[:, 1].min(), b[:, 1].max()))
        if domain is not None:
            xmin, xmax, zmin, zmax = domain[:4]
            inside = (b[:, 0] >= xmin) & (b[:, 0] <= xmax) & (b[:, 1] >= zmin) & (b[:, 1] <= zmax)
            n_in, n_out = int(inside.sum()), int((~inside).sum())
            print("    of these: {} inside physical domain, {} OUTSIDE (PML region)"
                  .format(n_in, n_out))
            if n_out and not n_in:
                print("    => non-finite values are PML-only (points to a PML setup/stability issue).")
            elif n_in:
                print("    => non-finite values reach the physical domain (not PML-confined).")


def run(wkd, only=None, locate=False, domain=None):
    only = only or ('geometry', 'snapshots', 'traces')
    any_bad = False

    if 'geometry' in only:
        g = check_geometry(wkd)
        if g:
            any_bad = True
            print("[geometry] NON-FINITE static fields:")
            for f, name, n_bad, absmax in g:
                print("  {}:/{}  {} non-finite  (finite |max|={:.3g})".format(f, name, n_bad, absmax))
        else:
            print("[geometry] clean (or no geometry*.h5 found)")

    if 'snapshots' in only:
        first_bad, details, n_files = check_snapshots(wkd)
        if first_bad is not None:
            any_bad = True
            print("[snapshots] FIRST non-finite at iteration Rsem{:04d}:".format(first_bad))
            for it, rg, name, n_bad, absmax in details:
                print("  rg={} /{}  {} non-finite  (finite |max|={:.3g})".format(rg, name, n_bad, absmax))
        elif n_files:
            print("[snapshots] clean ({} sem_field files scanned)".format(n_files))
        else:
            print("[snapshots] no Rsem*/sem_field*.h5 found")

    if 'traces' in only:
        t = check_traces(wkd)
        if t:
            any_bad = True
            print("[traces] NON-FINITE receiver columns:")
            for f, c, n_bad, first_row, absmax in t:
                print("  {}  col{}  {} non-finite (first at row {}, finite |max|={:.3g})"
                      .format(f, c, n_bad, first_row, absmax))
        else:
            print("[traces] clean (or no rec_*.vel found)")

    if locate:
        locate_nan(wkd, domain=domain)

    return any_bad


def _selftest():
    # nonfinite_report
    a = np.array([1.0, 2.0, np.nan, np.inf, -3.0])
    n_bad, n_tot, absmax = nonfinite_report(a)
    assert n_bad == 2 and n_tot == 5 and absmax == 3.0, (n_bad, n_tot, absmax)
    # clean array
    n_bad, _, absmax = nonfinite_report(np.arange(10.0))
    assert n_bad == 0 and absmax == 9.0
    print("check_nan selftest OK")


def main():
    p = argparse.ArgumentParser(description="Scan a SEM run for NaN/Inf.")
    p.add_argument('wkd', nargs='?', default='res', help="results directory (default: res)")
    p.add_argument('--only', choices=['geometry', 'snapshots', 'traces'], nargs='+',
                   help="restrict the scan to these categories")
    p.add_argument('--locate', action='store_true',
                   help="map non-finite snapshot nodes to coordinates")
    p.add_argument('--domain', type=float, nargs='+', metavar='B',
                   help="physical domain bounds for PML classification: 2D 'xmin xmax zmin zmax'")
    p.add_argument('--selftest', action='store_true', help="run internal checks and exit")
    a = p.parse_args()
    if a.selftest:
        _selftest(); return
    if not os.path.isdir(a.wkd):
        print("ERR: not a directory: {}".format(a.wkd)); sys.exit(2)
    any_bad = run(a.wkd, only=a.only, locate=a.locate, domain=a.domain)
    sys.exit(1 if any_bad else 0)


if __name__ == '__main__':
    main()
