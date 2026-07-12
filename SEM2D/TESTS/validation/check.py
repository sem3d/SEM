#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check.py CHECK ARGS...  -- per-item validators for the SEM2D validation suite.
Exit 0 = PASS, 1 = FAIL. Uses the material-property snapshot (geometry*.h5, item 8) and the
h5 station traces (capteurs.*.h5). Domain [0,500]x[0,300]; PML region is x<0 or x>500.
"""
import sys
import glob
import re
import numpy as np
import h5py

LX = 500.0


def _geo(d):
    """Return (coords Nx2, dict name->Nx array) from the first geometry*.h5 in dir d."""
    f = sorted(glob.glob(d + "/res/geometry*.h5"))[0]
    h = h5py.File(f, "r")
    nodes = h["Nodes"][()]
    xz = nodes[:, 0:2]
    fields = {k: h[k][()] for k in h if h[k].ndim == 1 and h[k].dtype.kind == "f"}
    return xz, fields


def _traces(d):
    out = {}
    for f in sorted(glob.glob(d + "/traces/capteurs.*.h5")):
        with h5py.File(f, "r") as h:
            for k, ds in h.items():
                if ds.ndim == 2 and not k.endswith("_pos") and k != "Variables" and not k.startswith("Energy"):
                    out[k] = ds[()]
    return out


def _totals_line(log):
    for ln in open(log):
        m = re.search(r"TOTAL\s+(\d+)\s+elem\s+(\d+)\s+gll\s+(\d+)\s+dof", ln)
        if m:
            return tuple(int(g) for g in m.groups())
    return None


# --- E: all-ranks domain totals identical np1 vs np4 ------------------------
def totals(d1, d4):
    t1, t4 = _totals_line(d1 + "/sem.log"), _totals_line(d4 + "/sem.log")
    print("  np1 TOTAL(elem,gll,dof)=%s ; np4=%s" % (t1, t4))
    ok = t1 is not None and t1 == t4
    return ok


# --- C: aniso->iso projection reaches the PML (Lambda/Mu non-zero & = interior) ----
def iso_proj(d):
    xz, F = _geo(d)
    lam, mu = F["Lambda"], F["Mu"]
    x = xz[:, 0]
    pml = (x < 0) | (x > LX)
    interior = ~pml
    li, mi = np.median(lam[interior]), np.median(mu[interior])
    lp, mp = lam[pml], mu[pml]
    print("  interior lambda=%.3e mu=%.3e ; PML lambda[min,max]=[%.3e,%.3e] mu=[%.3e,%.3e]"
          % (li, mi, lp.min(), lp.max(), mp.min(), mp.max()))
    # PML must carry the SAME iso moduli as the (homogeneous) interior -- not zero/placeholder
    ok = (lp.min() > 0.5 * li) and (abs(lp.max() - li) < 1e-3 * li) and \
         (mp.min() > 0.5 * mi) and (abs(mp.max() - mi) < 1e-3 * mi)
    return ok


# --- B: heterogeneous field inherited, frozen at the face in the PML axis ----
#     PML is isotropic -> the inherited field shows up in Lambda (iso projection), frozen at
#     the face along the PML axis (constant in the PML, = the value at the domain edge).
def frozen_face(d, lam_left, lam_right):
    xz, F = _geo(d)
    x, lam = xz[:, 0], F["Lambda"]
    left = lam[x < -1.0]          # strictly inside the left PML (x<0)
    right = lam[x > LX + 1.0]     # strictly inside the right PML (x>LX)
    inter = lam[(x > 1.0) & (x < LX - 1.0)]
    print("  expect frozen Lambda: left~%.3e right~%.3e" % (lam_left, lam_right))
    print("  left PML Lambda[min,max]=[%.3e,%.3e] ; right=[%.3e,%.3e] ; interior spread=%.3e..%.3e"
          % (left.min(), left.max(), right.min(), right.max(), inter.min(), inter.max()))
    tol = 0.02
    ok = (abs(left.mean() - lam_left) < tol * lam_left and left.ptp() < tol * lam_left and
          abs(right.mean() - lam_right) < tol * lam_right and right.ptp() < tol * lam_right and
          inter.ptp() > 0.1 * inter.mean())   # interior genuinely heterogeneous
    return ok


# --- D: iso-limit -- fluid-aniso(iso Cstar)+PML matches iso-fluid+PML at receivers ----
def iso_limit(dA, dB, rtol=1e-2):
    A, B = _traces(dA), _traces(dB)
    common = sorted(set(A) & set(B))
    if not common:
        print("  no common receivers"); return False
    worst = 0.0
    for k in common:
        a, b = A[k], B[k]
        n = min(len(a), len(b))
        for c in range(1, min(a.shape[1], b.shape[1])):
            pk = np.abs(a[:n, c]).max()
            if pk < 1e-12:
                continue
            worst = max(worst, np.abs(a[:n, c] - b[:n, c]).max() / pk)
    print("  worst receiver reldiff (aniso-iso-Cstar vs native iso) = %.3e" % worst)
    return worst < rtol


# --- D: rest state / stability -- no NaN, ran to completion --------------------
def rest_state(d):
    log = open(d + "/sem.log").read()
    bad = ("NaN" in log) or ("below lower" in log) or ("Infinity" in log)
    done = "Execution completed" in log
    print("  completed=%s  no-NaN=%s" % (done, not bad))
    return done and not bad


if __name__ == "__main__":
    name = sys.argv[1]
    fn = {"totals": totals, "iso_proj": iso_proj, "frozen_face": frozen_face,
          "iso_limit": iso_limit, "rest_state": rest_state}[name]
    args = []
    for a in sys.argv[2:]:
        try:
            args.append(float(a))
        except ValueError:
            args.append(a)
    ok = fn(*args)
    print("  ->", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
