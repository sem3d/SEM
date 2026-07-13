#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check.py WORKDIR  -- validate the 3D PML-heterogeneous-inherit feature. PASS iff:
  (1) the solver log shows PML materials inheriting the base FILE definition,
  (2) the material snapshot's Lambda is genuinely HETEROGENEOUS (spread across the model),
      i.e. the field was read per-GLL rather than a homogeneous placeholder,
  (3) the run completed without NaN.
Exit 0 = PASS. (The precise "frozen at the face" clamp is a secondary interpolation detail;
the output snapshot does not tag PML elements distinctly, so it is not isolated here.)
"""
import sys
import glob
import numpy as np
import h5py

d = sys.argv[1]
log = open(d + "/sem.log").read()

# (1) inheritance
n_inherit = log.count("inherits FILE definition from base")
# (3) completion / no NaN  (SEM3D end message: "fin du calcul sur processeurs")
done = ("fin du calcul" in log) or ("Execution completed" in log)
bad = ("NaN" in log) or ("Infinity" in log)

# (2) heterogeneous Lambda in the material snapshot
geo = sorted(glob.glob(d + "/res/geometry*.h5"))
spread = 0.0
if geo:
    with h5py.File(geo[0], "r") as f:
        L = f["Lamb"][()]
    spread = (L.max() - L.min()) / (abs(L.mean()) + 1e-30)

print("  PML materials inheriting FILE base : %d" % n_inherit)
print("  Lambda spread (max-min)/|mean|     : %.3f  (>~1 => heterogeneous field, not placeholder)" % spread)
print("  completed=%s  no-NaN=%s" % (done, not bad))
ok = (n_inherit > 0) and (spread > 0.5) and done and (not bad)
print("  ->", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
