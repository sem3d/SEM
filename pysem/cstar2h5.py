# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
cstar2h5.py - Convert a homo (homofft) Cstar binary into the regular-grid
anisotropic HDF5 material file consumed by SEM3D.

Pipeline
--------
The homogenisation code (`homo/src/`) writes the effective Cstar tensor on a
*GLL* element grid (icode = -82). SEM3D, when a material uses
MATDEF_HOOKE_ANISO (elastic) or MATDEF_FLUID_ANISO (acoustic), reads a
*regular* grid HDF5 file (one group per property, each with `samples` +
`xMinGlob`/`xMaxGlob`). This script does the GLL -> regular interpolation
(exactly as SEM3D/SRC/build_prop_files.F90 :: init_prop_file_field_Cstar does
internally) and writes that HDF5.

Cases (driven by the component count ncomp = Nd*(Nd+1)/2 + 1):
  * ncomp = 22  -> ELASTIC  (Nd=6): full stiffness tensor + density.
                   Groups: C11,C22,C33,C44,C55,C66,C12,...,C56,Rho
  * ncomp = 7   -> ACOUSTIC (Nd=3): Capdeville & Cances (2015) density
                   formulation. The 6-component tensor is the effective
                   inverse-density rho*^{-1}_ij and the 7th value is the
                   effective inverse bulk modulus 1/kappa*.
                   Groups: iRho11,iRho22,iRho33,iRho12,iRho13,iRho23,iKappa

Usage
-----
    python cstar2h5.py <cstar_file> <out.h5> [--kind auto|elastic|acoustic]
                       [--split] [--selftest]

  --split     write one HDF5 file per property (<out_base>_<name>.h5 with a
              root `samples` dataset), matching generate_h5_materials.py.
              Default is a single file with one group per property.
  --selftest  run internal interpolation/round-trip checks and exit.
"""
import argparse
import os
import sys
import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None

__author__ = "SEMN tooling"
__license__ = "GPL"

# ---------------------------------------------------------------------------
# Property naming and raw -> named component mapping.
#
# The Cstar byte layout is upper-triangle row-major. SEM stores the mapping
# from that raw layout to its internal prop_field order (build_prop_files.F90).
# mapping[j] (1-indexed) gives the prop_field index for raw component j, i.e.
# prop_field(mapping[j]) <- raw[j]. We reuse it to assign the group name.
# ---------------------------------------------------------------------------
ELASTIC_NAMES = [
    "C11", "C22", "C33", "C44", "C55", "C66",
    "C12", "C13", "C14", "C15", "C16",
    "C23", "C24", "C25", "C26",
    "C34", "C35", "C36",
    "C45", "C46",
    "C56",
    "Rho",
]
ELASTIC_MAPPING = [1, 7, 8, 9, 10, 11, 2, 12, 13, 14, 15, 3, 16, 17, 18, 4, 19, 20, 5, 21, 6, 22]

ACOUSTIC_NAMES = [
    "iRho11", "iRho22", "iRho33", "iRho12", "iRho13", "iRho23", "iKappa",
]
ACOUSTIC_MAPPING = [1, 4, 5, 2, 6, 3, 7]

ACOUSTIC_MEANING = {
    "iRho11": "effective inverse-density tensor component rho*^{-1}_11",
    "iRho22": "effective inverse-density tensor component rho*^{-1}_22",
    "iRho33": "effective inverse-density tensor component rho*^{-1}_33",
    "iRho12": "effective inverse-density tensor component rho*^{-1}_12",
    "iRho13": "effective inverse-density tensor component rho*^{-1}_13",
    "iRho23": "effective inverse-density tensor component rho*^{-1}_23",
    "iKappa": "effective inverse bulk modulus 1/kappa*",
}


def ordered_names(ncomp):
    """Return the group name for each raw component (index 0..ncomp-1)."""
    if ncomp == 22:
        names, mapping = ELASTIC_NAMES, ELASTIC_MAPPING
    elif ncomp == 7:
        names, mapping = ACOUSTIC_NAMES, ACOUSTIC_MAPPING
    else:
        raise ValueError("Unsupported ncomp=%d (expected 7 or 22)" % ncomp)
    # raw component j -> prop_field(mapping[j]) -> name = names[mapping[j]-1]
    return [names[m - 1] for m in mapping]


# ---------------------------------------------------------------------------
# GLL nodes and Lagrange interpolation (GLL -> uniform).
# ---------------------------------------------------------------------------
def gll_nodes(ndeg):
    """Gauss-Lobatto-Legendre nodes on [-1, 1], ascending (matches ZELEGL)."""
    if ndeg < 1:
        raise ValueError("ndeg must be >= 1")
    if ndeg == 1:
        return np.array([-1.0, 1.0])
    coef = np.zeros(ndeg + 1)
    coef[ndeg] = 1.0
    dcoef = np.polynomial.legendre.legder(coef)
    interior = np.sort(np.polynomial.legendre.legroots(dcoef))
    return np.concatenate(([-1.0], interior, [1.0]))


def lagrange_matrix(nodes, xeval):
    """H[ig, ie] = L_ig(xeval[ie]) for Lagrange basis on `nodes`."""
    n = len(nodes)
    H = np.ones((n, len(xeval)))
    for ig in range(n):
        for m in range(n):
            if m == ig:
                continue
            H[ig, :] *= (xeval - nodes[m]) / (nodes[ig] - nodes[m])
    return H


def gll_to_uniform_matrix(ndeg):
    """1D operator mapping GLL nodal values to uniform sub-grid values."""
    nodes = gll_nodes(ndeg)
    xu = -1.0 + 2.0 * np.arange(ndeg + 1) / ndeg
    return lagrange_matrix(nodes, xu)  # shape (ndeg+1, ndeg+1): H[ig, iu]


def interp_element(buf, H):
    """Tensor-product GLL -> uniform interpolation of one element.

    buf : (ncomp, ndeg+1, ndeg+1, ndeg+1)  GLL nodal values (x,y,z)
    H   : (ndeg+1, ndeg+1)                  H[ig, iu]
    returns (ncomp, ndeg+1, ndeg+1, ndeg+1) on the uniform sub-grid.
    """
    out = np.einsum("cijk,iu->cujk", buf, H, optimize=True)
    out = np.einsum("cujk,jv->cuvk", out, H, optimize=True)
    out = np.einsum("cuvk,kw->cuvw", out, H, optimize=True)
    return out


# ---------------------------------------------------------------------------
# Cstar binary reader (icode = -82, direct-access unformatted).
# ---------------------------------------------------------------------------
class CstarHeader(object):
    __slots__ = ("icode", "iheader", "stored_len", "Nd", "ndeg",
                 "nelx", "nely", "nelz", "xel", "yel", "zel",
                 "xs_whole", "ys_whole", "zs_whole", "endian", "ncomp",
                 "reclen_bytes")


HEADER_BYTES = 8 * 4 + 6 * 4  # 8 int32 + 6 float32


def _parse_header(raw, endian):
    ints = np.frombuffer(raw[:32], dtype=endian + "i4")
    reals = np.frombuffer(raw[32:56], dtype=endian + "f4")
    h = CstarHeader()
    (h.icode, h.iheader, h.stored_len, h.Nd, h.ndeg,
     h.nelx, h.nely, h.nelz) = (int(v) for v in ints)
    (h.xel, h.yel, h.zel, h.xs_whole, h.ys_whole, h.zs_whole) = \
        (float(v) for v in reals)
    h.endian = endian
    return h


def read_cstar(path):
    """Read a homofft Cstar (-82) file into (header, fields).

    fields : float64 array (ncomp, NNx, NNy, NNz) on the global *uniform* grid,
             where NN* = nel* * ndeg + 1.
    """
    filesize = os.path.getsize(path)
    with open(path, "rb") as fh:
        raw = fh.read(HEADER_BYTES)
        if len(raw) < HEADER_BYTES:
            raise ValueError("File too small to contain a Cstar header")

        header = None
        for endian in ("<", ">"):
            cand = _parse_header(raw, endian)
            if cand.icode != -82:
                continue
            ncomp = cand.Nd * (cand.Nd + 1) // 2 + 1
            if ncomp not in (7, 22) or cand.ndeg < 1:
                continue
            reclen = ncomp * (cand.ndeg + 1) ** 3 * 4
            nrec = cand.iheader + cand.nelx * cand.nely * cand.nelz
            if reclen * nrec == filesize:
                cand.ncomp = ncomp
                cand.reclen_bytes = reclen
                header = cand
                break
        if header is None:
            raise ValueError(
                "Could not parse %s as a Cstar (-82) file: header/size "
                "inconsistent for both endiannesses." % path)

        h = header
        npte = h.ndeg + 1
        comp_block = h.ncomp * npte ** 3
        NNx = h.nelx * h.ndeg + 1
        NNy = h.nely * h.ndeg + 1
        NNz = h.nelz * h.ndeg + 1

        H = gll_to_uniform_matrix(h.ndeg).astype(np.float64)
        fields = np.empty((h.ncomp, NNx, NNy, NNz), dtype=np.float64)
        dtype = np.dtype(h.endian + "f4")

        for iz in range(h.nelz):
            for iy in range(h.nely):
                for ix in range(h.nelx):
                    irec = h.iheader + (ix + 1) + iy * h.nelx + iz * h.nelx * h.nely
                    fh.seek((irec - 1) * h.reclen_bytes)
                    flat = np.frombuffer(fh.read(comp_block * 4), dtype=dtype)
                    buf = flat.reshape((h.ncomp, npte, npte, npte), order="F").astype(np.float64)
                    sub = interp_element(buf, H)
                    gx, gy, gz = ix * h.ndeg, iy * h.ndeg, iz * h.ndeg
                    fields[:, gx:gx + npte, gy:gy + npte, gz:gz + npte] = sub

    return header, fields


# ---------------------------------------------------------------------------
# HDF5 writer.
# ---------------------------------------------------------------------------
def write_h5(header, fields, out_path, split=False):
    if h5py is None:
        raise RuntimeError("h5py is required to write HDF5 files")

    h = header
    names = ordered_names(h.ncomp)
    xmin = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    xmax = np.array([h.nelx * h.xel, h.nely * h.yel, h.nelz * h.zel],
                    dtype=np.float64)

    def samples_view(arr):
        # SEM3D (Fortran HDF5) reads var(Nx,Ny,Nz); h5py/C stores reversed,
        # so the dataset must have numpy shape (Nz,Ny,Nx).
        return np.ascontiguousarray(np.transpose(arr, (2, 1, 0)))

    if split:
        base = out_path[:-3] if out_path.endswith(".h5") else out_path
        for j, name in enumerate(names):
            with h5py.File("%s_%s.h5" % (base, name), "w") as f:
                f.create_dataset("samples", data=samples_view(fields[j]))
                f.attrs["xMinGlob"] = xmin
                f.attrs["xMaxGlob"] = xmax
                if name in ACOUSTIC_MEANING:
                    f.attrs["physical_meaning"] = ACOUSTIC_MEANING[name]
        return

    with h5py.File(out_path, "w") as f:
        for j, name in enumerate(names):
            g = f.create_group(name)
            g.create_dataset("samples", data=samples_view(fields[j]))
            g.attrs["xMinGlob"] = xmin
            g.attrs["xMaxGlob"] = xmax
            if name in ACOUSTIC_MEANING:
                g.attrs["physical_meaning"] = ACOUSTIC_MEANING[name]


# ---------------------------------------------------------------------------
# Self-test.
# ---------------------------------------------------------------------------
def _selftest():
    # 1) ndeg = 1: GLL == uniform -> H is identity.
    H1 = gll_to_uniform_matrix(1)
    assert np.allclose(H1, np.eye(2)), "ndeg=1 operator must be identity"

    # 2) A polynomial of degree <= ndeg is reproduced exactly by GLL interp.
    for ndeg in (2, 4, 6):
        nodes = gll_nodes(ndeg)
        xu = -1.0 + 2.0 * np.arange(ndeg + 1) / ndeg
        H = gll_to_uniform_matrix(ndeg)
        coef = np.random.RandomState(0).randn(ndeg + 1)
        poly = np.polynomial.polynomial.Polynomial(coef)
        f_gll = poly(nodes)
        f_uni = H.T @ f_gll  # values at uniform points
        assert np.allclose(f_uni, poly(xu)), "degree-%d interp failed" % ndeg

    # 3) Constant 3D field round-trips through the element interpolation.
    ndeg = 3
    H = gll_to_uniform_matrix(ndeg)
    buf = np.full((2, ndeg + 1, ndeg + 1, ndeg + 1), 7.5)
    sub = interp_element(buf, H)
    assert np.allclose(sub, 7.5), "constant field not preserved"

    # 4) Name mappings are complete and unique.
    for ncomp in (7, 22):
        nm = ordered_names(ncomp)
        assert len(set(nm)) == ncomp
    print("cstar2h5 selftest: OK")


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description="Convert homo Cstar (-82) to SEM anisotropic HDF5")
    p.add_argument("cstar", nargs="?", help="input Cstar binary (icode=-82)")
    p.add_argument("out", nargs="?", help="output HDF5 file")
    p.add_argument("--kind", choices=("auto", "elastic", "acoustic"),
                   default="auto", help="expected material kind (default: auto)")
    p.add_argument("--split", action="store_true",
                   help="write one file per property instead of one grouped file")
    p.add_argument("--selftest", action="store_true",
                   help="run internal checks and exit")
    args = p.parse_args(argv)

    if args.selftest:
        _selftest()
        return 0

    if not args.cstar or not args.out:
        p.error("cstar and out are required (unless --selftest)")

    header, fields = read_cstar(args.cstar)

    expected = {"elastic": 22, "acoustic": 7}
    if args.kind != "auto" and header.ncomp != expected[args.kind]:
        p.error("file has ncomp=%d but --kind %s expects %d"
                % (header.ncomp, args.kind, expected[args.kind]))

    kind = "elastic" if header.ncomp == 22 else "acoustic"
    print("Cstar: kind=%s ncomp=%d ndeg=%d nel=(%d,%d,%d) el=(%.4g,%.4g,%.4g)"
          % (kind, header.ncomp, header.ndeg, header.nelx, header.nely,
             header.nelz, header.xel, header.yel, header.zel))
    print("Output uniform grid: %d x %d x %d"
          % (header.nelx * header.ndeg + 1,
             header.nely * header.ndeg + 1,
             header.nelz * header.ndeg + 1))

    write_h5(header, fields, args.out, split=args.split)
    print("Wrote %s%s" % (args.out, " (split)" if args.split else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
