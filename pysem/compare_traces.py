# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
compare_traces.py - Compare two SEM receiver traces (rec_*.vel text columns).

Reports, for a chosen column: finiteness, peak amplitude, waveform correlation,
scale ratio, and an optional TIME-INTEGRAL relationship test. The integral test
validates e.g. the fluid pressure source (i_type_source=7) against the fluidpulse
source (type 3): in SEM the type-7 response equals the running time-integral of the
type-3 response (source injected as int(f dt) vs f(t)), so
    corr( cumint(trace_A) , trace_B ) ~ 1  and  peak-ratio ~ 1.

    Ex.1 : compare column 1 (VelPhi = -pressure) of two runs
        python3 compare_traces.py runA/rec_0000.vel runB/rec_0000.vel --col 1

    Ex.2 : test that B is the time-integral of A (type-3 vs type-7 validation)
        python3 compare_traces.py type3/rec_0000.vel type7/rec_0000.vel --col 1 --integral

    Ex.3 : self-test (no data needed)
        python3 compare_traces.py --selftest
"""
import argparse
import sys

import numpy as np


def load_trace(path):
    """Return (t, data2d) from a rec_*.vel text file (col 0 = time)."""
    d = np.loadtxt(path)
    if d.ndim == 1:
        d = d[:, None]
    t = d[:, 0]
    return t, d


def _norm(x):
    m = np.max(np.abs(x))
    return x / m if m > 0 else x


def cumtrapz0(y, t):
    """Cumulative trapezoidal integral of y over t, starting at 0 (same length as y)."""
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])


def compare(tA, yA, tB, yB, integral=False):
    """Compare two 1-D signals. Returns a dict of metrics."""
    n = min(len(tA), len(tB))
    tA, yA, yB = tA[:n], yA[:n], yB[:n]
    res = {
        'n': n,
        'finite_A': bool(np.all(np.isfinite(yA))),
        'finite_B': bool(np.all(np.isfinite(yB))),
        'peak_A': float(np.max(np.abs(yA))),
        'peak_B': float(np.max(np.abs(yB))),
        'corr_direct': float(np.corrcoef(_norm(yA), _norm(yB))[0, 1]),
    }
    res['scale_B_over_A'] = res['peak_B'] / res['peak_A'] if res['peak_A'] else float('nan')
    if integral:
        IA = cumtrapz0(yA, tA)
        res['peak_int_A'] = float(np.max(np.abs(IA)))
        res['corr_integral'] = float(np.corrcoef(_norm(IA), _norm(yB))[0, 1])
        res['scale_B_over_intA'] = (res['peak_B'] / res['peak_int_A']
                                    if res['peak_int_A'] else float('nan'))
    return res


def print_report(res, integral=False):
    print("samples compared : {}".format(res['n']))
    print("finite A / B     : {} / {}".format(res['finite_A'], res['finite_B']))
    print("peak |A| / |B|   : {:.4g} / {:.4g}".format(res['peak_A'], res['peak_B']))
    print("scale |B|/|A|    : {:.4g}".format(res['scale_B_over_A']))
    print("corr(A, B)       : {:.5f}".format(res['corr_direct']))
    if integral:
        print("peak |cumint A|  : {:.4g}".format(res['peak_int_A']))
        print("scale |B|/|intA| : {:.4g}".format(res['scale_B_over_intA']))
        print("corr(cumint A,B) : {:.5f}   <- ~1.0 means B is the time-integral of A"
              .format(res['corr_integral']))


def _selftest():
    t = np.linspace(0, 1, 2000)
    # A = a Ricker-ish pulse; B = its exact cumulative integral
    f0 = 8.0
    a = (1 - 2 * (np.pi * f0 * (t - 0.4)) ** 2) * np.exp(-(np.pi * f0 * (t - 0.4)) ** 2)
    b = cumtrapz0(a, t)
    res = compare(t, a, t, b, integral=True)
    assert res['corr_integral'] > 0.999, res['corr_integral']
    assert abs(res['scale_B_over_intA'] - 1.0) < 1e-6, res['scale_B_over_intA']
    assert abs(res['corr_direct']) < 0.5, res['corr_direct']  # A and its integral are ~uncorrelated
    print("compare_traces selftest OK (corr_integral={:.5f})".format(res['corr_integral']))


def main():
    p = argparse.ArgumentParser(description="Compare two SEM receiver traces.")
    p.add_argument('traceA', nargs='?', help="first rec_*.vel")
    p.add_argument('traceB', nargs='?', help="second rec_*.vel")
    p.add_argument('--col', type=int, default=1, help="column to compare (default 1)")
    p.add_argument('--integral', action='store_true',
                   help="test whether B == time-integral of A (type-3 vs type-7)")
    p.add_argument('--selftest', action='store_true', help="run internal checks and exit")
    a = p.parse_args()
    if a.selftest:
        _selftest(); return
    if not a.traceA or not a.traceB:
        p.error("need traceA and traceB (or --selftest)")
    tA, dA = load_trace(a.traceA)
    tB, dB = load_trace(a.traceB)
    c = a.col
    if c >= dA.shape[1] or c >= dB.shape[1]:
        print("ERR: column {} out of range (A has {}, B has {})".format(c, dA.shape[1], dB.shape[1]))
        sys.exit(2)
    res = compare(tA, dA[:, c], tB, dB[:, c], integral=a.integral)
    print_report(res, integral=a.integral)


if __name__ == '__main__':
    main()
