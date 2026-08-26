#!/usr/bin/env python3
"""Derivative scaling: g = dR/dt |_{t_c} ~ L^{y_t}.

Extracts y_t from the scaling of a *derivative estimator*, not by
differentiating a fitted R(t) curve.
"""

import _path  # noqa: F401

import argparse

import numpy as np

import fss


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file")
    ap.add_argument("--obs", default="obs")
    ap.add_argument("--at", type=float, required=True, help="critical control value (e.g. t_c or p_c)")
    ap.add_argument("--lmin", type=float)
    ap.add_argument("--lmax", type=float)
    args = ap.parse_args()

    data = fss.FSSData.load(args.file)
    size, g, gerr = fss.derivative.control_derivative(data, args.obs, args.at)
    if len(size) < 3:
        raise SystemExit("fewer than 3 sizes bracket the control value; cannot fit derivative scaling")
    print("derivative estimator g = dR/dt at t =", args.at)
    print("  L        g          g_err")
    for s, gv, ge in zip(size, g, gerr):
        print(f"  {s:<8g} {gv:<10.4g} {ge:<8.3g}")

    res = fss.derivative.derivative_scaling_fit(size, g, gerr, lmin=args.lmin, lmax=args.lmax)
    print("\n" + res.summary())

    # cross-check: effective exponent drift of g
    eff = fss.diagnostics.effective_exponent_series(size, g, gerr)
    print("\neffective exponent of g between consecutive sizes:")
    for mid, y, ye in eff:
        print(f"  L_mid={mid:<8g} y_eff={y:<8.4g} +/- {ye:.4g}")
    y_eff = eff[:, 1]
    if len(y_eff) > 2:
        drift = abs(y_eff[-1] - y_eff[0])
        print(f"\ndrift across sizes = {drift:.3f}"
              + ("  (OK)" if drift < 0.3 else "  (large - check ansatz / L_min)"))


if __name__ == "__main__":
    main()
