#!/usr/bin/env python3
"""Critical power-law fit: O(L) = a L^y at a fixed control value.

The observable is taken at (or nearest to) a fixed control value such as
t = 0 or p = p_c.  With ``--correction`` the ansatz gains the leading
amplitude correction: O = a L^y + b L^(y-omega).
"""

import _path  # noqa: F401

import argparse

import numpy as np

import fss


def _select_fixed_control(data, value):
    if value is None:
        vals = np.unique(data.df[data.control_col].to_numpy(dtype=float))
        if len(vals) == 1:
            return data.select_control_value(vals[0])
        raise SystemExit(
            f"data has {len(vals)} control values; pass --control-value"
        )
    return data.select_control_value(value)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file")
    ap.add_argument("--obs", default="obs", help="observable name")
    ap.add_argument("--control-value", type=float, help="fixed control value (default: sole value)")
    ap.add_argument("--correction", action="store_true", help="add leading amplitude correction")
    ap.add_argument("--omega", type=float, default=1.0, help="leading irrelevant exponent (default 1.0)")
    ap.add_argument("--omega-free", action="store_true", help="leave omega free")
    ap.add_argument("--lmin", type=float)
    ap.add_argument("--lmax", type=float)
    args = ap.parse_args()

    data = fss.FSSData.load(args.file)
    data = _select_fixed_control(data, args.control_value)
    control, size, obs, err = data.xy(args.obs)
    if err is None or not np.isfinite(err).any():
        err = None

    res = fss.fit_critical_power(
        size, obs, err,
        correction=args.correction,
        omega=args.omega, omega_fixed=not args.omega_free,
        lmin=args.lmin, lmax=args.lmax,
    )
    print(res.summary())
    print()
    print("residual diagnostics:")
    for k, v in fss.diagnostics.residual_diagnostics(res).items():
        print(f"  {k:18s} = {v:.4g}")


if __name__ == "__main__":
    main()
