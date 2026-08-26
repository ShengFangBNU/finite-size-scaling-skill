#!/usr/bin/env python3
"""Leading-correction fit: O(L) = a L^y + b L^(y-omega).

omega defaults to 1.0 (pinned).  ``--omega-free`` lets it float; the
staged protocol is: fit free, check the value is consistent with theory,
then pin and refit (the pinned fit is what is reported).
"""

import _path  # noqa: F401

import argparse

import numpy as np

import fss


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file")
    ap.add_argument("--obs", default="obs")
    ap.add_argument("--control-value", type=float, help="fixed control value (default: sole value)")
    ap.add_argument("--omega", type=float, default=1.0)
    ap.add_argument("--omega-free", action="store_true")
    ap.add_argument("--lmin", type=float)
    ap.add_argument("--lmax", type=float)
    args = ap.parse_args()

    data = fss.FSSData.load(args.file)
    if args.control_value is not None:
        data = data.select_control_value(args.control_value)
    else:
        vals = np.unique(data.df[data.control_col].to_numpy(dtype=float))
        if len(vals) == 1:
            data = data.select_control_value(vals[0])
        else:
            raise SystemExit("multiple control values; pass --control-value")

    _, size, obs, err = data.xy(args.obs)

    if args.omega_free:
        free = fss.fit_critical_power(
            size, obs, err, correction=True,
            omega=args.omega, omega_fixed=False,
            lmin=args.lmin, lmax=args.lmax,
        )
        print("stage 1 (omega free) -- consistency check:")
        print(free.summary())
        w = free.value("omega")
        print(f"\n  omega = {w:.4g} +/- {free.error('omega'):.4g}")

    pinned = fss.fit_critical_power(
        size, obs, err, correction=True,
        omega=args.omega, omega_fixed=True,
        lmin=args.lmin, lmax=args.lmax,
    )
    print("\nstage 2 (omega pinned) -- reported fit:")
    print(pinned.summary())
    print("\ninterpretation: y =", f"{pinned.value('y'):.5g}",
          "+/-", f"{pinned.error('y'):.5g}",
          "   b/a (relative correction at L=1) =",
          f"{pinned.value('b') / pinned.value('a'):.4g}")


if __name__ == "__main__":
    main()
