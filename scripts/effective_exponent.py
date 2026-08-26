#!/usr/bin/env python3
"""Effective-exponent diagnostics: y_eff(L) = d ln O / d ln L.

Reports the local exponent between consecutive sizes (or between L and
ratio*L), which reveals corrections to scaling through a systematic drift.
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
    ap.add_argument("--ratio", type=float, default=2.0,
                    help="pairing ratio; use --ratio 0 for consecutive pairs")
    args = ap.parse_args()

    data = fss.FSSData.load(args.file)
    if args.control_value is not None:
        data = data.select_control_value(args.control_value)
    else:
        vals = np.unique(data.df[data.control_col].to_numpy(dtype=float))
        if len(vals) == 1:
            data = data.select_control_value(vals[0])

    control, size, obs, err = data.xy(args.obs)
    if not np.isfinite(err).any():
        err = None

    if args.ratio and args.ratio > 1:
        series = fss.diagnostics.effective_exponent_ratios(size, obs, err, ratio=args.ratio)
        label = f"L -> {args.ratio} L"
    else:
        series = fss.diagnostics.effective_exponent_series(size, obs, err)
        label = "consecutive sizes"
    if len(series) == 0:
        raise SystemExit("not enough sizes for an effective-exponent series")

    print(f"effective exponent, pairing = {label}")
    print("  L_mid    y_eff     y_err")
    for mid, y, ye in series:
        print(f"  {mid:<9g} {y:<9.5g} {ye:.5g}")

    ys = series[:, 1]
    drift = ys[-1] - ys[0] if len(ys) > 1 else float("nan")
    print(f"\ndrift y_eff(L_max) - y_eff(L_min) = {drift:.4g}")
    if len(ys) > 2 and np.isfinite(drift):
        if abs(drift) > 0.2:
            print("  large drift: corrections to scaling are significant;")
            print("  do not quote a single pure-power exponent without them")
        else:
            print("  small drift: consistent with a stable exponent over this size range")


if __name__ == "__main__":
    main()
