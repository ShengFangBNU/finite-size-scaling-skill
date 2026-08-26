#!/usr/bin/env python3
"""L_min and fitting-window stability scans.

The L_min rule (fitting-and-systematics.md): take the smallest L_min with
chi2/dof ~ 1 such that raising L_min further does not lower chi2 by much
more than about one unit per degree of freedom.

Modes:
  --kind lmin    scan the size cut at a fixed control value (critical
                 power law, optionally with leading correction);
  --kind window  scan nested control windows [-w, w] around a center for
                 the dimensionless near-critical ansatz.
"""

import _path  # noqa: F401

import argparse

import numpy as np

import fss


def _parse_list(text):
    return [float(x) for x in text.split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file")
    ap.add_argument("--obs", default="obs")
    ap.add_argument("--kind", choices=["lmin", "window", "both"], default="lmin")
    ap.add_argument("--control-value", type=float, help="fixed control value (lmin mode)")
    ap.add_argument("--correction", action="store_true", help="include leading correction (lmin mode)")
    ap.add_argument("--omega", type=float, default=1.0)
    ap.add_argument("--lmin-values", default="16,24,32,48,64",
                    help="comma list of L_min to scan")
    ap.add_argument("--lmax", type=float, help="upper size cut for lmin mode")
    ap.add_argument("--center", type=float, help="window center (window mode; default: control-range midpoint)")
    ap.add_argument("--half-widths", default="0.05,0.1,0.15,0.2",
                    help="comma list of half-widths for window mode")
    ap.add_argument("--degree", type=int, default=2, help="dimensionless polynomial degree")
    ap.add_argument("--lmin", type=float, help="fixed L_min for window mode")
    args = ap.parse_args()

    data = fss.FSSData.load(args.file)
    params = ["y"] if args.kind in ("lmin", "both") and not args.correction else ["y", "b"]

    if args.kind in ("lmin", "both"):
        sel = data
        if args.control_value is not None:
            sel = sel.select_control_value(args.control_value)
        else:
            vals = np.unique(sel.df[sel.control_col].to_numpy(dtype=float))
            if len(vals) == 1:
                sel = sel.select_control_value(vals[0])
            else:
                raise SystemExit("lmin mode needs --control-value or single-value data")
        _, size, obs, err = sel.xy(args.obs)

        def lmin_fn(lmin=None, lmax=None, **kw):
            return fss.fit_critical_power(
                size, obs, err, correction=args.correction, omega=args.omega,
                lmin=lmin, lmax=lmax,
            )

        rows = fss.diagnostics.lmin_scan(lmin_fn, _parse_list(args.lmin_values),
                                         lmax=args.lmax, param_names=params)
        print("L_min scan (critical power law):")
        print(fss.diagnostics.scan_table(rows))
        print("\nrule: smallest L_min with chi2/dof ~ 1 and no further drop "
              "of more than ~1/DF on increasing L_min.")

    if args.kind in ("window", "both"):
        control, size, obs, err = data.xy(args.obs)
        center = args.center
        if center is None:
            cmin, cmax = data.control_range()
            center = 0.5 * (cmin + cmax)
            print(f"\n(window center from control range midpoint: {center})")
        else:
            print()

        def win_fn(control_window=None, **kw):
            return fss.fit_dimensionless(
                control, size, obs, err, degree=args.degree, with_pc=True,
                control_window=control_window, lmin=args.lmin,
            )

        rows = fss.diagnostics.window_scan(
            win_fn, center, _parse_list(args.half_widths),
            param_names=["pc", "yt"],
        )
        print("control-window scan (dimensionless near-critical fit, center =", center, "):")
        print(fss.diagnostics.scan_table(rows))
        print("\nflag any window where pc or yt drifts by more than its error bar.")


if __name__ == "__main__":
    main()
