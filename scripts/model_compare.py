#!/usr/bin/env python3
"""Model comparison for competing ansatze on the same data.

Candidates (size-axis fits at a fixed control value):
  1. O = a L^y                          (pure power law)
  2. O = a L^y + b L^(y-omega)          (+ leading correction)
  3. O = c0 + a0 L^y                    (+ analytic background)
  4. O = c0 + a0 L^y + b L^(y-omega)    (+ correction and background)

Comparison uses chi2/dof plus AIC/AICc/BIC on the chi-square scale.  Add
a term only if chi2/dof improves meaningfully and the parameters stay
stable (fitting-and-systematics.md: prefer the fit whose parameters are
stable and whose physics is justified, not merely the smallest error).
"""

import _path  # noqa: F401

import argparse

import numpy as np

import fss
from fss.models import critical_power_spec, critical_power_correction_spec, scaling_observable_spec


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

    control, size, obs, err = data.xy(args.obs)
    if not np.isfinite(err).any():
        err = None

    specs = [
        ("pure", critical_power_spec()),
        ("+ correction", critical_power_correction_spec(
            omega=args.omega, omega_fixed=not args.omega_free)),
        ("+ background", scaling_observable_spec(background=True)),
        ("+ correction + background", scaling_observable_spec(
            correction_exponents=[(args.omega, "fixed")], background=True)),
    ]
    results = []
    print(f"{'ansatz':24s} {'chi2':>8s} {'dof':>4s} {'chi2/dof':>9s} {'AIC':>8s} {'AICc':>9s} {'BIC':>8s} {'y':>8s}")
    for name, spec in specs:
        res = fss.fit_spec(spec, np.zeros_like(size), size, obs, err,
                           lmin=args.lmin, lmax=args.lmax)
        results.append(res)
        y = res.value("y")
        yerr = res.error("y")
        print(f"{name:24s} {res.chi2:8.3f} {res.dof:4d} {res.chi2_reduced:9.3f} "
              f"{res.aic:8.2f} {res.aicc:9.2f} {res.bic:8.2f} {y:8.4g} +/- {yerr:.3g}")

    best = min(results, key=lambda r: r.aicc)
    print(f"\nbest by AICc: {best.model_name!r} ({best.describe})")
    print("note: prefer the ansatz with stable parameters and justified physics,"
          "\n      not merely the smallest nominal error.")


if __name__ == "__main__":
    main()
