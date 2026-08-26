#!/usr/bin/env python3
"""Pairwise crossing analysis for a dimensionless observable.

Finds the control value t_x(L, sL) where R(t, L) curves cross, fits the
trajectory t_x = t_c + a L^(-lambda), and *checks the linear amplitude*
a_1 before reporting t_c (an observable with a_1 = 0 does not give a
valid crossing estimate of the critical point).
"""

import _path  # noqa: F401

import argparse

import numpy as np

import fss


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file")
    ap.add_argument("--obs", default="obs")
    ap.add_argument("--control-window", type=float, nargs=2, metavar=("LO", "HI"))
    ap.add_argument("--no-fit", action="store_true", help="only list crossings, skip the fit")
    ap.add_argument("--lmin", type=float, help="L_min for the crossing fit")
    ap.add_argument("--lmax", type=float, help="L_max for the crossing fit")
    args = ap.parse_args()

    data = fss.FSSData.load(args.file)
    win = tuple(args.control_window) if args.control_window else None
    rows = fss.crossings(data, args.obs, control_window=win)
    if not rows:
        raise SystemExit("no crossings found in the given window")
    print("pairwise crossings (t_x = control where R(t,L1) = R(t,L2)):")
    print("  L1       L2       t_x        t_x_err    flat")
    for r in rows:
        print(f"  {r['L1']:<8g} {r['L2']:<8g} {r['t_x']:<10.5g} {r['t_x_err']:<10.4g} {r['flat']}")
    usable = [r for r in rows if not r["flat"]]
    if not usable:
        raise SystemExit("all crossings are flat; t_c cannot be estimated from these curves")

    sizes = np.unique([r["L1"] for r in usable])
    largest = sizes.max()
    amp = fss.crossing.linear_amplitude_check(data, args.obs, largest, win)
    print(f"\nlinear-amplitude check at L={largest} (a1 must be != 0 for crossings to estimate t_c):")
    print(f"  a1 = {amp.value('a1'):.5g} +/- {amp.error('a1'):.5g}"
          f"   (a1/sigma = {amp.value('a1') / max(amp.error('a1'), 1e-300):.2f})")

    if args.no_fit:
        return
    res = fss.crossing.crossing_fit(usable, lmin=args.lmin, lmax=args.lmax)
    print("\n" + res.summary())
    print("\ninterpretation: t_c =", f"{res.value('t_c'):.6g}", "+/-", f"{res.error('t_c'):.4g}",
          "  lambda =", f"{res.value('lam'):.4g}", "(~ y_t + |y_i| expected)")


if __name__ == "__main__":
    main()
