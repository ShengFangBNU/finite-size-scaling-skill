"""Minimal plotting helpers (matplotlib, Agg backend).

Kept deliberately small: the toolkit's outputs are numbers (fit results,
scan tables, diagnostics); plots are conveniences for the scripts.
All functions return the figure so scripts can save it.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

__all__ = [
    "plot_effective_exponents",
    "plot_lmin_trajectory",
    "plot_residuals",
    "plot_collapse",
]


def plot_effective_exponents(series, ax=None, yline=None, label=None, **kw):
    """``series`` is the (N, 3) output of
    :func:`fss.diagnostics.effective_exponent_series`."""
    if ax is None:
        _, ax = plt.subplots()
    series = np.asarray(series, float)
    errs = series[:, 2]
    ax.errorbar(series[:, 0], series[:, 1], yerr=errs, fmt="o", label=label, **kw)
    if yline is not None:
        ax.axhline(yline, color="tab:red", ls="--", lw=1, alpha=0.7)
    ax.set_xlabel(r"$L$")
    ax.set_ylabel(r"$y_\mathrm{eff}(L)$")
    ax.set_xscale("log")
    return ax


def plot_lmin_trajectory(rows, param, ax=None, **kw):
    """Trajectory of one parameter vs L_min from a :func:`lmin_scan`."""
    if ax is None:
        _, ax = plt.subplots()
    lmin = np.array([r["lmin"] for r in rows])
    val = np.array([r[param] for r in rows])
    err = np.array([r.get(f"{param}_err", np.nan) for r in rows])
    ax.errorbar(lmin, val, yerr=err, fmt="o", **kw)
    ax.set_xlabel(r"$L_\mathrm{min}$")
    ax.set_ylabel(param)
    ax.set_xscale("log")
    return ax


def plot_residuals(result, ax=None, **kw):
    """Unweighted residuals vs size (left) and, if present, vs control (right)."""
    if ax is None:
        _, ax = plt.subplots(1, 2, figsize=(10, 4))
    r = result.residuals_unweighted
    if result.size is not None:
        ax[0].axhline(0, color="grey", lw=1)
        ax[0].scatter(result.size, r, **kw)
        ax[0].set_xscale("log")
        ax[0].set_xlabel(r"$L$")
        ax[0].set_ylabel("unweighted residual")
    if result.control is not None and np.unique(result.control).size > 1:
        ax[1].axhline(0, color="grey", lw=1)
        ax[1].scatter(result.control, r, **kw)
        ax[1].set_xlabel("control")
        ax[1].set_ylabel("unweighted residual")
    return ax


def plot_collapse(x, y, err, ax=None, **kw):
    if ax is None:
        _, ax = plt.subplots()
    ax.errorbar(x, y, yerr=err, fmt="o", ms=4, **kw)
    ax.set_xlabel(r"$x = t\, L^{y_t}$")
    ax.set_ylabel(r"$R$")
    return ax
