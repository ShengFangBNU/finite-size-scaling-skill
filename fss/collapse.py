"""Ordinary data-collapse quality metric.

Collapse is a *diagnostic*, not the estimator itself (the toolkit uses
fits on the scaling fields as the primary estimators; see
`references/scaling-field-framework.md`).  This module provides a single
deterministic number for a hypothesized collapse: the reduced chi2 of the
vertical spread of the collapsed points about a locally smooth master
curve, obtained by binning in the scaling variable x.

If the data collapse onto one curve and the error bars are correct, the
reduced chi2 is ~1.  Larger values flag a bad collapse / wrong exponents;
smaller values flag overestimated errors.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

__all__ = ["collapse_quality"]


def collapse_quality(x, y, err, n_bins: int = 8):
    """Binned vertical-spread chi2 for a candidate collapse.

    Parameters
    ----------
    x : array-like
        scaling variable (e.g. ``t L^yt``).
    y : array-like
        collapsed observable (e.g. ``R`` or ``O L^-yO``).
    err : array-like
        uncertainty of ``y`` (same shape).
    n_bins : int
        number of equal-count bins in x.

    Returns a dict with ``chi2``, ``dof`` (points minus bins), reduced
    chi2, and the per-bin point counts.  Points with non-finite x/y/err
    are dropped.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    err = np.asarray(err, float)
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(err) & (err > 0)
    x, y, err = x[keep], y[keep], err[keep]
    if len(x) < 2 or n_bins < 1:
        return {"chi2": float("nan"), "dof": 0, "chi2_reduced": float("nan"),
                "n_points": int(len(x)), "bin_counts": []}
    n_bins = min(n_bins, len(x))
    order = np.argsort(x)
    x, y, err = x[order], y[order], err[order]

    edges = np.round(np.linspace(0, len(x), n_bins + 1)).astype(int)
    edges[-1] = len(x)
    chi2 = 0.0
    counts = []
    for b in range(n_bins):
        lo, hi = edges[b], edges[b + 1]
        n = hi - lo
        if n == 0:
            counts.append(0)
            continue
        counts.append(int(n))
        w = 1.0 / err[lo:hi] ** 2
        ybar = float(np.sum(w * y[lo:hi]) / np.sum(w))
        chi2 += float(np.sum(((y[lo:hi] - ybar) / err[lo:hi]) ** 2))
    dof = max(len(x) - n_bins, 1)
    return {
        "chi2": chi2,
        "dof": int(dof),
        "chi2_reduced": chi2 / dof,
        "n_points": int(len(x)),
        "bin_counts": counts,
    }
