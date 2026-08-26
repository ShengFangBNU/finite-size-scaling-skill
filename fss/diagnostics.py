"""Statistical diagnostics: effective exponents, residual checks, stability scans.

Implements the `L_min`-selection rule, the fitting-window scan, residual
diagnostics, and the identifiability report described in
`references/fitting-and-systematics.md` and the paper notes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

from . import statistics as _stats
from .fitting import FitResult

__all__ = [
    "effective_exponent_series",
    "effective_exponent_ratios",
    "residual_diagnostics",
    "identifiability_report",
    "lmin_scan",
    "window_scan",
    "scan_table",
]


# ----------------------------------------------------------------------
# effective exponents
# ----------------------------------------------------------------------

def effective_exponent_series(size, obs, err=None):
    """Local exponent between consecutive sizes: y(L) = d ln O / d ln L.

    Returns an (n-1, 3) array of ``(L_mid, y, yerr)``, where ``L_mid`` is
    the geometric mean of each pair.  Errors are NaN when the data carry
    no uncertainties.
    """
    size = np.asarray(size, float)
    obs = np.asarray(obs, float)
    if err is not None:
        err = np.asarray(err, float)
    order = np.argsort(size)
    size, obs = size[order], obs[order]
    if err is not None:
        err = err[order]
    out = []
    for i in range(len(size) - 1):
        l1, l2 = size[i], size[i + 1]
        y, ye = _stats.effective_exponent_pair(
            obs[i], obs[i + 1], l1, l2,
            s1=err[i] if err is not None else None,
            s2=err[i + 1] if err is not None else None,
        )
        out.append((float(np.sqrt(l1 * l2)), y, ye))
    return np.array(out)


def effective_exponent_ratios(size, obs, err=None, ratio=2):
    """Local exponent between L and ``ratio*L``: y = ln O(sL)/ln s.

    Only pairs whose ``ratio*L`` is present in the data are reported.
    """
    size = np.asarray(size, float)
    obs = np.asarray(obs, float)
    if err is not None:
        err = np.asarray(err, float)
    lut = {round(l, 10): (o, e) for l, o, e in zip(
        size, obs, err if err is not None else [None] * len(size))}
    out = []
    for i, l in enumerate(size):
        sl = ratio * l
        hit = lut.get(round(sl, 10))
        if hit is None:
            continue
        o1, e1 = obs[i], (err[i] if err is not None else None)
        o2, e2 = hit
        y, ye = _stats.effective_exponent_pair(o1, o2, l, sl, s1=e1, s2=e2)
        out.append((float(np.sqrt(l * sl)), y, ye))
    out.sort(key=lambda r: r[0])
    return np.array(out)


# ----------------------------------------------------------------------
# residual diagnostics
# ----------------------------------------------------------------------

def _pearson(x, y):
    """Pearson correlation, NaN (no warning) when either input has zero variance."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    x = x - x.mean()
    y = y - y.mean()
    sx = np.sqrt(np.sum(x * x))
    sy = np.sqrt(np.sum(y * y))
    if sx == 0.0 or sy == 0.0:
        return float("nan")
    return float(np.dot(x, y) / (sx * sy))


def residual_diagnostics(result: FitResult):
    """Check the residuals of a fit for remaining structure.

    Returns a dict with:
      - ``runs_z``, ``runs_p``: Wald-Wolfowitz runs test on the sign of the
        unweighted residuals ordered by size (small p = structured residuals);
      - ``corr_with_logsize``, ``corr_with_control``: Pearson correlations of
        the unweighted residuals with log L and with the control variable.
        (NaN when a variable has zero variance -- e.g. residuals that are
        exactly constant across the data.)
    """
    out = {}
    size = result.size
    if size is not None and len(size) > 2:
        order = np.argsort(size)
        r = result.residuals_unweighted[order]
        z, p = _stats.runs_test(r)
        out["runs_z"] = z
        out["runs_p"] = p
        out["corr_with_logsize"] = _pearson(np.log(size[order]), r)
    else:
        out["runs_z"] = out["runs_p"] = float("nan")
        out["corr_with_logsize"] = float("nan")
    control = result.control
    if control is not None and len(control) > 2:
        out["corr_with_control"] = _pearson(control, result.residuals_unweighted)
    else:
        out["corr_with_control"] = float("nan")
    return out


def identifiability_report(result: FitResult) -> dict:
    """Condensed identifiability diagnostics for one fit."""
    corr = result.corr
    pairs = []
    if corr.ndim == 2 and corr.shape[0] > 1:
        names = result.param_order
        for i in range(corr.shape[0]):
            for j in range(i + 1, corr.shape[1]):
                c = float(corr[i, j])
                if abs(c) > 0.9:
                    pairs.append((names[i], names[j], c))
    return {
        "condition_number": result.condition_number,
        "max_abs_corr": result.max_abs_corr,
        "high_corr_pairs": pairs,
        "n_warnings": len(result.warnings),
        "warnings": list(result.warnings),
    }


# ----------------------------------------------------------------------
# stability scans
# ----------------------------------------------------------------------

def lmin_scan(fit_fn: Callable[..., FitResult], lmin_values, lmax=None,
              param_names: Optional[Sequence[str]] = None):
    """Run ``fit_fn(lmin=..., lmax=...)`` for each candidate L_min.

    ``fit_fn`` must accept keyword arguments ``lmin`` and ``lmax`` and
    return a :class:`FitResult`.  Returns a list of row dicts, one per
    L_min, with keys ``lmin``, ``chi2``, ``dof``, ``chi2_reduced``,
    ``p_value``, ``n_points`` and one ``<name>`` / ``<name>_err`` entry
    per parameter of interest.
    """
    rows = []
    for lmin in lmin_values:
        res = fit_fn(lmin=lmin, lmax=lmax)
        row = {
            "lmin": float(lmin),
            "chi2": res.chi2,
            "dof": res.dof,
            "chi2_reduced": res.chi2_reduced,
            "p_value": res.p_value,
            "n_points": res.n_points,
            "model": res.model_name,
        }
        names = param_names or res.param_order
        for name in names:
            row[name] = res.params.get(name, float("nan"))
            row[f"{name}_err"] = res.stderr.get(name, float("nan"))
        rows.append(row)
    return rows


def window_scan(fit_fn: Callable[..., FitResult], center: float, half_widths,
                param_names: Optional[Sequence[str]] = None):
    """Fit over nested control windows ``[center-w, center+w]``.

    ``fit_fn`` must accept keyword arguments ``control_window=(lo, hi)``
    and return a :class:`FitResult`.  Returns row dicts keyed by
    ``half_width`` plus the same statistics and parameters as
    :func:`lmin_scan`.
    """
    rows = []
    for w in half_widths:
        res = fit_fn(control_window=(center - w, center + w))
        row = {
            "half_width": float(w),
            "chi2": res.chi2,
            "dof": res.dof,
            "chi2_reduced": res.chi2_reduced,
            "p_value": res.p_value,
            "n_points": res.n_points,
            "model": res.model_name,
        }
        names = param_names or res.param_order
        for name in names:
            row[name] = res.params.get(name, float("nan"))
            row[f"{name}_err"] = res.stderr.get(name, float("nan"))
        rows.append(row)
    return rows


def scan_table(rows, headers=None) -> str:
    """Format scan rows as a readable aligned table (no pandas required)."""
    if not rows:
        return "(empty scan)"
    keys = headers or list(rows[0].keys())
    widths = {k: max(len(k), *(len(f"{r[k]:.4g}") if isinstance(r[k], (int, float)) and np.isfinite(r[k]) else len(str(r[k])) for r in rows)) for k in keys}
    lines = ["  ".join(k.ljust(widths[k]) for k in keys)]
    lines.append("  ".join("-" * widths[k] for k in keys))
    for r in rows:
        cells = []
        for k in keys:
            v = r[k]
            if isinstance(v, (int, float)) and np.isfinite(v):
                cells.append(f"{v:.4g}".ljust(widths[k]))
            else:
                cells.append(str(v).ljust(widths[k]))
        lines.append("  ".join(cells))
    return "\n".join(lines)
