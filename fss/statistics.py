"""Statistical quantities for fit results, model comparison, and diagnostics.

Conventions used across the toolkit:

- ``chi2`` is the weighted sum of squared residuals, ``sum((obs-pred)^2/sigma^2)``.
- Information criteria are computed on the chi-square scale (the Gaussian
  constant term is dropped, as is standard in FSS practice):
    AIC  = chi2 + 2*k
    AICc = AIC + 2*k*(k+1)/(n - k - 1)
    BIC  = chi2 + k*ln(n)
  They are therefore meaningful only for comparing fits of the same data
  with correctly specified errors, which is exactly the intended use.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy import stats

__all__ = [
    "chi2_p_value",
    "aic",
    "aicc",
    "bic",
    "correlation_from_cov",
    "condition_number",
    "runs_test",
    "effective_exponent_pair",
    "FitReport",
]


def chi2_p_value(chi2: float, dof: int) -> float:
    """Right-tail probability of ``chi2`` with ``dof`` degrees of freedom."""
    if dof <= 0:
        return float("nan")
    return float(stats.chi2.sf(chi2, dof))


def aic(chi2: float, k: int) -> float:
    return float(chi2 + 2.0 * k)


def aicc(chi2: float, k: int, n: int) -> float:
    if n - k - 1 <= 0:
        return float("nan")
    return float(aic(chi2, k) + 2.0 * k * (k + 1) / (n - k - 1))


def bic(chi2: float, k: int, n: int) -> float:
    if n <= 0:
        return float("nan")
    return float(chi2 + k * np.log(n))


def correlation_from_cov(cov: np.ndarray) -> np.ndarray:
    """Pearson correlation matrix from a covariance matrix."""
    cov = np.asarray(cov, dtype=float)
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1]:
        raise ValueError("covariance must be a square matrix")
    d = np.sqrt(np.maximum(np.diag(cov), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.outer(d, d)
    corr[np.isnan(corr)] = 0.0
    corr[np.isinf(corr)] = 0.0
    return corr


def condition_number(cov: np.ndarray) -> float:
    """Condition number of a covariance matrix (>= 1)."""
    cov = np.asarray(cov, dtype=float)
    if cov.size == 0:
        return float("nan")
    if cov.shape[0] == 1:
        return 1.0 if cov[0, 0] > 0 else float("nan")
    w = np.linalg.eigvalsh(cov)
    w = w[w > 0]
    if len(w) == 0:
        return float("inf")
    return float(w.max() / w.min())


def runs_test(residuals: np.ndarray) -> tuple:
    """Runs (Wald-Wolfowitz) test for randomness of a signed sequence.

    Returns (z, p) where ``p`` is the two-sided probability under the
    null of a random sign sequence.  A small ``p`` flags structure
    (e.g. a remaining size trend in the residuals).
    """
    r = np.asarray(residuals, dtype=float)
    r = r[np.isfinite(r)]
    if len(r) < 3:
        return float("nan"), float("nan")
    signs = np.sign(r)
    signs = signs[signs != 0]
    n = len(signs)
    if n == 0:
        return float("nan"), float("nan")
    n_pos = int((signs > 0).sum())
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan"), 0.0
    runs = 1 + int(np.count_nonzero(np.diff(signs)))
    mu = 1 + 2 * n_pos * n_neg / n
    var = (2 * n_pos * n_neg * (2 * n_pos * n_neg - n)) / (n * n * (n - 1))
    if var <= 0:
        return float("nan"), float("nan")
    z = (runs - mu) / np.sqrt(var)
    p = 2.0 * stats.norm.sf(abs(z))
    return float(z), float(p)


def effective_exponent_pair(o1: float, o2: float, l1: float, l2: float,
                            s1: Optional[float] = None, s2: Optional[float] = None,
                            cov12: Optional[float] = None) -> tuple:
    """Local exponent between two sizes: y = ln(O2/O1) / ln(L2/L1).

    Returns (y, yerr).  Error propagates the reported uncertainties
    ``s1``/``s2`` and an optional covariance ``cov12``.  When the
    uncertainties are absent, ``yerr`` is NaN.
    """
    if l2 <= 0 or l1 <= 0 or l2 == l1:
        return float("nan"), float("nan")
    o1, o2 = float(o1), float(o2)
    if o1 <= 0 or o2 <= 0:
        return float("nan"), float("nan")
    dl = np.log(l2 / l1)
    y = np.log(o2 / o1) / dl
    if s1 is None or s2 is None:
        return y, float("nan")
    s1, s2 = float(s1), float(s2)
    if s1 <= 0 or s2 <= 0:
        return y, float("nan")
    c = 0.0 if cov12 is None else float(cov12)
    var = (s1 / o1) ** 2 + (s2 / o2) ** 2 - 2.0 * c / (o1 * o2)
    var = max(var, 0.0)
    return y, float(np.sqrt(var) / abs(dl))


@dataclass
class FitReport:
    """Lightweight per-parameter report used by scans and comparisons."""

    params: dict
    stderr: dict
    chi2: float
    dof: int
    chi2_reduced: float
    p_value: float
    aic: float
    aicc: float
    bic: float
    n_points: int
    n_params: int
    warnings: list
