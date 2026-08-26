"""Derivative estimators of the thermal scaling dimension.

Two estimators of the same object (a quantity scaling as ``L^{y_t}`` at
the critical point):

1. the continuous derivative of a dimensionless observable,
   ``g(t) = dR/dt``, evaluated at ``t = t_c``; and
2. the (discrete) covariance estimator used in bond percolation,
   ``g(p) = cov(R, N_b)`` (proportional to ``dR/dp``; Wang2013, Eq. 4).

The exponent ``y_t`` is extracted by fitting ``g(t_c, L) ~ a L^{y_t}``;
it is *not* obtained by differentiating a fitted R(p) curve
(see `references/paper-notes/01-3d-percolation.md`).
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from . import fitting as _fit
from .io import FSSData

__all__ = ["control_derivative", "covariance_estimator", "derivative_scaling_fit"]


def control_derivative(data: FSSData, observable: str, at: float):
    """Estimate ``g = dR/dt |_{t=at}`` at every size.

    At a control grid point the symmetric difference
    ``(R(t+h) - R(t-h)) / 2h`` is used (so the derivative is evaluated at
    ``at`` itself, not at a segment midpoint); at the edge of the control
    range a one-sided difference is used.  Between grid points the slope is
    taken from the parabola through the three bracketing points, evaluated
    at ``at`` (a plain secant would return the slope at the *segment
    midpoint*, biasing the estimate when the curve is curved there).
    Returns ``(size, g, gerr)`` arrays, one point per size whose control
    range brackets ``at``.  The uncertainty propagates the reported point
    errors through the estimator.
    """
    control, size, obs, err = data.xy(observable)
    out_s, out_g, out_ge = [], [], []
    for L in np.unique(size):
        m = size == L
        c, o, e = control[m], obs[m], err[m]
        order = np.argsort(c)
        c, o, e = c[order], o[order], e[order]
        if at < c[0] or at > c[-1]:
            continue
        i = int(np.argmin(np.abs(c - at)))
        if np.isclose(c[i], at):
            if 0 < i < len(c) - 1:
                # symmetric central difference at the grid point
                dc = c[i + 1] - c[i - 1]
                g = (o[i + 1] - o[i - 1]) / dc
                ge = float(np.hypot(e[i + 1], e[i - 1]) / abs(dc))
            else:
                # boundary grid point: one-sided difference
                j = i + 1 if i == 0 else i - 1
                dc = c[j] - c[i]
                if dc <= 0:
                    continue
                g = (o[j] - o[i]) / dc
                ge = float(np.hypot(e[i], e[j]) / abs(dc))
        elif len(c) >= 3:
            # between grid points: parabola through the three bracketing
            # points, slope evaluated at `at`
            i_hi = int(np.searchsorted(c, at))   # first c >= at
            i_lo = i_hi - 1
            if i_lo >= 1 and i_hi <= len(c) - 2:
                win = slice(i_lo - 1, i_hi + 1)
            elif i_lo == 0:
                win = slice(0, 3)
            else:
                win = slice(len(c) - 3, len(c))
            g, ge = _quad_deriv(c[win], o[win], e[win], at)
        else:
            # only two points: fall back to the secant across them
            i_hi = int(np.searchsorted(c, at))
            i_lo = i_hi - 1
            dc = c[i_hi] - c[i_lo]
            if dc <= 0:
                continue
            g = (o[i_hi] - o[i_lo]) / dc
            ge = float(np.hypot(e[i_lo], e[i_hi]) / abs(dc))
        out_s.append(float(L))
        out_g.append(float(g))
        out_ge.append(ge)
    return np.array(out_s), np.array(out_g), np.array(out_ge)


def _quad_deriv(x, y, e, at):
    """Slope at ``at`` of the parabola through three consecutive points.

    The control derivative of a smooth observable is better estimated by the
    parabola through the three points bracketing ``at`` than by the secant
    across the containing segment (which returns the slope at the *segment
    midpoint*, not at ``at``).  Errors propagate linearly through the
    Lagrange basis of the interpolation.
    """
    x, y, e = (np.asarray(a, float) for a in (x, y, e))
    x0, x1, x2 = x
    y0, y1, y2 = y
    e0, e1, e2 = e
    Lp0 = (2 * at - x1 - x2) / ((x0 - x1) * (x0 - x2))
    Lp1 = (2 * at - x0 - x2) / ((x1 - x0) * (x1 - x2))
    Lp2 = (2 * at - x0 - x1) / ((x2 - x0) * (x2 - x1))
    g = float(Lp0 * y0 + Lp1 * y1 + Lp2 * y2)
    ge = float(np.hypot(Lp0 * e0, np.hypot(Lp1 * e1, Lp2 * e2)))
    return g, ge


def covariance_estimator(R, Nb):
    """``g = cov(R, N_b)`` over samples for one size, with its standard error.

    ``R`` and ``Nb`` are per-sample arrays (length = number of samples).
    This is the bond-percolation estimator ``g = p(1-p) dR/dp`` up to a
    constant factor (Wang2013).  Returns ``(g, gerr)``.
    """
    R = np.asarray(R, float)
    Nb = np.asarray(Nb, float)
    n = len(R)
    if n < 2:
        return float("nan"), float("nan")
    z = (R - R.mean()) * (Nb - Nb.mean())
    g = float(z.mean())
    ge = float(z.std(ddof=1) / np.sqrt(n))
    return g, ge


def derivative_scaling_fit(size, g, gerr=None, **fit_kwargs) -> _fit.FitResult:
    """Fit ``g(L) = a L^yt`` on the derivative estimator."""
    return _fit.fit_critical_power(size, g, gerr, y_init=1.0, **fit_kwargs)
