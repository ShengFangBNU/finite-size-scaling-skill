"""Pairwise-crossing analysis for dimensionless observables.

A critical point is located from crossings only when the linear amplitude
``a_1`` of the observable is non-zero (``R(t,L) = Rc + a1 t L^yt + ...``);
an observable whose leading field coefficient vanishes has crossings that
*do not* converge to t_c.  Callers should fit the leading slope (see
:func:`linear_amplitude_check`) before trusting crossings.

For a genuine transition the crossing trajectory satisfies
``t_x(L, sL) = t_c + a L^{-lambda}`` with ``lambda ~ y_t + |y_i|``
(see `references/dimensionless-crossings.md`).
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from . import fitting as _fit
from .io import FSSData
from .models import dimensionless_near_critical_spec, crossing_spec

__all__ = ["crossings", "crossing_fit", "linear_amplitude_check"]


def _linear_interp(cs, vs, xs):
    cs = np.asarray(cs, float)
    vs = np.asarray(vs, float)
    xs = np.asarray(xs, float)
    idx = np.clip(np.searchsorted(cs, xs), 1, len(cs) - 1)
    c0, c1 = cs[idx - 1], cs[idx]
    v0, v1 = vs[idx - 1], vs[idx]
    denom = c1 - c0
    w = np.where(denom > 0, (xs - c0) / denom, 0.0)
    return v0 + w * (v1 - v0)


def _crossing_root(c1, y1, e1, c2, y2, e2, lo=None, hi=None, ref=None, n_grid=4001):
    """Intersection of two smooth curves, with propagated error.

    The difference ``D = R1 - R2`` is evaluated on a dense grid of
    ``n_grid`` points spanning the common control range (piecewise-linear
    interpolation between the two curves), because the crossing of two
    curved lines on a sparse grid is badly biased.  Returns
    ``(t_x, t_x_err)`` or ``(nan, nan)`` if the curves do not cross inside
    ``[lo, hi]``.  When several crossings exist (e.g. a quadratic-in-t
    observable), the one nearest ``ref`` is returned.
    """
    c1, c2 = np.asarray(c1, float), np.asarray(c2, float)
    cmin = max(c1.min(), c2.min())
    cmax = min(c1.max(), c2.max())
    if lo is not None:
        cmin = max(cmin, lo)
    if hi is not None:
        cmax = min(cmax, hi)
    if cmax <= cmin:
        return float("nan"), float("nan")
    grid = np.linspace(cmin, cmax, n_grid)
    d = _linear_interp(c1, y1, grid) - _linear_interp(c2, y2, grid)
    candidates = []

    def _in_window(t):
        return (lo is None or t >= lo) and (hi is None or t <= hi)

    def _add(t, terr):
        dist = abs(t - ref) if ref is not None else 0.0
        candidates.append((dist, t, terr))

    # exact zeros at grid nodes
    for i in np.where(d == 0.0)[0]:
        t = float(grid[i])
        if not _in_window(t):
            continue
        # a touch point (curves meet but do not cross) is not a valid
        # crossing: its location is ill-determined -> flat
        if 0 < i < len(d) - 1 and np.sign(d[i - 1]) == np.sign(d[i + 1]) \
                and d[i - 1] != 0.0 and d[i + 1] != 0.0:
            _add(t, float("nan"))
            continue
        eR1 = _linear_interp(c1, e1, [t])[0]
        eR2 = _linear_interp(c2, e2, [t])[0]
        slope = _local_slope(grid, d, t)
        terr = float(np.hypot(eR1, eR2) / abs(slope)) if abs(slope) > 1e-12 else float("nan")
        _add(t, terr)

    sgn = np.sign(d)
    crossing = np.where((np.diff(sgn) != 0) & (d[:-1] * d[1:] < 0))[0]
    for i in crossing:
        c0, c1_ = grid[i], grid[i + 1]
        d0, d1 = d[i], d[i + 1]
        t = c0 - d0 * (c1_ - c0) / (d1 - d0)
        if not _in_window(t):
            continue
        slope = (d1 - d0) / (c1_ - c0)
        if abs(slope) < 1e-12:
            continue
        eR1 = _linear_interp(c1, e1, [t])[0]
        eR2 = _linear_interp(c2, e2, [t])[0]
        _add(t, float(np.hypot(eR1, eR2) / abs(slope)))

    if not candidates:
        return float("nan"), float("nan")
    candidates.sort(key=lambda c: c[0])
    return candidates[0][1], candidates[0][2]


def _local_slope(grid, d, t):
    i = int(np.searchsorted(grid, t))
    i = min(max(i, 1), len(grid) - 1)
    if i < len(grid):
        return (d[i] - d[i - 1]) / (grid[i] - grid[i - 1])
    return 1.0


def _by_size(data: FSSData, observable: str, lo=None, hi=None) -> dict:
    """{L: (control_sorted, obs_sorted, err_sorted)} for one observable."""
    control, size, obs, err = data.xy(observable)
    out = {}
    for L in np.unique(size):
        m = size == L
        c, o, e = control[m], obs[m], err[m]
        order = np.argsort(c)
        c, o, e = c[order], o[order], e[order]
        if lo is not None:
            keep = c >= lo
            c, o, e = c[keep], o[keep], e[keep]
        if hi is not None:
            keep = c <= hi
            c, o, e = c[keep], o[keep], e[keep]
        if len(c) < 2:
            continue
        out[float(L)] = (c, o, e)
    return out


def crossings(data: FSSData, observable: str, size_pairs: Optional[Sequence[tuple]] = None,
              control_window: Optional[tuple] = None) -> list:
    """Compute pairwise crossing control values for a dimensionless observable.

    ``size_pairs`` defaults to consecutive sizes in ascending order.
    Each returned row is a dict: ``L1``, ``L2``, ``t_x``, ``t_x_err``,
    ``flat`` (True when the crossing is nearly parallel, i.e. useless).

    Notes
    -----
    The crossing of two *near-parallel* curves carries a huge error and is
    not an estimator of t_c; rows are flagged ``flat`` instead of silently
    trusted.
    """
    lo = hi = None
    if control_window:
        lo, hi = control_window
    by = _by_size(data, observable, lo, hi)
    if control_window:
        ref = 0.5 * (lo + hi)
    else:
        cmin, cmax = data.control_range()
        ref = 0.5 * (cmin + cmax)
    sizes = sorted(by)
    if len(sizes) < 2:
        return []
    if size_pairs is None:
        size_pairs = list(zip(sizes, sizes[1:]))
    rows = []
    # process pairs in ascending (L1, L2) so the crossing of a smaller pair
    # can serve as the reference that disambiguates multiple roots of the
    # next pair (a quadratic-in-t observable crosses twice; crossings converge
    # to t_c as L -> inf, so the previous estimate picks the physical one)
    prev_t = None
    for L1, L2 in sorted(set(size_pairs)):
        if L1 not in by or L2 not in by:
            continue
        c1, y1, e1 = by[L1]
        c2, y2, e2 = by[L2]
        t, terr = _crossing_root(c1, y1, e1, c2, y2, e2, lo, hi,
                                 prev_t if prev_t is not None else ref)
        if not np.isfinite(t):
            continue
        if np.isfinite(terr):
            prev_t = float(t)
        # dense-grid slope near the root, for the flat/near-parallel flag
        cmin = max(float(c1.min()), float(c2.min()), lo if lo is not None else -np.inf)
        cmax = min(float(c1.max()), float(c2.max()), hi if hi is not None else np.inf)
        dgrid = np.linspace(cmin, cmax, 4001)
        d = _linear_interp(c1, y1, dgrid) - _linear_interp(c2, y2, dgrid)
        slope = _local_slope(dgrid, d, t)
        control_span = max(float(np.ptp(c1)), float(np.ptp(c2)))
        flat = (not np.isfinite(terr)) or abs(slope) < 1e-9 or terr > 2.0 * control_span
        rows.append({
            "L1": float(L1), "L2": float(L2),
            "t_x": float(t), "t_x_err": float(terr),
            "flat": bool(flat),
        })
    return rows


def crossing_fit(rows: list, with_errors: bool = True, tc_init=None, **fit_kwargs) -> _fit.FitResult:
    """Fit ``t_x(L) = t_c + a L^{-lambda}`` to crossing rows.

    Rows flagged ``flat`` or with non-finite errors are excluded.
    Returns a :class:`FitResult`.
    """
    rows = [r for r in rows if not r.get("flat", False)]
    if not rows:
        raise ValueError("no usable crossing points (all flat or empty)")
    size = np.array([r["L1"] for r in rows], float)
    tx = np.array([r["t_x"] for r in rows], float)
    err = np.array([r["t_x_err"] for r in rows], float) if with_errors else None
    if tc_init is None:
        tc_init = float(np.median(tx))
    spec = crossing_spec(tc_init=tc_init)
    return _fit.fit_spec(spec, np.zeros_like(size), size, tx, err, **fit_kwargs)


def linear_amplitude_check(data: FSSData, observable: str, size: float,
                           pc: Optional[float] = None,
                           control_window: Optional[tuple] = None,
                           yt: float = 1.0) -> _fit.FitResult:
    """Fit ``R(t) = Rc + a1 t L^yt + a2 t^2 L^(2yt)`` at a single size to
    check a1 != 0.

    The amplitude ``a1`` is the leading field coefficient; if it is
    compatible with zero, crossings of this observable do not estimate t_c
    (see `references/dimensionless-crossings.md`).  ``yt`` and ``pc`` are
    *pinned*: on one size the leading exponent is not separable from the
    amplitudes, and a free ``pc`` leaks into ``a1`` through the ``a2`` cross
    term (``a1`` is only determined up to ``2 a2 pc L``).  Pass the values
    known from the transition class or from a prior derivative/crossing
    estimate.
    """
    sel = data.select_sizes(size, size)
    if control_window:
        sel = sel.select_control(*control_window)
    control, sz, obs, err = sel.xy(observable)
    # degree 2 so a quadratic-in-t observable does not alias its a1 onto
    # the missing a2 term
    spec = dimensionless_near_critical_spec(
        degree=2, with_pc=True, yt_fixed=yt, pc_fixed=pc)
    return _fit.fit_spec(spec, control, sz, obs, err, control_window=control_window)
