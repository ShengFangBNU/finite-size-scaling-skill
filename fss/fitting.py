"""Weighted nonlinear least-squares fitting with a common result structure.

Every fit returns a :class:`FitResult` carrying (Phase H /
`references/fitting-and-systematics.md`):

- parameter estimates and standard errors, covariance and correlation matrices;
- chi2, dof, reduced chi2, p-value, AIC/AICc/BIC;
- weighted and unweighted residuals;
- bookkeeping (L_min, L_max, control window, ansatz description);
- identifiability flags (condition number, max |correlation|, warnings).

Statistical convention: residuals are weighted by the *reported* errors,
``r = (obs - pred) / sigma``, and the covariance is ``inv(J^T J)`` with
``J`` the weighted Jacobian at the solution.  A reduced chi2 far from one
is reported as a warning, not silently absorbed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
from scipy.optimize import least_squares

from . import statistics as _stats
from .models import ModelSpec, critical_power_spec, critical_power_correction_spec

__all__ = [
    "FitResult",
    "fit_spec",
    "fit_critical_power",
    "fit_dimensionless",
    "fit_scaling_observable",
]


@dataclass
class FitResult:
    model_name: str
    describe: str
    params: dict
    stderr: dict
    cov: np.ndarray
    corr: np.ndarray
    param_order: list
    chi2: float
    dof: int
    chi2_reduced: float
    p_value: float
    aic: float
    aicc: float
    bic: float
    residuals: np.ndarray              # weighted, (obs-pred)/sigma
    residuals_unweighted: np.ndarray   # obs - pred
    n_points: int
    n_params: int                      # free parameters
    n_total_params: int
    pinned: dict
    lmin: Optional[float] = None
    lmax: Optional[float] = None
    control_window: Optional[tuple] = None
    success: bool = True
    message: str = ""
    warnings: list = field(default_factory=list)
    condition_number: float = float("nan")
    max_abs_corr: float = float("nan")
    control: Optional[np.ndarray] = None
    size: Optional[np.ndarray] = None
    obs: Optional[np.ndarray] = None

    def value(self, name: str) -> float:
        return self.params[name]

    def error(self, name: str) -> float:
        return self.stderr.get(name, float("nan"))

    def summary(self) -> str:
        lines = [f"fit: {self.model_name}  |  {self.describe}"]
        lines.append(f"  chi2 = {self.chi2:.3f}  dof = {self.dof}  chi2/dof = {self.chi2_reduced:.3f}"
                     f"  p = {self.p_value:.4g}")
        lines.append(f"  AIC = {self.aic:.3f}  AICc = {self.aicc:.3f}  BIC = {self.bic:.3f}")
        lines.append(f"  points = {self.n_points}  free params = {self.n_params}"
                     f"  cond(cov) = {self.condition_number:.2g}"
                     f"  max|corr| = {self.max_abs_corr:.3f}")
        for name in self.param_order:
            if name in self.pinned:
                continue
            se = self.stderr.get(name)
            se_s = "       " if se is None or not np.isfinite(se) else f"{se:.5g}"
            lines.append(f"    {name:8s} = {self.params[name]: .6g}  +/- {se_s}")
        for name, val in self.pinned.items():
            lines.append(f"    {name:8s} = {val:.6g}   (pinned)")
        if self.warnings:
            lines.append("  warnings:")
            for w in self.warnings:
                lines.append(f"    - {w}")
        return "\n".join(lines)


def _mask_valid(control, size, obs, err):
    control = np.asarray(control, float)
    size = None if size is None else np.asarray(size, float)
    obs = np.asarray(obs, float)
    if err is None:
        err = np.ones_like(obs)
    else:
        err = np.asarray(err, float)
    keep = np.isfinite(obs) & np.isfinite(err) & (err > 0)
    if size is not None:
        keep &= np.isfinite(size)
    keep &= np.isfinite(control)
    return control[keep], (size[keep] if size is not None else None), obs[keep], err[keep]


def fit_spec(
    spec: ModelSpec,
    control,
    size,
    obs,
    err=None,
    *,
    lmin: Optional[float] = None,
    lmax: Optional[float] = None,
    control_window: Optional[tuple] = None,
    p0: Optional[dict] = None,
    scale_covariance: bool = False,
    max_nfev: int = 100_000,
    xtol: float = 1e-12,
) -> FitResult:
    """Fit ``spec`` to ``obs(control, size)`` with weights ``1/err**2``.

    Rows with non-finite values or non-positive errors are dropped.
    """
    control, size, obs, err = _mask_valid(control, size, obs, err)
    if size is None:
        raise ValueError("no size values available; FSS fits require a size axis")

    # apply the requested size window and control window (these subset the
    # data, so lmin_scan / window_scan drive genuine re-fits, not bookkeeping)
    if lmin is not None or lmax is not None:
        keep = np.isfinite(size)
        if lmin is not None:
            keep &= size >= lmin
        if lmax is not None:
            keep &= size <= lmax
        control, size, obs, err = control[keep], size[keep], obs[keep], err[keep]
    if control_window is not None:
        lo, hi = control_window
        keep = (control >= lo) & (control <= hi)
        control, size, obs, err = control[keep], size[keep], obs[keep], err[keep]

    free_names = spec.free_names()
    pinned = spec.pinned()
    n_free = len(free_names)
    n_points = len(obs)

    if n_free == 0:
        # nothing to fit: evaluate the model with all parameters pinned
        pred = spec.evaluate({}, control, size)
        resid = (obs - pred) / err
        chi2 = float(np.sum(resid**2))
        dof = n_points
        return _assemble(
            spec, free_names, {}, {}, np.zeros((0, 0)), chi2, dof, resid,
            control, size, obs, lmin=lmin, lmax=lmax, control_window=control_window,
        )

    init = spec.init_values()
    if p0 is not None:
        init.update(p0)
    x0 = np.array([init[n] for n in free_names], float)
    lo, hi = spec.bounds_array(free_names)
    lo = np.maximum(lo, x0 - 1e6)  # keep finite for trf; loose
    hi = np.minimum(hi, x0 + 1e6)

    def unpack(x):
        return dict(zip(free_names, x))

    def resid_fun(x):
        pred = spec.evaluate(unpack(x), control, size)
        return (obs - pred) / err

    result = least_squares(
        resid_fun, x0, bounds=(lo, hi), method="trf", jac="2-point",
        # 'jac' scales each parameter by the norm of its Jacobian column:
        # essential when parameters live on very different scales (e.g.
        # pc ~ 1e-2 alongside amplitudes ~ 1), otherwise trf can get
        # trapped in a poor local minimum.
        x_scale="jac",
        max_nfev=max_nfev, xtol=xtol, ftol=xtol, gtol=xtol,
    )
    x = result.x
    params = dict(pinned)
    params.update(unpack(x))

    resid = resid_fun(x)
    chi2 = float(np.sum(resid**2))
    dof = max(n_points - n_free, 0)

    # covariance from the weighted Jacobian
    J = result.jac
    JtJ = J.T @ J
    cond = float(np.linalg.cond(JtJ)) if JtJ.size else float("nan")
    if scale_covariance:
        scale = chi2 / dof if dof > 0 else 1.0
    else:
        scale = 1.0
    try:
        cov = np.linalg.inv(JtJ) * scale
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(JtJ) * scale
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    stderr = {name: float(s) for name, s in zip(free_names, se)}
    corr = _stats.correlation_from_cov(cov)
    max_corr = 0.0
    if corr.size > 1:
        off = np.abs(corr[~np.eye(corr.shape[0], dtype=bool)])
        if off.size:
            max_corr = float(off.max())

    return _assemble(
        spec, free_names, params, stderr, cov, chi2, dof, resid,
        control, size, obs, corr=corr, cond=cond, max_corr=max_corr,
        lmin=lmin, lmax=lmax, control_window=control_window,
        success=result.success, message=result.message, n_free=n_free,
    )


def _assemble(spec, free_names, params, stderr, cov, chi2, dof, resid,
              control, size, obs, corr=None, cond=float("nan"),
              max_corr=float("nan"), lmin=None, lmax=None,
              control_window=None, success=True, message="", n_free=None):
    n = len(obs)
    k = n_free if n_free is not None else len(free_names)
    reduced = chi2 / dof if dof > 0 else float("inf")
    p_val = _stats.chi2_p_value(chi2, dof)
    res = FitResult(
        model_name=spec.name,
        describe=spec.describe,
        params=params,
        stderr=stderr,
        cov=cov,
        corr=np.zeros((k, k)) if corr is None else corr,
        param_order=list(free_names),
        chi2=chi2,
        dof=dof,
        chi2_reduced=reduced,
        p_value=p_val,
        aic=_stats.aic(chi2, k),
        aicc=_stats.aicc(chi2, k, n),
        bic=_stats.bic(chi2, k, n),
        residuals=resid,
        residuals_unweighted=obs - spec.evaluate(params, control, size),
        n_points=n,
        n_params=k,
        n_total_params=len(spec.param_names),
        pinned=spec.pinned(),
        lmin=lmin,
        lmax=lmax,
        control_window=control_window,
        success=success,
        message=message,
        condition_number=cond,
        max_abs_corr=max_corr,
        control=np.array(control),
        size=np.array(size),
        obs=np.array(obs),
    )
    res.warnings = _identifiability_warnings(res)
    return res


def _identifiability_warnings(res: FitResult) -> list:
    warnings = []
    if not res.success:
        warnings.append(f"optimizer did not converge: {res.message}")
    if np.isfinite(res.condition_number) and res.condition_number > 1e8:
        warnings.append(
            f"ill-conditioned covariance (cond = {res.condition_number:.2g}); "
            "parameters are not separately identifiable"
        )
    if np.isfinite(res.max_abs_corr) and res.max_abs_corr > 0.99:
        warnings.append(
            f"extremely high parameter correlation (max|corr| = {res.max_abs_corr:.3f}); "
            "parameters may be degenerate - pin one and refit"
        )
    if np.isfinite(res.chi2_reduced) and res.dof > 0 and res.chi2_reduced > 3.0:
        warnings.append(
            f"reduced chi2 = {res.chi2_reduced:.2f} >> 1: missing corrections or "
            "mis-specified errors - increase L_min or extend the ansatz"
        )
    if np.isfinite(res.chi2_reduced) and res.dof > 0 and res.chi2_reduced < 0.3:
        warnings.append(
            f"reduced chi2 = {res.chi2_reduced:.2f} << 1: errors may be overestimated "
            "or the model overparameterized"
        )
    return warnings


# ----------------------------------------------------------------------
# convenience entry points
# ----------------------------------------------------------------------

def fit_critical_power(
    size,
    obs,
    err=None,
    *,
    correction: bool = False,
    omega: float = 1.0,
    omega_fixed: bool = True,
    y_init: float = 1.5,
    lmin: Optional[float] = None,
    lmax: Optional[float] = None,
    control_window: Optional[tuple] = None,
) -> FitResult:
    spec = (
        critical_power_correction_spec(omega=omega, omega_fixed=omega_fixed, y_init=y_init)
        if correction
        else critical_power_spec(y_init=y_init)
    )
    return fit_spec(
        spec, np.zeros_like(size), size, obs, err,
        lmin=lmin, lmax=lmax, control_window=control_window,
    )


def fit_dimensionless(
    control,
    size,
    obs,
    err=None,
    *,
    degree: int = 2,
    correction_exponents=(),
    mixed_exponents=(),
    with_pc: bool = True,
    yt_init: float = 1.5,
    pc_init: Optional[float] = None,
    lmin: Optional[float] = None,
    lmax: Optional[float] = None,
    control_window: Optional[tuple] = None,
) -> FitResult:
    from .models import dimensionless_near_critical_spec

    if pc_init is None and with_pc and len(control):
        pc_init = float(np.median(np.asarray(control, float)))
    spec = dimensionless_near_critical_spec(
        degree=degree,
        correction_exponents=correction_exponents,
        mixed_exponents=mixed_exponents,
        with_pc=with_pc,
        yt_init=yt_init,
        pc_init=pc_init or 0.0,
    )
    return fit_spec(
        spec, control, size, obs, err,
        lmin=lmin, lmax=lmax, control_window=control_window,
    )


def fit_scaling_observable(
    size,
    obs,
    err=None,
    *,
    correction_exponents=(),
    background: bool = False,
    y_init: float = 1.5,
    lmin: Optional[float] = None,
    lmax: Optional[float] = None,
    control_window: Optional[tuple] = None,
) -> FitResult:
    from .models import scaling_observable_spec

    spec = scaling_observable_spec(
        correction_exponents=correction_exponents,
        background=background,
        y_init=y_init,
    )
    return fit_spec(
        spec, np.zeros_like(size), size, obs, err,
        lmin=lmin, lmax=lmax, control_window=control_window,
    )
