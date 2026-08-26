"""Effective-exponent series, residual diagnostics, identifiability report."""

import numpy as np
import pytest

import fss


def test_effective_exponent_series_pure_power_is_flat():
    sizes = np.array([16.0, 32.0, 64.0, 128.0])
    obs = 2.0 * sizes ** 2.5
    series = fss.diagnostics.effective_exponent_series(sizes, obs)
    assert series.shape == (3, 3)
    np.testing.assert_allclose(series[:, 1], 2.5, rtol=1e-9)


def test_effective_exponent_series_detects_correction_drift():
    data = fss.synthetic.correction_power_data(
        [8, 16, 32, 64, 128, 256], y=2.0, a=1.0, b=-2.0, omega=1.0,
        err_fraction=1e-12, noise=0.0, seed=0)
    _, size, obs, err = data.xy("obs")
    series = fss.diagnostics.effective_exponent_series(size, obs)
    drift = series[-1, 1] - series[0, 1]
    # negative correction -> effective exponent approaches y from above
    assert series[0, 1] > 2.05
    assert series[-1, 1] < series[0, 1]


def test_effective_exponent_ratios():
    sizes = np.array([8.0, 16.0, 32.0, 64.0])
    obs = sizes ** 3.0
    series = fss.diagnostics.effective_exponent_ratios(sizes, obs, ratio=2)
    np.testing.assert_allclose(series[:, 1], 3.0)


def test_residual_diagnostics_random_residuals():
    rng = np.random.default_rng(4)
    size = np.array([16.0, 32.0, 64.0, 128.0])
    obs = 1.0 * size ** 2.0
    err = np.full(4, 0.02 * obs)
    res = fss.fit_critical_power(size, obs + rng.normal(scale=err), err)
    diag = fss.diagnostics.residual_diagnostics(res)
    assert set(diag) == {"runs_z", "runs_p", "corr_with_logsize", "corr_with_control"}
    assert np.isfinite(diag["runs_p"])


def test_identifiability_report_fields():
    # a constant plus an L^2 trend is well separated -> modest correlation
    sizes = 2 ** np.arange(4, 9)  # 16..256
    rng = np.random.default_rng(9)
    pred = 0.5 + 0.1 * sizes ** 2.0
    err = 0.01 * pred
    obs = pred + rng.normal(scale=err)
    from fss.models import ModelSpec, ParamSpec
    spec = ModelSpec(
        name="const_plus_L2",
        param_names=["a0", "c0"],
        predict=lambda params, c, s: params["c0"] + params["a0"] * s ** 2.0,
        params=[ParamSpec("a0", 0.1), ParamSpec("c0", 0.5)],
    )
    res = fss.fit_spec(spec, np.zeros_like(sizes), sizes, obs, err)
    rep = fss.diagnostics.identifiability_report(res)
    assert rep["max_abs_corr"] < 0.9
    assert rep["high_corr_pairs"] == []
    assert "condition_number" in rep
    assert rep["n_warnings"] == 0


def test_runs_test_detects_oscillation_in_residuals():
    # a structured W shape in residuals against size
    rng = np.random.default_rng(5)
    size = np.array([8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0])
    resid = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]) * rng.uniform(0.5, 1.5, 8)
    z, p = fss.statistics.runs_test(resid)
    assert z > 2.0
