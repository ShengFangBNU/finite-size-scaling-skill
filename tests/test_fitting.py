"""Fit engine: recovery of exact data, filtering, corrections, diagnostics."""

import numpy as np
import pytest

import fss
from fss import fitting as ft


def test_exact_pure_power_recovery():
    sizes = np.array([8.0, 16.0, 32.0, 64.0, 128.0])
    y, a = 2.5, 0.7
    obs = a * sizes ** y
    res = fss.fit_critical_power(sizes, obs)
    assert res.value("y") == pytest.approx(y, abs=1e-5)
    assert res.value("a") == pytest.approx(a, abs=1e-5)
    assert res.n_points == 5
    assert res.dof == 3
    assert res.chi2 < 1e-8


def test_power_fit_with_noise_has_reasonable_chi2():
    data = fss.synthetic.pure_power_data([16, 32, 64, 128], y=2.5, a=0.7,
                                         err_fraction=0.02, noise=1.0, seed=42)
    _, size, obs, err = data.xy("obs")
    res = fss.fit_critical_power(size, obs, err)
    assert res.value("y") == pytest.approx(2.5, abs=0.1)
    assert 0.3 < res.chi2_reduced < 3.0
    assert res.max_abs_corr < 0.99


def test_lmin_filters_points():
    sizes = np.array([8.0, 16.0, 32.0, 64.0, 128.0])
    obs = 1.0 * sizes ** 2.0
    res = fss.fit_critical_power(sizes, obs, lmin=32.0)
    assert res.n_points == 3
    assert res.lmin == 32.0
    assert res.size.min() >= 32.0


def test_control_window_filters():
    data = fss.synthetic.dimensionless_data(
        [16, 32, 64], [-0.3, -0.15, 0.0, 0.15, 0.3], pc=0.0,
        err_fraction=1e-6, noise=0.0, seed=0)
    control, size, obs, err = data.xy("obs")
    res = fss.fit_dimensionless(control, size, obs, err, degree=2,
                                control_window=(-0.15, 0.15))
    assert res.control.min() >= -0.15
    assert res.control.max() <= 0.15
    assert res.control_window == (-0.15, 0.15)


def test_dimensionless_fit_recovers_parameters():
    data = fss.synthetic.dimensionless_data(
        [16, 32, 64, 128], [-0.2, -0.1, 0.0, 0.1, 0.2],
        Rc=0.592, yt=1.0, a1=1.0, a2=0.5, pc=0.05, b1=0.2, yi=-1.0,
        err_fraction=1e-6, noise=0.0, seed=1)
    control, size, obs, err = data.xy("obs")
    res = fss.fit_dimensionless(control, size, obs, err, degree=2,
                                with_pc=True, correction_exponents=[(-1.0, "fixed")],
                                pc_init=0.0)
    assert res.value("pc") == pytest.approx(0.05, abs=1e-3)
    assert res.value("yt") == pytest.approx(1.0, abs=1e-3)
    assert res.value("Rc") == pytest.approx(0.592, abs=1e-3)


def test_correction_fit_recovers_y():
    data = fss.synthetic.correction_power_data(
        [8, 16, 32, 64, 128], y=2.0, a=1.0, b=-0.5, omega=1.0,
        err_fraction=1e-6, noise=0.0, seed=2)
    _, size, obs, err = data.xy("obs")
    # without correction, chi2 is large
    bare = fss.fit_critical_power(size, obs, err)
    assert bare.chi2_reduced > 5.0
    # with correction (omega pinned), y recovers
    res = fss.fit_critical_power(size, obs, err, correction=True, omega=1.0)
    assert res.value("y") == pytest.approx(2.0, abs=1e-3)


def test_fit_result_structure_complete():
    data = fss.synthetic.pure_power_data([16, 32, 64, 128, 256], y=1.5,
                                         err_fraction=0.01, seed=3)
    _, size, obs, err = data.xy("obs")
    res = fss.fit_critical_power(size, obs, err)
    for attr in ("params", "stderr", "cov", "corr", "chi2", "dof",
                 "chi2_reduced", "p_value", "aic", "aicc", "bic",
                 "residuals", "residuals_unweighted", "n_points", "n_params",
                 "n_total_params", "pinned", "warnings", "condition_number",
                 "max_abs_corr"):
        assert hasattr(res, attr)
    assert res.p_value > 0.0
    assert set(res.params) == {"a", "y"}
    assert set(res.pinned) == set()


def test_pinned_all_params_evaluates_chi2_only():
    from fss.models import ModelSpec, ParamSpec
    p = [
        ParamSpec("a", 1.0, fixed=1.0),
        ParamSpec("y", 1.5, fixed=2.0),
    ]
    spec = ModelSpec(name="pinned_power", param_names=["a", "y"],
                     predict=lambda params, c, s: params["a"] * s ** params["y"],
                     params=p)
    assert spec.free_names() == []
    assert spec.pinned() == {"a": 1.0, "y": 2.0}
    sizes = np.array([8.0, 16.0])
    obs = 1.0 * sizes ** 2.0
    res = ft.fit_spec(spec, np.zeros(2), sizes, obs)
    assert res.n_params == 0
    assert res.dof == 2
    assert res.chi2 == pytest.approx(0.0)


def test_illconditioned_warning_appears():
    # two nearly identical exponential bases produce a large correlation
    sizes = np.array([8.0, 16.0, 32.0, 64.0])
    obs = 1.0 * sizes ** 2.0
    from fss.models import scaling_observable_spec
    # try fitting two almost-degenerate correction exponents
    spec = scaling_observable_spec(
        correction_exponents=[(-1.0, "fixed"), (-1.001, "fixed")])
    res = ft.fit_spec(spec, np.zeros(4), sizes, obs, err=np.full(4, 1e-3))
    assert res.max_abs_corr > 0.9 or any("correlation" in w for w in res.warnings)
