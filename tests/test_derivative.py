"""Derivative estimators of y_t."""

import numpy as np
import pytest

import fss


def test_control_derivative_recovers_yt():
    data = fss.synthetic.dimensionless_data(
        [16, 32, 64, 128], [-0.2, -0.1, 0.0, 0.1, 0.2],
        Rc=0.5, yt=1.0, a1=1.0, a2=0.5, pc=0.0, b1=0.1, yi=-1.0,
        err_fraction=1e-6, noise=0.0, seed=0)
    size, g, gerr = fss.derivative.control_derivative(data, "obs", at=0.0)
    assert len(size) == 4
    # g = a1 L^yt exactly at t=pc
    np.testing.assert_allclose(g, 1.0 * size ** 1.0, rtol=1e-3)

    res = fss.derivative.derivative_scaling_fit(size, g, gerr)
    assert res.value("y") == pytest.approx(1.0, abs=1e-3)


def test_derivative_fit_with_noise():
    data = fss.synthetic.dimensionless_data(
        [16, 32, 64, 128, 256], [-0.2, -0.1, 0.0, 0.1, 0.2],
        Rc=0.5, yt=1.25, a1=1.0, a2=0.5, pc=0.0, b1=0.1, yi=-1.0,
        err_fraction=0.01, noise=1.0, seed=3)
    size, g, gerr = fss.derivative.control_derivative(data, "obs", at=0.0)
    res = fss.derivative.derivative_scaling_fit(size, g, gerr)
    assert res.value("y") == pytest.approx(1.25, abs=0.15)


def test_covariance_estimator():
    R = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    Nb = 2.0 * R
    g, ge = fss.derivative.covariance_estimator(R, Nb)
    assert g > 0
    assert np.isfinite(ge)
    assert ge > 0
    # uncorrelated -> g ~ 0
    rng = np.random.default_rng(0)
    g0, _ = fss.derivative.covariance_estimator(rng.normal(size=2000), rng.normal(size=2000))
    assert abs(g0) < 0.1


def test_control_derivative_edge_at_bounds():
    data = fss.synthetic.dimensionless_data(
        [16, 32], [0.0, 0.1, 0.2], Rc=0.5, yt=1.0, a1=1.0, a2=0.0,
        pc=0.0, err_fraction=1e-6, noise=0.0, seed=5)
    # at the left edge the derivative uses the one-sided difference
    size, g, _ = fss.derivative.control_derivative(data, "obs", at=0.0)
    assert len(size) == 2
    np.testing.assert_allclose(g, size, rtol=1e-3)
