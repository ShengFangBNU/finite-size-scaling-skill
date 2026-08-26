"""ModelSpec ansatz library: evaluation, pinning, bounds."""

import numpy as np
import pytest

from fss.models import (
    critical_power_spec,
    critical_power_correction_spec,
    scaling_observable_spec,
    dimensionless_near_critical_spec,
    crossing_spec,
    ModelSpec,
)


def test_critical_power_predict():
    spec = critical_power_spec()
    pred = spec.evaluate({"a": 2.0, "y": 1.5}, np.zeros(3), np.array([4.0, 16.0, 64.0]))
    np.testing.assert_allclose(pred, 2.0 * np.array([4.0, 16.0, 64.0]) ** 1.5)
    assert spec.free_names() == ["a", "y"]
    assert spec.pinned() == {}


def test_critical_power_correction_pinning():
    spec = critical_power_correction_spec(omega=1.0, omega_fixed=True)
    assert spec.pinned() == {"omega": 1.0}
    assert spec.free_names() == ["a", "y", "b"]
    free = critical_power_correction_spec(omega_fixed=False)
    assert "omega" in free.free_names()


def test_scaling_observable_with_background():
    spec = scaling_observable_spec(correction_exponents=[(-1.0, "fixed")], background=True)
    size = np.array([8.0, 16.0])
    pred = spec.evaluate({"y": 2.0, "a0": 1.0, "b0": 0.5, "c0": -2.0}, np.zeros(2), size)
    expected = -2.0 + size ** 2.0 * (1.0 + 0.5 * size ** (-1.0))
    np.testing.assert_allclose(pred, expected)
    assert "c0" in spec.param_names


def test_free_correction_exponent_is_optimized():
    spec = scaling_observable_spec(correction_exponents=[(-1.0, "free")])
    names = spec.free_names()
    assert "yi0" in names
    lo, hi = spec.bounds_array(names)
    assert hi[names.index("yi0")] == 0.0  # corrections have negative exponents


def test_dimensionless_near_critical():
    spec = dimensionless_near_critical_spec(
        degree=2, with_pc=True,
        correction_exponents=[(-1.0, "fixed")], mixed_exponents=[(-1.0, "fixed")],
    )
    control = np.array([0.1, -0.1])
    size = np.array([16.0, 16.0])
    params = {"pc": 0.0, "Rc": 0.5, "yt": 1.0, "a1": 1.0, "a2": 0.5,
              "cb0": 0.1, "mx0": 0.2}
    pred = spec.evaluate(params, control, size)
    t = control - 0.0
    expected = (0.5 + 1.0 * t * 16.0 + 0.5 * t ** 2 * 256.0
                + 0.1 * 16.0 ** (-1.0) + 0.2 * t * 16.0 ** (1.0 - 1.0))
    np.testing.assert_allclose(pred, expected)


def test_dimensionless_without_pc_uses_control_directly():
    spec = dimensionless_near_critical_spec(degree=1, with_pc=False)
    assert "pc" not in spec.param_names
    pred = spec.evaluate({"Rc": 0.5, "yt": 1.0, "a1": 2.0}, np.array([0.1]), np.array([8.0]))
    assert pred[0] == pytest.approx(0.5 + 2.0 * 0.1 * 8.0)


def test_crossing_spec():
    spec = crossing_spec(tc_init=0.5)
    pred = spec.evaluate({"t_c": 0.5, "a": 2.0, "lam": 1.0}, np.zeros(2), np.array([4.0, 16.0]))
    np.testing.assert_allclose(pred, [0.5 + 2.0 / 4.0, 0.5 + 2.0 / 16.0])


def test_missing_paramspec_raises():
    with pytest.raises(ValueError):
        ModelSpec(name="bad", param_names=["a", "b"],
                  predict=lambda p, c, s: s, params=[])
