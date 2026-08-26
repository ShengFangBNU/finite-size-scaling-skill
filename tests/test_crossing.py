"""Pairwise-crossing analysis and the a1 != 0 requirement."""

import numpy as np
import pytest

import fss


def _dim(a1, pc, seed=0, noise=0.0, err_fraction=1e-6, b1=0.1):
    # a fine control mesh (step 0.02): crossing accuracy requires that the
    # piecewise-linear interpolation of each R(t) curve resolve its curvature.
    # Sizes are in FIXED ratio 2: the crossing trajectory is a clean
    # t_x = t_c + a L^{-lambda} only for fixed-ratio pairs -- with varying
    # ratios the prefactor L1^yt (1 - s^yt) is itself size dependent and the
    # effective power law is washed out.
    controls = np.linspace(-0.2, 0.2, 21)
    return fss.synthetic.dimensionless_data(
        [16, 32, 64, 128, 256], controls,
        Rc=0.5, yt=1.0, a1=a1, a2=0.2, pc=pc, b1=b1, yi=-1.0,
        err_fraction=err_fraction, noise=noise, seed=seed)


def test_crossings_found_and_converge_to_pc():
    pc = 0.02
    data = _dim(a1=1.0, pc=pc, noise=0.5, err_fraction=1e-3, seed=7)
    rows = fss.crossings(data, "obs")
    assert len(rows) == 4  # 5 sizes -> 4 consecutive pairs
    usable = [r for r in rows if not r["flat"]]
    assert len(usable) == 4
    for r in usable:
        assert abs(r["t_x"] - pc) < 0.05

    res = fss.crossing.crossing_fit(usable)
    assert res.value("t_c") == pytest.approx(pc, abs=0.05)
    assert 0.5 < res.value("lam") < 3.0  # ~ yt + |yi| = 2


def test_crossing_fit_exact_data():
    data = _dim(a1=1.0, pc=0.02, noise=0.0, err_fraction=1e-6, seed=8)
    rows = fss.crossings(data, "obs")
    res = fss.crossing.crossing_fit([r for r in rows if not r["flat"]])
    assert res.value("t_c") == pytest.approx(0.02, abs=1e-3)


def test_a1_zero_gives_no_usable_crossings():
    # R(t,L) = Rc + a2 t^2 L^(2yt): curves touch only at t=0 -> flat
    data = _dim(a1=0.0, pc=0.0, b1=0.0, seed=9)
    rows = fss.crossings(data, "obs")
    assert rows, "expected at least a flat crossing at the touch point"
    assert all(r["flat"] for r in rows)
    with pytest.raises(ValueError):
        fss.crossing.crossing_fit(rows)


def test_linear_amplitude_check_detects_a1_zero():
    # positive amplitude (pc is pinned: known from the crossing estimate)
    data = _dim(a1=1.0, pc=0.02, noise=0.0, err_fraction=1e-5, seed=10)
    res = fss.crossing.linear_amplitude_check(data, "obs", size=64.0, pc=0.02)
    assert res.value("a1") > 0.2
    assert abs(res.value("a1") / max(res.error("a1"), 1e-300)) > 5.0

    # vanishing amplitude -> compatible with zero
    data0 = _dim(a1=0.0, pc=0.0, b1=0.0, noise=0.0, err_fraction=1e-4, seed=11)
    res0 = fss.crossing.linear_amplitude_check(data0, "obs", size=64.0, pc=0.0)
    assert abs(res0.value("a1")) < 3.0 * max(res0.error("a1"), 1e-6)
