"""L_min and fitting-window stability scans."""

import numpy as np
import pytest

import fss


def test_lmin_scan_trajectory_and_rule():
    data = fss.synthetic.correction_power_data(
        [8, 16, 24, 32, 48, 64, 96, 128], y=2.0, a=1.0, b=-0.5, omega=1.0,
        err_fraction=0.02, noise=1.0, seed=0)
    _, size, obs, err = data.xy("obs")

    def fn(lmin=None, lmax=None, **kw):
        return fss.fit_critical_power(size, obs, err, correction=True, omega=1.0,
                                      lmin=lmin, lmax=lmax)

    rows = fss.diagnostics.lmin_scan(fn, [16, 24, 32, 48], lmax=128,
                                     param_names=["y"])
    assert len(rows) == 4
    # chi2/dof should settle as L_min grows (missing-correction bias shrinks)
    r0, r1 = rows[0], rows[-1]
    assert r1["chi2_reduced"] < r0["chi2_reduced"]
    # accepted y approaches the true 2.0
    assert abs(rows[-1]["y"] - 2.0) < 0.2
    # n_points shrinks with L_min
    assert rows[0]["n_points"] > rows[-1]["n_points"]
    # every row carries the required bookkeeping
    for k in ("lmin", "chi2", "dof", "chi2_reduced", "p_value", "n_points"):
        assert k in rows[0]


def test_lmin_rule_selects_good_cut():
    # a strong correction: the pure power law only becomes good at large L.
    # (b=-1.5 with tight errors makes chi2/dof clearly > 3 at the smallest
    # L_min; the weaker b=-0.8 / 1% combination was already acceptable at
    # lmin=16 and the rule rightly accepted it.)
    data = fss.synthetic.correction_power_data(
        [8, 16, 24, 32, 48, 64, 96, 128], y=2.0, a=1.0, b=-1.5, omega=1.0,
        err_fraction=0.005, noise=1.0, seed=2)
    _, size, obs, err = data.xy("obs")

    def fn(lmin=None, lmax=None, **kw):
        return fss.fit_critical_power(size, obs, err, correction=False,
                                      lmin=lmin, lmax=lmax)

    rows = fss.diagnostics.lmin_scan(fn, [16, 24, 32, 48, 64, 96], lmax=128)
    reduced = np.array([r["chi2_reduced"] for r in rows])
    # smallest L_min with chi2/dof ~ 1 (within, say, 2) is the accepted cut
    ok = np.where(reduced <= 3.0)[0]
    assert ok.size > 0
    accepted = rows[ok[0]]["lmin"]
    assert accepted >= 32.0


def test_window_scan_pc_stable():
    data = fss.synthetic.dimensionless_data(
        [16, 32, 64, 128], [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3],
        Rc=0.592, yt=1.0, a1=1.0, a2=0.5, pc=0.02, b1=0.2, yi=-1.0,
        err_fraction=0.005, noise=1.0, seed=2)
    control, size, obs, err = data.xy("obs")

    def fn(control_window=None, **kw):
        return fss.fit_dimensionless(control, size, obs, err, degree=2,
                                     correction_exponents=[(-1.0, "fixed")],
                                     control_window=control_window)

    rows = fss.diagnostics.window_scan(fn, 0.02, [0.1, 0.15, 0.2, 0.3],
                                       param_names=["pc", "yt"])
    assert len(rows) == 4
    pcs = np.array([r["pc"] for r in rows])
    # pc must be stable within its errors across the windows
    assert abs(pcs.mean() - 0.02) < 0.05
    assert np.ptp(pcs) < 0.05
    assert rows[0]["half_width"] == pytest.approx(0.1)
