"""End-to-end: synthetic data through the whole ordinary-FSS pipeline,
plus smoke tests that the CLI scripts run (validating the sys.path
bootstrap and argparse wiring).
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import fss

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_full_pipeline_dimensionless():
    seed = 123
    data = fss.synthetic.dimensionless_data(
        [16, 24, 32, 48, 64, 96, 128], [-0.25, -0.125, 0.0, 0.125, 0.25],
        Rc=0.592, yt=1.0, a1=1.2, a2=0.4, pc=0.015, b1=0.25, yi=-1.0,
        err_fraction=0.01, noise=1.0, seed=seed)
    control, size, obs, err = data.xy("obs")

    # 1) dimensionless near-critical fit
    fit = fss.fit_dimensionless(control, size, obs, err, degree=2,
                                correction_exponents=[(-1.0, "fixed")], pc_init=0.0)
    assert fit.value("pc") == pytest.approx(0.015, abs=0.01)
    assert fit.value("yt") == pytest.approx(1.0, abs=0.1)
    assert fit.value("Rc") == pytest.approx(0.592, abs=0.02)

    # 2) crossing analysis
    rows = fss.crossings(data, "obs")
    usable = [r for r in rows if not r["flat"]]
    xfit = fss.crossing.crossing_fit(usable)
    assert xfit.value("t_c") == pytest.approx(0.015, abs=0.05)

    # 3) derivative scaling
    gsize, g, gerr = fss.derivative.control_derivative(data, "obs", at=0.015)
    gfit = fss.derivative.derivative_scaling_fit(gsize, g, gerr)
    assert gfit.value("y") == pytest.approx(1.0, abs=0.15)

    # 4) stability: L_min scan of the dimensionless fit must be flat in pc
    def fn(lmin=None, lmax=None, **kw):
        return fss.fit_dimensionless(control, size, obs, err, degree=2,
                                     correction_exponents=[(-1.0, "fixed")],
                                     lmin=lmin, lmax=lmax)

    lrows = fss.diagnostics.lmin_scan(fn, [24, 32, 48, 64], lmax=128,
                                      param_names=["pc", "yt"])
    pcs = np.array([r["pc"] for r in lrows])
    assert np.ptp(pcs) < 0.05


def test_full_pipeline_critical_power():
    data = fss.synthetic.correction_power_data(
        [16, 24, 32, 48, 64, 96, 128], y=2.5, a=0.6, b=-0.3, omega=1.0,
        err_fraction=0.01, noise=1.0, seed=7)
    _, size, obs, err = data.xy("obs")

    res = fss.fit_critical_power(size, obs, err, correction=True, omega=1.0)
    assert res.value("y") == pytest.approx(2.5, abs=0.1)

    # collapse-quality metric on the corrected observable vs size-collapse x
    q = fss.collapse.collapse_quality(size, obs / size ** res.value("y"), err)
    assert q["chi2_reduced"] < 3.0


def _write_synthetic_csv(tmp_path):
    data = fss.synthetic.pure_power_data([16, 32, 64, 128], y=2.0, a=1.0,
                                         err_fraction=0.01, seed=1)
    p = tmp_path / "data.csv"
    data.df.to_csv(p, index=False)
    return p


def _run_script(name, *args, tmp_path):
    script = REPO_ROOT / "scripts" / name
    cmd = [sys.executable, str(script), *map(str, args)]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT)


def test_script_critical_power_fit_runs(tmp_path):
    p = _write_synthetic_csv(tmp_path)
    r = _run_script("critical_power_fit.py", p, "--obs", "obs", tmp_path=tmp_path)
    assert r.returncode == 0, r.stderr
    assert "critical_power" in r.stdout
    assert "y" in r.stdout


def test_script_crossing_runs(tmp_path):
    data = fss.synthetic.dimensionless_data(
        [16, 32, 64], [-0.2, -0.1, 0.0, 0.1, 0.2], Rc=0.5, yt=1.0,
        a1=1.0, a2=0.2, pc=0.0, b1=0.1, yi=-1.0, err_fraction=1e-4,
        noise=0.0, seed=2)
    p = tmp_path / "dim.csv"
    data.df.to_csv(p, index=False)
    r = _run_script("crossing.py", p, "--obs", "obs", tmp_path=tmp_path)
    assert r.returncode == 0, r.stderr
    assert "pairwise crossings" in r.stdout


def test_script_stability_scan_runs(tmp_path):
    p = _write_synthetic_csv(tmp_path)
    r = _run_script("stability_scan.py", p, "--obs", "obs",
                    "--lmin-values", "16,32", tmp_path=tmp_path)
    assert r.returncode == 0, r.stderr
    assert "L_min scan" in r.stdout


def test_script_inspect_data_runs(tmp_path):
    p = _write_synthetic_csv(tmp_path)
    r = _run_script("inspect_data.py", p, tmp_path=tmp_path)
    assert r.returncode == 0, r.stderr
    assert "rows:" in r.stdout
