"""Canonical data I/O: schema, guessing, selection, extraction."""

import numpy as np
import pandas as pd
import pytest

import fss
from fss import FSSData, FSSDataError, guess_columns


def _write_csv(path, text):
    path.write_text(text)
    return path


def test_read_csv_sniffs_commas(tmp_path):
    p = _write_csv(tmp_path / "d.csv", (
        "L,control,obs,obs_err\n"
        "8,0.0,1.0,0.01\n"
        "16,0.0,2.0,0.02\n"
    ))
    data = FSSData.load(p)
    assert data.size_col == "L"
    assert data.control_col == "control"
    assert data.observable_names() == ["obs"]
    assert data.err_cols == {"obs": "obs_err"}
    assert list(data.sizes()) == [8.0, 16.0]
    assert data.control_range() == (0.0, 0.0)


def test_read_tsv(tmp_path):
    p = _write_csv(tmp_path / "d.tsv", (
        "L\tcontrol\tR\tR_err\n"
        "8\t0.1\t1.0\t0.01\n"
    ))
    data = FSSData.load(p)
    assert data.observable_names() == ["R"]


def test_guess_columns_excludes_N_from_size(tmp_path):
    p = _write_csv(tmp_path / "d.csv", (
        "N,L,p,R,R_err\n"
        "26,8,0.5,1.0,0.01\n"
    ))
    df = pd.read_csv(p)
    g = guess_columns(df)
    assert g["size"] == "L", "N must never be guessed as the size"
    assert g["control"] == "p"
    assert g["observables"].get("R") == "R"
    assert g["errors"].get("R") == "R_err"


def test_volume_not_assumed_Ld():
    df = pd.DataFrame({"V": [512.0, 4096.0], "control": [0.0, 0.0],
                       "obs": [1.0, 2.0]})
    data = FSSData.from_frame(df, control="control", volume="V",
                              observables={"obs": "obs"})
    assert data.size_is_volume()
    assert data.volume_from_size() is None, "V = L^d must not be assumed without d"


def test_volume_from_size_with_dimension():
    df = pd.DataFrame({"L": [8.0, 16.0], "control": [0.0, 0.0], "obs": [1.0, 2.0]})
    data = FSSData.from_frame(df, control="control", size="L", dimension=3,
                              observables={"obs": "obs"})
    np.testing.assert_allclose(data.volume_from_size(), [512.0, 4096.0])


def test_selection_and_xy():
    df = pd.DataFrame({
        "L": [8, 8, 16, 16],
        "control": [0.0, 0.1, 0.0, 0.1],
        "obs": [1.0, 1.2, 2.0, 2.4],
        "obs_err": [0.01] * 4,
    })
    data = FSSData.from_frame(df, control="control", size="L",
                              observables={"obs": "obs"}, errors={"obs": "obs_err"})
    assert data.select_sizes(16, 16).sizes().tolist() == [16.0]
    sel = data.select_control_value(0.1)
    c, s, o, e = sel.xy("obs")
    np.testing.assert_allclose(np.sort(s), [8.0, 16.0])
    np.testing.assert_allclose(np.sort(c), [0.1, 0.1])
    # errors are dropped when non-positive / missing
    c, s, o, e = data.xy("obs")
    assert e.min() > 0


def test_missing_control_raises():
    df = pd.DataFrame({"L": [8.0], "obs": [1.0]})
    with pytest.raises(FSSDataError):
        FSSData.from_frame(df, control="nope")


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        FSSData.load(tmp_path / "nope.csv")


def test_info_dict(tmp_path):
    p = _write_csv(tmp_path / "d.csv", "L,control,obs,obs_err\n8,0.0,1.0,0.01\n16,0.0,2.0,0.02\n")
    info = FSSData.load(p).info()
    assert info["n_sizes"] == 2
    assert info["with_errors"] is True
