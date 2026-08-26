"""Statistical utilities: p-values, information criteria, correlations,
runs test, effective exponents."""

import numpy as np
import pytest

from fss import statistics as st


def test_chi2_p_value_known():
    # chi2=3.841 with 1 dof -> p ~ 0.05
    assert abs(st.chi2_p_value(3.841, 1) - 0.05) < 1e-3
    assert np.isnan(st.chi2_p_value(5.0, 0))


def test_information_criteria():
    assert st.aic(10.0, 2) == 14.0
    assert st.aic(10.0, 2) == 10.0 + 4.0
    # BIC = chi2 + k ln n
    assert abs(st.bic(10.0, 2, 50) - (10.0 + 2 * np.log(50))) < 1e-12
    assert np.isnan(st.aicc(10.0, 5, 5))  # n - k - 1 <= 0


def test_correlation_from_cov():
    cov = np.array([[4.0, 1.0], [1.0, 1.0]])
    corr = st.correlation_from_cov(cov)
    assert corr[0, 0] == pytest.approx(1.0)
    assert abs(corr[0, 1] - 0.5) < 1e-12
    assert abs(corr[1, 0] - 0.5) < 1e-12


def test_condition_number():
    assert st.condition_number(np.array([[2.0]])) == pytest.approx(1.0)
    assert np.isfinite(st.condition_number(np.array([[4.0, 0], [0, 1.0]])))


def test_runs_test_alternating():
    # perfectly alternating signs -> many runs -> large positive z
    r = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0])
    z, p = st.runs_test(r)
    assert z > 2.0
    assert p < 0.05


def test_runs_test_structured_trend():
    # all same sign -> extreme
    z, p = st.runs_test(np.array([1.0, 1.0, 1.0, 1.0, 1.0]))
    assert p == pytest.approx(0.0)


def test_effective_exponent_pair_exact():
    y, ye = st.effective_exponent_pair(1.0, 8.0, 4.0, 16.0)
    # ln(8)/ln(4) = 1.5
    assert y == pytest.approx(1.5)
    assert np.isnan(ye)  # no errors supplied


def test_effective_exponent_pair_errors():
    y, ye = st.effective_exponent_pair(1.0, 8.0, 4.0, 16.0, s1=0.05, s2=0.4)
    assert np.isfinite(ye) and ye > 0


def test_fit_report_dataclass():
    r = st.FitReport({}, {}, 1.0, 5, 0.2, 0.9, 1.0, 1.0, 1.0, 5, 1, [])
    assert r.chi2_reduced == pytest.approx(0.2)
