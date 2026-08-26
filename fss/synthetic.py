"""Deterministic synthetic FSS data for validating the core pipeline.

All generators return :class:`FSSData` in the canonical long format and
accept a ``seed`` for a numpy ``default_rng`` so tests are reproducible.
With ``noise=0.0`` the data are exact; the fits must then recover the
input parameters to within numerical tolerance.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .io import FSSData

__all__ = [
    "pure_power_data",
    "correction_power_data",
    "dimensionless_data",
]


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def _frame(control, size, obs, err) -> FSSData:
    df = pd.DataFrame({
        "control": np.asarray(control, float),
        "L": np.asarray(size, float),
        "obs": np.asarray(obs, float),
        "obs_err": np.asarray(err, float),
    })
    return FSSData.from_frame(
        df, control="control", size="L",
        observables={"obs": "obs"}, errors={"obs": "obs_err"},
        dimension=3, geometry="hypercubic",
    )


def pure_power_data(sizes, y=2.5, a=1.0, err_fraction=0.01, noise=1.0, seed=0) -> FSSData:
    """O(L) = a L^y with multiplicative noise.

    The observable is taken at a fixed control value (t = 0), so every row
    shares ``control = 0`` and the data support a critical power-law fit.
    """
    sizes = np.asarray(sizes, float)
    pred = a * sizes ** y
    rng = _rng(seed)
    err = err_fraction * pred
    obs = pred * (1.0 + noise * err_fraction * rng.standard_normal(sizes.size))
    return _frame(np.zeros(sizes.size), sizes, obs, err)


def correction_power_data(sizes, y=2.5, a=1.0, b=-0.5, omega=1.0,
                          err_fraction=0.01, noise=1.0, seed=0) -> FSSData:
    """O(L) = a L^y + b L^{y - omega}, the amplitude-correction form."""
    sizes = np.asarray(sizes, float)
    pred = a * sizes ** y + b * sizes ** (y - omega)
    rng = _rng(seed)
    err = err_fraction * np.abs(pred)
    obs = pred + noise * err * rng.standard_normal(sizes.size)
    return _frame(np.zeros(sizes.size), sizes, obs, err)


def dimensionless_data(
    sizes,
    controls,
    Rc=0.592,
    yt=1.0,
    a1=1.0,
    a2=0.5,
    pc=0.0,
    b1=0.0,
    yi=-1.0,
    err_fraction=0.01,
    noise=1.0,
    seed=0,
) -> FSSData:
    """R(t,L) = Rc + a1 t L^yt + a2 t^2 L^{2 yt} + b1 L^yi, t = control - pc.

    The workhorse of the near-critical dimensionless-observable tests.
    """
    sizes = np.asarray(sizes, float)
    controls = np.asarray(controls, float)
    L, T = np.meshgrid(sizes, controls, indexing="ij")
    t = T - pc
    pred = (
        Rc
        + a1 * t * L ** yt
        + a2 * t ** 2 * L ** (2 * yt)
        + b1 * L ** yi
    )
    rng = _rng(seed)
    err = err_fraction * np.abs(pred) + 1e-6
    obs = pred + noise * err * rng.standard_normal(pred.shape)
    # control is the thermal field T, size is the linear size L
    return _frame(T.ravel(), L.ravel(), obs.ravel(), err.ravel())
