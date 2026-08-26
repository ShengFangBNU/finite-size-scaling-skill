"""Parameterized FSS ansatz library.

Each model is a :class:`ModelSpec`: a name, an ordered parameter list
(with init values, bounds, and pinned/fixed values), and a ``predict``
callable.  Exponents and amplitudes may be free or pinned; pinned values
are removed from the optimization but reported in the result.

Model conventions (matching the reference papers):

- ``t`` denotes the (possibly uncentered) control variable; when a model
  has a ``pc`` parameter, the thermal field is ``t = control - pc``.
- Correction exponents ``y_i`` are negative; ``-omega`` is used for the
  leading irrelevant exponent of a power-law amplitude.
- Amplitudes are always in the *product* form ``b * L^yi`` (i.e. an
  ``a L^y (1 + b L^-omega)`` amplitude correction is equivalent to
  ``a L^y + a b L^{y-omega}``; the second amplitude absorbs the ``a``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "ParamSpec",
    "ModelSpec",
    "critical_power_spec",
    "critical_power_correction_spec",
    "scaling_observable_spec",
    "dimensionless_near_critical_spec",
    "crossing_spec",
]


@dataclass(frozen=True)
class ParamSpec:
    name: str
    init: float
    bounds: Tuple[Optional[float], Optional[float]] = (None, None)
    fixed: Optional[float] = None


@dataclass
class ModelSpec:
    name: str
    param_names: list
    predict: Callable[[dict, np.ndarray, np.ndarray], np.ndarray]
    params: list = field(default_factory=list)  # list[ParamSpec]
    describe: str = ""

    def __post_init__(self):
        by_name = {p.name: p for p in self.params}
        missing = [n for n in self.param_names if n not in by_name]
        if missing:
            raise ValueError(f"ModelSpec {self.name}: no ParamSpec for {missing}")

    def free_names(self) -> list:
        return [p.name for p in self.params if p.fixed is None]

    def pinned(self) -> dict:
        return {p.name: p.fixed for p in self.params if p.fixed is not None}

    def init_values(self) -> dict:
        return {p.name: p.init for p in self.params}

    def bounds_array(self, names: Sequence[str]):
        lo = []
        hi = []
        by_name = {p.name: p for p in self.params}
        for n in names:
            p = by_name[n]
            lo.append(-np.inf if p.bounds[0] is None else p.bounds[0])
            hi.append(np.inf if p.bounds[1] is None else p.bounds[1])
        return np.array(lo), np.array(hi)

    def evaluate(self, params: dict, control: np.ndarray, size: np.ndarray) -> np.ndarray:
        """Predict with a partial parameter dict; pinned/absent params filled
        from the spec so submodels may be evaluated on free parameters only."""
        full = dict(self.pinned())
        full.update(params)
        return np.asarray(self.predict(full, np.asarray(control, float), np.asarray(size, float)))


# ----------------------------------------------------------------------
# power-law models
# ----------------------------------------------------------------------

def _critical_power(params, control, size):
    # O(L) = a L^y   (control unused: the observable is taken at a fixed control)
    return params["a"] * size ** params["y"]


def critical_power_spec(y_init: float = 1.5, a_init: float = 1.0) -> ModelSpec:
    p = [
        ParamSpec("a", a_init, (None, None)),
        ParamSpec("y", y_init, (0.0, 100.0)),
    ]
    return ModelSpec(
        name="critical_power",
        param_names=["a", "y"],
        predict=_critical_power,
        params=p,
        describe="O(L) = a L^y",
    )


def _critical_power_correction(params, control, size):
    # O(L) = a L^y + b L^{y - omega}   (amplitude-correction form)
    y = params["y"]
    return params["a"] * size ** y + params["b"] * size ** (y - params["omega"])


def critical_power_correction_spec(
    omega: float = 1.0,
    omega_fixed: bool = True,
    y_init: float = 1.5,
) -> ModelSpec:
    """O(L) = a L^y + b L^{y - omega}.

    Equivalent to the amplitude form ``a L^y (1 + (b/a) L^{-omega})``.
    ``omega`` is the leading irrelevant exponent; by default it is pinned
    (the reference papers pin it after a free consistency check).
    """
    p = [
        ParamSpec("a", 1.0),
        ParamSpec("y", y_init, (0.0, 100.0)),
        ParamSpec("b", 0.1),
        ParamSpec("omega", omega, (0.0, 10.0), fixed=omega if omega_fixed else None),
    ]
    return ModelSpec(
        name="critical_power_correction",
        param_names=["a", "y", "b", "omega"],
        predict=_critical_power_correction,
        params=p,
        describe=f"O(L) = a L^y + b L^(y-omega), omega={'fixed=%g' % omega if omega_fixed else 'free'}",
    )


# ----------------------------------------------------------------------
# general scaling observable with corrections and background
# ----------------------------------------------------------------------
#   O(L) = c0 + L^y ( a0 + sum_i b_i L^y_i )
#


def _make_scaling_observable_predict(correction_terms):
    # correction_terms: list of dicts with keys 'name' (amplitude param),
    # 'exponent_name' and 'exponent' (float when pinned)
    def predict(params, control, size):
        y = params["y"]
        out = params["a0"] * size ** y
        for term in correction_terms:
            yi = params.get(term["exponent_name"], term["exponent"])
            out = out + params[term["name"]] * size ** (y + yi)
        if "c0" in params:
            out = out + params["c0"]
        return out

    return predict


def scaling_observable_spec(
    correction_exponents: Sequence[Tuple[float, str]] = (),
    background: bool = False,
    y_init: float = 1.5,
    a0_init: float = 1.0,
    c0_init: float = 0.0,
) -> ModelSpec:
    """O(L) = c0 + L^y (a0 + sum_i b_i L^y_i).

    ``correction_exponents`` is a sequence of ``(y_i, mode)`` with mode
    ``'fixed'`` or ``'free'``.  With ``background=True`` an analytic
    constant ``c0`` is added (see `analytic-backgrounds.md`).
    """
    params = [
        ParamSpec("y", y_init, (0.0, 100.0)),
        ParamSpec("a0", a0_init),
    ]
    names = ["y", "a0"]
    terms = []
    for i, (yi, mode) in enumerate(correction_exponents):
        if mode not in ("fixed", "free"):
            raise ValueError(f"mode must be 'fixed' or 'free', got {mode!r}")
        bname = f"b{i}"
        ename = f"yi{i}" if mode == "free" else None
        params.append(ParamSpec(bname, 0.1))
        names.append(bname)
        if mode == "free":
            params.append(ParamSpec(ename, yi, (-20.0, 0.0)))
            names.append(ename)
        terms.append({"name": bname, "exponent_name": ename, "exponent": yi if mode == "fixed" else None})
    if background:
        params.append(ParamSpec("c0", c0_init))
        names.append("c0")
    desc = "O(L) = " + ("c0 + " if background else "") + f"L^y (a0 + sum_i b_i L^y_i)"
    return ModelSpec(
        name="scaling_observable",
        param_names=names,
        predict=_make_scaling_observable_predict(terms),
        params=params,
        describe=desc,
    )


# ----------------------------------------------------------------------
# dimensionless near-critical observable
# ----------------------------------------------------------------------
#   R(t,L) = Rc + sum_{k=1..K} a_k t^k L^{k y_t}
#              + sum_i b_i L^y_i
#              + sum_j c_j t L^{y_t + y_j}          (mixed thermal-irrelevant)
#   t = control - pc   (when with_pc is True)
#


def _make_dimensionless_predict(degree, correction_terms, mixed_terms, with_pc):
    # correction_terms: list of dicts {'name','exponent_name','exponent'}
    # mixed_terms:      list of dicts {'name','exponent_name','exponent'}
    def predict(params, control, size):
        t = control
        if with_pc:
            t = control - params["pc"]
        yt = params["yt"]
        out = np.full_like(size, params["Rc"], dtype=float)
        for k in range(1, degree + 1):
            out = out + params[f"a{k}"] * (t ** k) * size ** (k * yt)
        for term in correction_terms:
            yi = params.get(term["exponent_name"], term["exponent"])
            out = out + params[term["name"]] * size ** yi
        for term in mixed_terms:
            yi = params.get(term["exponent_name"], term["exponent"])
            out = out + params[term["name"]] * t * size ** (yt + yi)
        return out

    return predict


def dimensionless_near_critical_spec(
    degree: int = 2,
    correction_exponents: Sequence[Tuple[float, str]] = (),
    mixed_exponents: Sequence[Tuple[float, str]] = (),
    with_pc: bool = True,
    yt_init: float = 1.5,
    yt_fixed: Optional[float] = None,
    Rc_init: float = 0.5,
    pc_init: float = 0.0,
    pc_fixed: Optional[float] = None,
) -> ModelSpec:
    """Dimensionless near-critical observable.

    ``degree`` is the polynomial order in the thermal field.  Corrections
    and mixed terms use the same ``(y_i, mode)`` convention as
    :func:`scaling_observable_spec`.  With ``with_pc=True`` the control is
    uncentered and ``pc`` is a fitted parameter; otherwise ``control`` is
    the centered field ``t`` and the model uses it directly.

    ``yt_fixed`` and ``pc_fixed`` pin the leading thermal exponent and the
    critical-point location.  On a single size the exponent is not separable
    from the amplitudes ``a_k`` (any ``yt`` can be absorbed into them), and a
    free ``pc`` leaks into ``a1`` through the ``a2`` cross term
    (``a1`` is then only determined up to ``2 a2 pc L``), so one-size checks
    like :func:`fss.crossing.linear_amplitude_check` must hold both fixed.
    """
    if degree < 1:
        raise ValueError("degree must be >= 1")
    params = []
    names = []
    if with_pc:
        params.append(ParamSpec(
            "pc", pc_fixed if pc_fixed is not None else pc_init,
            (None, None), fixed=pc_fixed))
        names.append("pc")
    params.append(ParamSpec("Rc", Rc_init))
    names.append("Rc")
    params.append(ParamSpec(
        "yt", yt_fixed if yt_fixed is not None else yt_init,
        (0.0, 100.0), fixed=yt_fixed))
    names.append("yt")
    for k in range(1, degree + 1):
        params.append(ParamSpec(f"a{k}", 0.1))
        names.append(f"a{k}")
    correction_terms = []
    for i, (yi, mode) in enumerate(correction_exponents):
        if mode not in ("fixed", "free"):
            raise ValueError(f"mode must be 'fixed' or 'free', got {mode!r}")
        bname = f"cb{i}"
        ename = f"cy{i}" if mode == "free" else None
        params.append(ParamSpec(bname, 0.0))
        names.append(bname)
        if mode == "free":
            params.append(ParamSpec(ename, yi, (-20.0, 0.0)))
            names.append(ename)
        correction_terms.append(
            {"name": bname, "exponent_name": ename, "exponent": yi if mode == "fixed" else None}
        )
    mixed_terms = []
    for j, (yj, mode) in enumerate(mixed_exponents):
        if mode not in ("fixed", "free"):
            raise ValueError(f"mode must be 'fixed' or 'free', got {mode!r}")
        cname = f"mx{j}"
        ename = f"my{j}" if mode == "free" else None
        params.append(ParamSpec(cname, 0.0))
        names.append(cname)
        if mode == "free":
            params.append(ParamSpec(ename, yj, (-20.0, 0.0)))
            names.append(ename)
        mixed_terms.append(
            {"name": cname, "exponent_name": ename, "exponent": yj if mode == "fixed" else None}
        )
    desc = (
        f"R(t,L) = Rc + sum_k a_k t^k L^(k yt)  [K={degree}]"
        + (", t=control-pc" if with_pc else ", t=control")
        + f", corrections={len(correction_terms)}, mixed={len(mixed_terms)}"
    )
    return ModelSpec(
        name="dimensionless_near_critical",
        param_names=names,
        predict=_make_dimensionless_predict(degree, correction_terms, mixed_terms, with_pc),
        params=params,
        describe=desc,
    )


# ----------------------------------------------------------------------
# pairwise crossing trajectory
# ----------------------------------------------------------------------
#   t_x(L) = t_c + a L^{-lambda}
#
#   t_x(L, sL) is the control value where two R(t, L) curves cross.  For a
#   genuine critical point it converges to t_c as L -> inf with exponent
#   lambda ~ y_t + |y_i| (see `references/dimensionless-crossings.md`).


def _crossing_predict(params, control, size):
    # control is unused; the trajectory lives in the size axis only
    return params["t_c"] + params["a"] * size ** (-params["lam"])


def crossing_spec(tc_init: float = 0.0, lam_init: float = 1.0) -> ModelSpec:
    p = [
        ParamSpec("t_c", tc_init),
        ParamSpec("a", 1.0),
        ParamSpec("lam", lam_init, (0.0, 10.0)),
    ]
    return ModelSpec(
        name="crossing",
        param_names=["t_c", "a", "lam"],
        predict=_crossing_predict,
        params=p,
        describe="t_x(L) = t_c + a L^(-lambda)",
    )
