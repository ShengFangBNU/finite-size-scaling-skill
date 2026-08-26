"""Canonical data I/O for the general FSS toolkit.

Schema (long format, one measurement per row):

    size       linear size L (the size variable actually controlled)
    volume     number of sites / system volume V (optional)
    control    thermal/field-like control parameter (t, beta, p, T, ...)
    <obs>      one column per observable
    <obs>_err  one column per observable uncertainty (optional)
    sample_id  replicate / run / seed id (optional, for raw-sample use)

Rules enforced here:

- The linear size and the volume are distinct fields.  Never silently
  identify V with L, and never assume V = L^d unless ``dimension`` and
  ``geometry`` are explicitly supplied.
- ``control`` is a single scalar field; scripts that need a control
  window or a critical value operate on this column.
- Columns are canonicalized once at load time; the rest of the toolkit
  consumes :class:`FSSData`, never raw dataframes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "FSSData",
    "load_fss_csv",
    "fss_from_frame",
    "guess_columns",
    "FSSDataError",
]

# Column-name guesses, in priority order.  ``N`` is deliberately excluded
# from size guesses: in many models N is a control parameter (e.g. the
# number of equivalent neighbors), not a size.
SIZE_GUESSES = ["L", "l", "size", "linear", "linear_size", "length"]
VOLUME_GUESSES = ["V", "vol", "volume", "sites", "Nsite"]
CONTROL_GUESSES = [
    "t", "T", "beta", "p", "r", "K", "J", "h", "control", "x",
    "lambda", "eps", "epsilon", "temperature", "param", "g",
]
ERROR_SUFFIXES = ("_err", "_error", "_std", "_sigma", "err", "error", "std", "sigma")
SAMPLE_ID_GUESSES = ["sample", "sample_id", "run", "replicate", "repl", "seed", "rep"]


class FSSDataError(ValueError):
    """Raised when data do not fit the canonical schema."""


def _is_number_like(s: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(s) or pd.api.types.is_bool_dtype(s)


def guess_columns(df: pd.DataFrame) -> dict:
    """Best-effort guess of canonical columns from a raw dataframe.

    Returns a dict with keys ``size``, ``volume``, ``control``,
    ``observables``, ``errors``, ``sample_id`` (the last three mapping
    canonical observable names to source column names, or empty).
    Guessing is only a convenience: every call to
    :func:`fss_from_frame` may override every entry.
    """
    result: dict = {
        "size": None,
        "volume": None,
        "control": None,
        "observables": {},
        "errors": {},
        "sample_id": None,
    }
    cols = list(df.columns)
    numeric_cols = [c for c in cols if _is_number_like(df[c])]
    nonnumeric = [c for c in cols if c not in numeric_cols]

    for guess in SIZE_GUESSES:
        if guess in numeric_cols:
            result["size"] = guess
            break
    for guess in VOLUME_GUESSES:
        if guess in numeric_cols:
            result["volume"] = guess
            break
    for guess in CONTROL_GUESSES:
        if guess in numeric_cols:
            result["control"] = guess
            break

    for c in cols:
        if c in (result["size"], result["volume"], result["control"]):
            continue
        low = c.lower()
        if any(low.endswith(s) or low == s for s in ERROR_SUFFIXES) and _is_number_like(df[c]):
            # error column; attach to the observable whose name it decorates
            base = c
            for suf in ("_err", "_error", "_std", "_sigma"):
                if base.lower().endswith(suf):
                    base = base[: -len(suf)]
                    break
            else:
                continue
            if base in numeric_cols:
                result["errors"][base] = c
        elif low in SAMPLE_ID_GUESSES or (c in nonnumeric and low.startswith("sample")):
            if result["sample_id"] is None:
                result["sample_id"] = c

    # everything numeric left over is a candidate observable
    for c in numeric_cols:
        if c in (result["size"], result["volume"], result["control"]):
            continue
        if c in result["errors"].values():
            continue  # an error column itself is not an observable
        result["observables"][c] = c

    return result


@dataclass
class FSSData:
    """Canonical container for finite-size data in long format."""

    df: pd.DataFrame
    control_col: str
    obs_cols: Mapping[str, str] = field(default_factory=dict)
    err_cols: Mapping[str, str] = field(default_factory=dict)
    size_col: Optional[str] = None
    volume_col: Optional[str] = None
    sample_id_col: Optional[str] = None
    dimension: Optional[float] = None
    geometry: Optional[str] = None

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    @classmethod
    def from_frame(
        cls,
        df: pd.DataFrame,
        *,
        control: str,
        size: Optional[str] = None,
        volume: Optional[str] = None,
        observables: Optional[Mapping[str, str]] = None,
        errors: Optional[Mapping[str, str]] = None,
        sample_id: Optional[str] = None,
        dimension: Optional[float] = None,
        geometry: Optional[str] = None,
    ) -> "FSSData":
        if control not in df.columns:
            raise FSSDataError(f"control column {control!r} not present")
        if size is not None and size not in df.columns:
            raise FSSDataError(f"size column {size!r} not present")
        if volume is not None and volume not in df.columns:
            raise FSSDataError(f"volume column {volume!r} not present")
        obs = dict(observables) if observables else {}
        for name, col in obs.items():
            if col not in df.columns:
                raise FSSDataError(f"observable column {col!r} (for {name!r}) not present")
        errs = dict(errors) if errors else {}
        for name, col in errs.items():
            if col not in df.columns:
                raise FSSDataError(f"error column {col!r} (for {name!r}) not present")
            if name not in obs:
                raise FSSDataError(f"error column {col!r} has no matching observable {name!r}")
        if dimension is not None and geometry is None:
            geometry = "unknown"
        return cls(
            df=df.copy(),
            control_col=control,
            size_col=size,
            volume_col=volume,
            obs_cols=obs,
            err_cols=errs,
            sample_id_col=sample_id,
            dimension=dimension,
            geometry=geometry,
        )

    @classmethod
    def load(
        cls,
        path,
        *,
        size: Optional[str] = None,
        volume: Optional[str] = None,
        control: Optional[str] = None,
        observables: Optional[Mapping[str, str]] = None,
        errors: Optional[Mapping[str, str]] = None,
        sample_id: Optional[str] = None,
        dimension: Optional[float] = None,
        geometry: Optional[str] = None,
    ) -> "FSSData":
        df = read_table(path)
        guess = guess_columns(df)
        return cls.from_frame(
            df,
            control=control if control is not None else guess["control"],
            size=size if size is not None else guess["size"],
            volume=volume if volume is not None else guess["volume"],
            observables=observables if observables is not None else guess["observables"],
            errors=errors if errors is not None else guess["errors"],
            sample_id=sample_id if sample_id is not None else guess["sample_id"],
            dimension=dimension,
            geometry=geometry,
        )

    # ------------------------------------------------------------------
    # accessors
    # ------------------------------------------------------------------
    def observable_names(self) -> list:
        return list(self.obs_cols)

    def sizes(self) -> np.ndarray:
        if self.size_col is None:
            return np.array([])
        return np.array(sorted(self.df[self.size_col].dropna().unique()))

    def control_range(self) -> tuple:
        vals = self.df[self.control_col].dropna()
        return (float(vals.min()), float(vals.max()))

    def has_size(self) -> bool:
        return self.size_col is not None

    def size_is_volume(self) -> bool:
        """True when the single size field is explicitly a volume, not a length."""
        return self.size_col is None and self.volume_col is not None

    def volume_from_size(self) -> Optional[np.ndarray]:
        """Return V = L^d only when d and geometry are explicitly known.

        Returns None (instead of guessing) whenever ``dimension`` is
        missing or the geometry is not stated as translationally regular.
        """
        if self.size_col is None or self.dimension is None:
            return None
        return self.df[self.size_col].to_numpy() ** self.dimension

    # ------------------------------------------------------------------
    # selection
    # ------------------------------------------------------------------
    def select_control(self, low: Optional[float] = None, high: Optional[float] = None) -> "FSSData":
        mask = np.ones(len(self.df), dtype=bool)
        if low is not None:
            mask &= self.df[self.control_col] >= low
        if high is not None:
            mask &= self.df[self.control_col] <= high
        return self._subset(mask)

    def select_sizes(self, lmin: Optional[float] = None, lmax: Optional[float] = None) -> "FSSData":
        if self.size_col is None:
            if lmin is None and lmax is None:
                return self
            raise FSSDataError("cannot select by size: no size column")
        mask = np.ones(len(self.df), dtype=bool)
        if lmin is not None:
            mask &= self.df[self.size_col] >= lmin
        if lmax is not None:
            mask &= self.df[self.size_col] <= lmax
        return self._subset(mask)

    def select_control_value(self, control_value: float, atol: float = 1e-9) -> "FSSData":
        """Rows whose control equals ``control_value`` within ``atol``."""
        mask = np.isclose(self.df[self.control_col].to_numpy(dtype=float), control_value, atol=atol)
        return self._subset(mask)

    def _subset(self, mask: np.ndarray) -> "FSSData":
        return FSSData(
            df=self.df[mask].copy(),
            control_col=self.control_col,
            obs_cols=dict(self.obs_cols),
            err_cols=dict(self.err_cols),
            size_col=self.size_col,
            volume_col=self.volume_col,
            sample_id_col=self.sample_id_col,
            dimension=self.dimension,
            geometry=self.geometry,
        )

    # ------------------------------------------------------------------
    # extraction for fitting
    # ------------------------------------------------------------------
    def xy(
        self, observable: str
    ) -> tuple:
        """Return (control, size, values, errors) arrays for one observable.

        Rows with missing values or non-positive errors are dropped.
        ``size`` is None when the data have no size column.
        """
        if observable not in self.obs_cols:
            raise FSSDataError(
                f"unknown observable {observable!r}; have {list(self.obs_cols)}"
            )
        col = self.obs_cols[observable]
        errcol = self.err_cols.get(observable)
        control = self.df[self.control_col].to_numpy(dtype=float)
        values = self.df[col].to_numpy(dtype=float)
        if errcol is not None:
            errors = self.df[errcol].to_numpy(dtype=float)
        else:
            errors = np.full_like(values, np.nan)
        if self.size_col is not None:
            size = self.df[self.size_col].to_numpy(dtype=float)
        elif self.volume_col is not None:
            size = self.df[self.volume_col].to_numpy(dtype=float)
        else:
            size = None

        keep = np.isfinite(values)
        if size is not None:
            keep &= np.isfinite(size)
        keep &= np.isfinite(control)
        keep &= np.isfinite(errors)
        if errcol is not None:
            keep &= errors > 0
        return control[keep], (size[keep] if size is not None else None), values[keep], errors[keep]

    # ------------------------------------------------------------------
    # reporting
    # ------------------------------------------------------------------
    def info(self) -> dict:
        sizes = self.sizes()
        return {
            "n_rows": len(self.df),
            "n_sizes": int(len(sizes)),
            "size_min": float(sizes.min()) if len(sizes) else None,
            "size_max": float(sizes.max()) if len(sizes) else None,
            "sizes": sizes.tolist(),
            "size_col": self.size_col,
            "volume_col": self.volume_col,
            "control_col": self.control_col,
            "control_range": self.control_range(),
            "observables": dict(self.obs_cols),
            "errors": dict(self.err_cols),
            "with_errors": bool(self.err_cols),
            "sample_id_col": self.sample_id_col,
            "dimension": self.dimension,
            "geometry": self.geometry,
            "has_missing": bool(self.df.isna().any().any()),
            "n_control_values": int(self.df[self.control_col].nunique()),
            "volume_from_size": None if self.volume_from_size() is None else "L^d (dimension given)",
        }


def read_table(path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"data file not found: {path}")
    sep = "\t" if path.suffix in (".tsv", ".tab") else None
    if sep is None:
        # sniff: prefer comma, fall back to whitespace/tab
        sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
        sep = "," if "," in sample.splitlines()[0] else None
    try:
        df = pd.read_csv(path, sep=sep)
    except Exception as exc:  # pragma: no cover - defensive
        raise FSSDataError(f"failed to read {path}: {exc}") from exc
    if df is None or df.empty:
        raise FSSDataError(f"empty table: {path}")
    return df


def load_fss_csv(path, **kwargs) -> FSSData:
    return FSSData.load(path, **kwargs)


def fss_from_frame(df: pd.DataFrame, **kwargs) -> FSSData:
    return FSSData.from_frame(df, **kwargs)
