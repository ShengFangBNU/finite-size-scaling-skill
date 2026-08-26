"""General finite-size-scaling toolkit.

Milestone 2 scope (ordinary FSS core):

- canonical data I/O            -> :mod:`fss.io`
- dimensionless-observable fits -> :mod:`fss.fitting`, :mod:`fss.models`
- crossing analysis             -> :mod:`fss.crossing`
- derivative scaling            -> :mod:`fss.derivative`
- critical power-law fits       -> :mod:`fss.fitting`
- leading corrections           -> :mod:`fss.models`
- analytic backgrounds          -> :mod:`fss.models`
- stability scans & diagnostics -> :mod:`fss.diagnostics`
- common fit-result structure   -> :mod:`fss.fitting`
- deterministic synthetic data  -> :mod:`fss.synthetic`

BKT, logarithmic-FSS, crossover, and distribution modules are Milestone 3
scope and deliberately not started.
"""

from . import io
from . import statistics
from . import models
from . import fitting
from . import diagnostics
from . import collapse
from . import plotting
from . import synthetic
from . import crossing
from . import derivative

from .io import FSSData, load_fss_csv, fss_from_frame, guess_columns, FSSDataError
from .fitting import FitResult, fit_spec, fit_critical_power, fit_dimensionless, fit_scaling_observable
from .crossing import crossings, crossing_fit, linear_amplitude_check
from .derivative import control_derivative, covariance_estimator, derivative_scaling_fit
from .models import (
    critical_power_spec,
    critical_power_correction_spec,
    scaling_observable_spec,
    dimensionless_near_critical_spec,
    crossing_spec,
)

__version__ = "0.1.0"

__all__ = [
    "io", "statistics", "models", "fitting", "diagnostics", "collapse",
    "plotting", "synthetic", "crossing", "derivative",
    "FSSData", "load_fss_csv", "fss_from_frame", "guess_columns", "FSSDataError",
    "FitResult", "fit_spec", "fit_critical_power", "fit_dimensionless",
    "fit_scaling_observable",
    "crossings", "crossing_fit", "linear_amplitude_check",
    "control_derivative", "covariance_estimator", "derivative_scaling_fit",
    "critical_power_spec", "critical_power_correction_spec",
    "scaling_observable_spec", "dimensionless_near_critical_spec",
    "crossing_spec",
    "__version__",
]
