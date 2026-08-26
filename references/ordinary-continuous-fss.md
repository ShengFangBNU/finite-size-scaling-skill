# Ordinary Continuous FSS

The core case for Milestone 2. All other modules build on this structure.

## Assumptions

- Correlation length diverges as a power law `ξ ~ |t|^{-ν}`; scaling variable `x = t L^{y_t}` with `y_t = 1/ν > 0`.
- One relevant thermal field `t`, one (or more) irrelevant field(s) `u_i` with `y_i < 0`.
- Observable with scaling dimension `y_O` and analytic background `O_reg`:

```
O(t, L) = O_reg + L^{y_O} [ F(t L^{y_t}) + corrections ]
```

## Standard fit forms

Dimensionless observable `R` (universal at criticality):

```
R(t, L) = R^* + Σ_k a_k t^k L^{k y_t} + b_1 L^{y_i} + b_2 L^{-2} + ...
```

The `k = 1` term gives the crossing behavior; the `k = 2` term gives the curvature in the intersection region. This is the exact structure used in Wang (2013), Hou (2019).

Scaling (thermodynamic) observable at criticality:

```
O(L) = L^{y_O} (a_0 + b_1 L^{y_i} + b_2 L^{-2} + ...)
```

The exponent `y_O` is estimated by fitting the full corrected form, not by a bare two-parameter log-log slope.

Derivative (covariance) estimator for a dimensionless ratio with control `p` and sampled bond count `N_b`:

```
g = cov(R, N_b) = p(1-p) ∂_p R,   g|_{p_c} ~ L^{y_t}
```

## Procedure

1. Classify observables (dimensionless / scaling / derivative / geometric).
2. Locate `t_c` from dimensionless crossings + fits (see `dimensionless-crossings.md`).
3. Estimate `y_t` from the derivative observable and/or from dimensionless fits.
4. Estimate `y_O` from corrected power-law fits.
5. `L_min` and window stability scans (see `fitting-and-systematics.md`).
6. Corrections + background tests (see `corrections-to-scaling.md`, `analytic-backgrounds.md`).
7. Report with statistical + systematic errors.

## Diagnostics

- `χ²/DF ≈ 1` after the `L_min` cutoff; stable estimates as `L_min` grows.
- Consistency of `t_c` across observables and across ansätze.
- Residuals without systematic curvature or size-trend.
- Agreement between derivative-estimated `y_t` and fit-estimated `y_t`.

## Failure modes

- Omitting the `L^{y_i}` correction when the leading irrelevant exponent is weak (small `|y_i|`).
- Omitting the analytic background for specific-heat-like or ratio-like observables.
- Estimating `y_t` by differentiating a fitted `R(p)` instead of using the covariance estimator.
- Relying on a single observable for `t_c`.

## Implementation guidance

- `fss/fitting.py`: weighted nonlinear least squares over `(t, L)` data with `L_min`/`L_max`/window cuts, per-parameter pinning, covariance + correlation output.
- `fss/models.py`: ansatz library with the forms above; exponents may be free or pinned.
- `scripts/critical_power_fit.py`, `scripts/crossing.py`, `scripts/derivative_scaling.py`, `scripts/correction_fit.py`, `scripts/effective_exponent.py`, `scripts/stability_scan.py`.
