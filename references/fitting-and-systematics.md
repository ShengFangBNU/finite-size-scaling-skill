# Fitting and Systematics

The statistical machinery and the stability discipline every fit must pass.

## Fit results structure

Every fit returns a single structured result:

- parameters (estimates), standard errors, covariance matrix, correlation matrix;
- `χ²`, degrees of freedom, reduced `χ²`, p-value;
- AIC / AICc / BIC for model comparison;
- per-point residuals (weighted) and residual diagnostics;
- bookkeeping: `L_min`, `L_max`, fitting window in the control variable, ansatz string, pinned parameters and their pinned values;
- identifiability flags (condition number, max |correlation|).

## Statistical estimation

- Weighted nonlinear least squares with the full data covariance when available; otherwise diagonal weights `w = 1/σ²`.
- Report one-sigma margins from the covariance matrix. Bootstrap/jackknife over replicates when the data provide them, as a cross-check of the covariance errors.
- Propagate the uncertainty of pinned (fixed) quantities into the final error. If a universal constant was fixed at a value with its own error, the final error must reflect it — a fit that pins `p_c` or `Q` and quotes only the conditional error under-reports uncertainty.

## `L_min` scan

- Criterion: the smallest `L_min` for which `χ²/DF ≈ 1`, and for which raising `L_min` does not lower `χ²` by much more than ~1 unit per DF. Prefer the smallest such `L_min`.
- Report the full `L_min` trajectory of each parameter (estimate vs `L_min`), not just the accepted value, so the stability is auditable.
- A parameter that drifts monotonically with `L_min` indicates a missing correction; a parameter that jitters within noise indicates stability.

## Fitting-window scan

- For near-critical fits in the control variable `t`, repeat the fit over a set of nested windows `[-w, w]` (or `[t_min, t_max]`) around `t_c`.
- Flag windows where `t_c` or the exponents drift by more than the error bar.

## Residual diagnostics

- Plot weighted residuals vs `L` and vs `t`; check for curvature or size-trends.
- Runs test / visual randomness check for autocorrelation in `L`.
- A structured residual (e.g. remaining `L^{y_i}` trend) means the ansatz is missing a term.

## Model comparison

- Competing ansätze: pure power law, + leading correction, + second correction, + mixed term, + background, + log term. Add terms only if `χ²/DF` improves meaningfully and parameters stay stable.
- Use AIC/AICc/BIC to penalize parameter count when the improvement in `χ²` is marginal.
- Do not prefer a fit because its nominal error is smallest; prefer the fit whose parameters are stable and whose physics is justified.
- Systematic error = spread of the accepted estimate across the reasonable ansätze and across the reasonable `L_min`/window choices.

## Identifiability

- Report the correlation matrix; flag `|ρ| > 0.9` pairs.
- Report the covariance condition number; flag ill-conditioned fits.
- If two parameters are degenerate (shifted-log constant vs amplitude, or two exponents), pin one, report the correlation, and state that the data cannot separate them.
- Never report artificial precision from an overparameterized fit.

## Implementation guidance

- `fss/statistics.py`: chi-square, dof, reduced chi-square, p-value (regularized upper incomplete gamma), AIC/AICc/BIC, correlation matrix, condition number.
- `fss/diagnostics.py`: residual computation, residual diagnostics, `L_min` trajectory, window trajectory, identifiability flags.
- `scripts/stability_scan.py`: run the `L_min` and window scans and emit the trajectory tables.
- `scripts/model_compare.py`: compare candidate ansätze with the information criteria and the stability evidence.
