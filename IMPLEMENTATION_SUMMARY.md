# Implementation Summary

**Package**: `fss` (general finite-size-scaling toolkit) — v0.1.0
**Status**: Milestones 1–2 complete. The ordinary-FSS core passes its deterministic synthetic validation suite: **58/58 tests, 0 failures, 0 warnings**.
**Milestone 3** (BKT, logarithmic FSS, crossover, distribution modules) is explicitly *not started*, per the RESEARCH_PLAN milestone constraint — it begins only once the ordinary-FSS core passes.

---

## What was built

### 1. Canonical data I/O — `fss/io.py`
`FSSData` is the single interchange type: a pandas frame with `control`, `size`, and one or more `obs_*` / `err_*` columns.
- `read_csv` / `read_dataframe` with **column guessing** (observable/error columns auto-matched; `N` is never misread as a size axis).
- Selection helpers: `select_control`, `select_sizes`, `control_range`, `xy(obs)`.
- Missing/duplicate/typo robust; keeps the full frame so downstream windows can re-subset.

### 2. Common fit-result structure — `fss/fitting.py`, `fss/statistics.py`
Every fit — power-law, correction, dimensionless near-critical, crossing, derivative — returns the same `FitResult`:
- parameter values + standard errors, covariance (`inv(JᵀJ)`, honest *unscaled*), correlation matrix;
- `chi2`, `dof`, `chi2_reduced`, p-value, AIC/AICc/BIC (chi2-scale);
- weighted and unweighted residuals;
- bookkeeping (`lmin`, `lmax`, `control_window`, ansatz description);
- identifiability flags: condition number, max |correlation|, and automatic warnings (non-convergence, cond > 1e8, |corr| > 0.99, reduced chi2 > 3 or < 0.3).

**Fitting engine** (`fit_spec`): scipy `least_squares`, `trf`, `x_scale="jac"` — essential when parameters live on different scales (`pc ~ 1e-2` alongside amplitudes ~ 1; without it `trf` trapped in a bad local minimum, pc = 0.0027 with chi2/dof = 319 vs. the correct 0.015 / 0.99). `lmin`/`lmax`/`control_window` genuinely *subset* the data so stability scans re-fit, not just re-report.

### 3. Ansatz library — `fss/models.py`
`ModelSpec`/`ParamSpec` with free-or-pinned exponents and amplitudes:
- `critical_power_spec` — `O = a L^y`
- `critical_power_correction_spec` — `O = a L^y + b L^(y-omega)`
- `scaling_observable_spec` — `O = c0 + L^y(a0 + Σ b_i L^y_i)` (corrections + analytic background)
- `dimensionless_near_critical_spec` — `R(t,L) = Rc + Σ a_k t^k L^{k·yt} + corrections + mixed terms`, with `yt_fixed`/`pc_fixed` pins (see §8)
- `crossing_spec` — `t_x(L) = t_c + a L^{-λ}`

### 4. Dimensionless-observable fitting
`fit_dimensionless` with polynomial degree in the thermal field, correction exponents `(y_i, mode)` (fixed/free), and mixed thermal×irrelevant terms.

### 5. Crossing analysis — `fss/crossing.py`
- Dense-grid piecewise-linear root finding (4001 points) so the crossing of two *curved* R(t) lines is not biased by the sparse control mesh.
- **Touch points** (curves meet but do not cross, e.g. a1 = 0) are flagged `flat` — never silently trusted.
- **Multiple roots** (quadratic-in-t observables cross twice) disambiguated by seeding each pair with the previous (smaller-L) crossing — crossings converge to t_c as L → ∞.
- `crossing_fit` fits `t_x = t_c + a L^{-λ}`; **fixed-ratio size pairs** are required for a clean power-law trajectory (varying ratios make the prefactor size-dependent and wash out the law).
- `linear_amplitude_check` — verifies a1 ≠ 0 before crossings are trusted (an observable whose leading field coefficient vanishes does *not* have crossings that converge).

### 6. Derivative scaling — `fss/derivative.py`
- `control_derivative`: symmetric central difference at control grid points; one-sided at the boundary; and — new this session — a **parabola through the three bracketing points** for off-grid `at` (a plain secant returns the slope at the segment *midpoint*, biasing y_t when R(t) is curved there).
- `covariance_estimator`: the Wang2013 bond-percolation `g = cov(R, N_b)` estimator.
- `derivative_scaling_fit`: `g ~ a L^{y_t}`.

### 7. Critical power-law fitting, corrections, analytic backgrounds
`fit_critical_power(correction=..., omega=..., omega_fixed=...)` and `fit_scaling_observable(correction_exponents=..., background=...)`, following the staged protocol in the paper notes (fit free → pin consistent-with-zero/theoretical params → refit).

### 8. L_min and fitting-window stability scans — `fss/diagnostics.py`
- `lmin_scan` / `window_scan` over genuinely re-fit subsets; `scan_table` formatting.
- **L_min selection rule**: smallest L_min with reduced chi2 ≈ 1, no later drop > ~1 unit/DF.
- `effective_exponent_series` / `effective_exponent_ratios` for detecting correction drift.
- `identifiability_report` and `residual_diagnostics` (runs test + Pearson correlations, with NaN-returning guarded correlation so zero-variance residuals do not emit warnings).

### 9. Deterministic synthetic data — `fss/synthetic.py`
Fixed-seed generators (`pure_power_data`, `correction_power_data`, `dimensionless_data`) producing `FSSData` frames, so every test is reproducible.

---

## Test suite (58 tests)

| File | Coverage |
|---|---|
| `test_io.py` | schema, column guessing, selections, N-not-size |
| `test_statistics.py` | chi2 p-value, AIC/BIC, effective-exponent pair, runs test, correlation-from-cov |
| `test_models.py` | spec construction, pinning, free/pinned bookkeeping |
| `test_fitting.py` | exact recovery, lmin subsetting, correction fit, result structure, ill-conditioning warnings |
| `test_crossing.py` | crossings converge to pc, a1=0 → flat touch points, crossing trajectory fit, amplitude check |
| `test_derivative.py` | grid-point + off-grid derivative, covariance estimator, edge handling |
| `test_diagnostics.py` | effective exponents, residual diagnostics, identifiability report |
| `test_stability.py` | L_min trajectory, L_min rule selects good cut, window scan pc stability |
| `test_end_to_end.py` | full dimensionless pipeline (fit → crossings → derivative → L_min scan), critical-power pipeline, CLI script smoke tests |

---

## Bugs found and fixed during validation

1. **synthetic `dimensionless_data` swapped control/size** (root cause of ~7 failures) — sizes were written into the control column, producing negative sizes and NaN residuals.
2. **`trf` trapped in a bad local minimum for pc** — fixed with `x_scale="jac"`.
3. **Off-grid derivative evaluated at the segment midpoint, not at `at`** — gave y = 1.59 instead of 1.0; fixed with a quadratic-interpolation derivative (exact for a quadratic-in-t observable).
4. **Crossing roots biased on a coarse control grid** — fixed with dense-grid interpolation; tests additionally use a fine control mesh.
5. **Touch-point crossings not flagged flat** — exact-zero-with-no-sign-change roots now get `terr = NaN` → flagged `flat`.
6. **Crossing trajectory fit with varying-ratio pairs** gave an unphysical λ = 4.3 — switched to fixed-ratio sizes (the trajectory is a clean `t_c + a L^{-λ}` only then), plus previous-crossing seeding to disambiguate the second (spurious) root of quadratic observables.
7. **`linear_amplitude_check` ill-conditioned on a single size** — with yt free, a1 was degenerate (a1 = 0.385, se = 52); pinning yt alone still left pc leaking into a1 through the a2 cross term (a1 = 0.485 = A1/L + 2·a2·pc·L). Both `yt` and `pc` are now pinned in one-size checks.
8. **Zero-variance residuals triggered `np.corrcoef` divide-by-zero warnings** — replaced with a guarded Pearson helper.
9. **`guess_columns` skipped a valid observable** when its name collided with an errors-dict key.

---

## Files

```
fss/                io, statistics, models, fitting, diagnostics, synthetic,
                    crossing, derivative, collapse, plotting, __init__
scripts/            inspect_data, critical_power_fit, correction_fit, crossing,
                    derivative_scaling, effective_exponent, stability_scan,
                    model_compare  (all functional)
scripts/            bkt_fit, bootstrap, crossover_fit, distribution_collapse,
                    log_corrected_fit  (Milestone-3 stubs, raise SystemExit)
tests/              9 test files, 58 tests
references/         framework + ordinary-FSS docs + 5 paper notes
SKILL.md            the skill document
```

## Next (Milestone 3 — not started)

BKT / logarithmic-FSS / crossover / distribution modules, plus the stub scripts. Blocked on nothing in the ordinary core; begins only after the current green suite (per the milestone constraint).
