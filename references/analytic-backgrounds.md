# Analytic Backgrounds

The analytic (regular) part of an observable is not noise; it must be modeled explicitly when it is not negligible.

## Where backgrounds appear

A specific-heat-like, magnetization-like, or ratio-like observable can be written as

```
O(t, L) = O_reg(t) + O_sing(t, L)
```

with `O_reg` an analytic function of `t` (a constant plus powers). Common explicit forms from the reference papers:

- Constant background: `O(L) = c_0 + L^{y_A}(a_0 + b_1 L^{y_1} + b_2 L^{y_2})` (Hou 2019 Eq. A3 for `C_e`).
- Ratio-like correlation-length observable: `ρ = ρ_0 + L^{y_t - d}(a + b L^{y_1})` (Hou 2019 Eq. 8) — a constant `ρ_0` plus the scaling part.
- Binder-ratio background corrections: `b_3 L^{d - 2 y_h}` and `b_4 L^{y_t - 2 y_h}` (Ouyang 2018 Eq. 27), the latter arising because the temperature field contains a `ρ H²` term (field-squared dependence). For 2D percolation these are `L^{-43/24}` and `L^{-5/48}`.

## Rules

- The background is analytic in the control variable, so it does NOT scale with `L`; it survives in the `L → ∞` limit only for quantities whose singular part vanishes or for specific-heat-like quantities.
- For a dimensionless ratio whose singular part is universal, the background shows up as the `L`-independent offset that corrections approach.
- For a scaling observable, folding the background into the `L^{y_O}` power law biases the exponent, typically toward smaller `|y_O|`.
- Backgrounds couple to magnetic-field-squared terms for moment-ratio observables (`b_3, b_4` above) — do not omit them in high-precision Binder fits.

## Fitting strategy

- Include the background only when it is measurable: fit with and without `c_0` (or `b_3, b_4`) and keep the simpler model unless `χ²/DF` and residual structure justify the extra terms.
- `L_min` scan: a missing background is usually the cause of a residual `L`-trend that grows at small `L`.
- When fitting at fixed `t_c`, compare `c_0` estimates across `L_min`; a stable `c_0` indicates the background is modeled correctly.

## Diagnostics

- Residuals vs `L` show a flat offset that does not vanish with increasing `L_min`.
- Fit without background: `χ²/DF` systematically high, exponent biased; with background: `χ²/DF ≈ 1`, exponent stable.
- Estimated background amplitude consistent across windows and observables.

## Failure modes

- Treating a quantity with an analytic background as a pure power law → wrong exponent.
- Including a background that isn't there → overfitting, inflated error bars.
- Omitting `b_4 L^{y_t - 2 y_h}` in Binder fits → biased `t_c`/`y_t`.

## Implementation guidance

- `fss/models.py`: every scaling ansatz accepts a `background` switch: `none | constant | analytic-series`.
- The correction-fit and critical-power-fit scripts must expose the background terms as selectable.
