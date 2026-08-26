# Corrections to Scaling

Every finite-`L` estimate needs corrections. This module defines which corrections are legitimate and how to fit them.

## Sources

1. Irrelevant scaling fields: `b_i L^{y_i}`, `y_i < 0`. Amplitudes `b_i` are observable-dependent and analytic in the scaling fields; an exponent shared by all observables, an amplitude that differs per observable.
2. Multiple irrelevant fields: several `y_i` values; only include those the data support.
3. Analytic corrections `L^{-1}, L^{-2}, ...`: analytic background of the free energy / observable expansion.
4. Analytic backgrounds `c_0` (a constant) — see `analytic-backgrounds.md`.
5. Mixed thermal-irrelevant terms `c_1 t L^{y_t + y_i}`: cross product of the relevant field with an irrelevant field; required when the fit window extends away from `t_c`.
6. Logarithmic corrections from degenerate irrelevant exponents: if `y_i = y_j` (or an exponent coincides with `-2`), the scaling fields couple and produce `b_l ln(L) L^{y_i}` terms (Ouyang 2018). Also multiplicative logs from marginal operators — but don't assume every log is from an upper critical dimension.
7. Amplitude-of-power corrections: when the leading power `L^{y_O}` is divided out, the amplitude itself carries corrections, `A(L) = L^{y_O}(a_0 + b_1 L^{y_1} + b_2 L^{y_2})` (Wang 2013, Hou 2019).

## Canonical fit forms (from the reference papers)

Dimensionless (Wang 2013 Eq. 4 / Hou 2019 Eq. 5):

```
R(t, L) = R^* + Σ_{k=1..K} q_k t^k L^{k y_t} + c_1 t L^{y_t + y_i} + b_1 L^{y_i} + b_2 L^{-2} + b_3 L^{-3}
```

Scaling amplitude (Hou 2019 Eq. 6):

```
O(L) = L^{y_O} (a_0 + b_1 L^{y_1} + b_2 L^{y_2})
```

Specific-heat-like with background (Hou 2019 Eq. A3):

```
C_e(L) = c_0 + L^{y_A}(a_0 + b_1 L^{y_1} + b_2 L^{y_2})
```

Degenerate-irrelevant case (Ouyang 2018 Eq. 39/40): corrections `b_i L^{-2} + b_l ln(L) L^{-2} + b_3 L^{-3}`.

## Fitting strategy

- Start with the leading correction only; add terms one at a time and demand a real improvement in `χ²/DF` (about one unit per DF is not improvement).
- Staged protocol: leave all candidate corrections free → pin amplitudes consistent with zero → refit → compare.
- Fix a well-established correction exponent (e.g. `y_i = -0.83` for 3D Ising, `-1`/`-2` in BKT/log contexts) rather than fitting it, to stabilize; a free-`y_i` fit is a consistency check only.
- Remember: a universal correction exponent with observable-dependent amplitude — never force all observables to share an amplitude.

## Diagnostics

- `L_min` scan: correction terms are visible as the `χ²/DF` dropping below ~1 only once `L_min` is large enough, and rising steeply when `L_min` is too small.
- Fitted correction amplitudes stable across `L_min` and window.
- Residuals free of a remaining `L^{y_i}` trend.

## Failure modes

- Omitting the mixed term `t L^{y_t + y_i}` when off-critical data are included → biased `t_c` and `y_t`.
- Omitting a degenerate `ln(L)L^{y_i}` term → biased amplitudes and, sometimes, exponents.
- Fitting too many corrections on too few `L` values → overparameterization (amplitudes consistent with zero → pin them).
- Treating a correction amplitude as universal.

## Implementation guidance

- `fss/models.py`: correction library with per-term on/off and per-term amplitude; support pinned exponents and pinned amplitudes.
- `fss/statistics.py`: chi-square, dof, reduced chi-square, p-value, AIC/AICc/BIC for model comparison.
- `scripts/correction_fit.py`: fit `O(L) = c0 + L^{y_O}(a0 + Σ b_i L^{y_i})` with selectable corrections.
