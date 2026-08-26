# Dimensionless Crossings

How to locate a critical point from dimensionless observables, and when not to.

## Definition

A dimensionless observable `R(t, L)` (Binder ratio, wrapping probability, correlation-length ratio) is universal at criticality: `R(t_c, L) → R^*` plus corrections. Curves of `R` vs the control `t` for different `L` cross near `t_c`.

Near-critical expansion:

```
R(t, L) = R^* + a_1 (t - t_c) L^{y_t} + a_2 (t - t_c)^2 L^{2 y_t} + b_1 L^{y_i} + ...
```

The pair-wise crossing `t_x(L, sL)` — the control value where curves of sizes `L` and `sL` intersect — shifts toward `t_c` as

```
t_x(L, sL) = t_c + a L^{-λ},   λ = y_t + y_i  (to leading order)
```

so a plot of `t_x` vs `L^{-λ}` (or vs `1/L` to see the trend) extrapolates `t_c`.

## Requirements

- The linear amplitude `a_1` must be nonzero. If `a_1 = 0` (e.g. forced by a self-duality symmetry), crossings do NOT converge to `t_c` in the naive way — the leading intersection displacement comes from a higher term. Check the fitted slope before trusting crossings (see `05-equivalent-neighbor-crossover.md`, where `R_1` has `a_1 = 0`).
- Several `(L, sL)` pairs, ideally with different `s`, to detect the `L^{-λ}` shift.
- Corrections (the `b_1 L^{y_i}` term) shift all crossings together; extrapolate in `L`, don't quote a single-pair crossing as `t_c`.

## Procedure

1. Plot `R` vs control for all `L`; visually confirm well-separated curves.
2. Interpolate (cubic) or fit each curve and compute `t_x(L, sL)` for each pair.
3. Check that `t_x` is monotone in `L` and extrapolate to `L → ∞` (e.g. fit `t_x = t_c + a L^{-λ}` with `λ ≈ y_t + y_i`, or use `1/L` trend).
4. Cross-validate with a full simultaneous fit `R(t,L) = R^* + Σ a_k t^k L^{k y_t} + corrections`; report the fit-based `t_c`.
5. Never rely on one observable: check consistency of `t_c` across at least two dimensionless observables.

## Diagnostics

- Crossing shifts consistent with a single `t_c` as `L → ∞`.
- Fit-based `t_c` inside the error bar of the extrapolated crossings.
- `χ²/DF ≈ 1` for the simultaneous fit.

## Failure modes

- Zero-`a_1` observable: crossings appear to exist but extrapolate to the wrong value (or do not converge).
- Corrections not accounted for: `t_x` shifts with `L`; quoting the smallest-`L` crossing biases `t_c`.
- Interpolation error when `R(t,L)` is noisy or the window is too wide.

## Implementation guidance

- `scripts/crossing.py`: accept `R(t, L)` data, interpolate per `L`, compute `t_x(L, sL)` for all pairs, report the `t_x` table and the extrapolation fit with `λ`.
- Return a warning flag when the fitted `a_1` of any observable is consistent with zero within `~2σ`.
