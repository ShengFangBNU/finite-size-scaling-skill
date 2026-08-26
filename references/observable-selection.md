# Observable Selection

Which observables to simulate, sample, and fit — and what each is good for.

## Classification

- **Dimensionless** (Binder ratio, wrapping probabilities, correlation-length ratio `ξ/L`, connectivity per monomer): universal at `t_c`; best for locating `t_c` and for data collapse. `R(t,L) → R^*` plus corrections.
- **Scaling / thermodynamic** (susceptibility, second moments, specific heat, cluster sizes): power laws `O ~ L^{y_O}`; best for exponents. Must be corrected for background and irrelevant fields.
- **Derivative / covariance** (`g = cov(R,N_b) = p(1-p) ∂_p R`): linear-response estimators of `∂R/∂t ~ L^{y_t}`; the cleanest route to `y_t`. Never obtained by differentiating a fitted curve.
- **Geometric** (largest-cluster mass, backbone, shortest path, hull): sector-specific exponents `y_O`; same machinery, different fixed-point values.
- **Distributional**: `P(X,L) = L^{-y_X} P̃(X/L^{y_X})`, `n(s,L) = s^{-τ} ñ(s/L^{d_F})`; a power law must never be inferred from a straight log-log segment alone.

## Selection criteria

1. **Strong leading amplitude**: dimensionless observables with a large linear amplitude `a_1` give good crossings; a zero-`a_1` observable (self-duality forced, e.g. `R_1` in percolation) is useless for locating `t_c`.
2. **Weak corrections**: observables whose leading correction amplitude vanishes are preferred for `t_c` (Wang 2013: `R(x)/R(a)` cleaner than `Q1/Q2`; Chen 2022: `χ_l` cleaner than `χ_h` for `β_c2`).
3. **Sector coverage**: choose at least one dimensionless + one scaling + (when possible) one derivative observable per sector so exponents and `t_c` are independently constrained.
4. **Known universal target**: if the universal value `R^*`/`Q^*` is known exactly, fixing it turns the fit into a measurement of nonuniversal amplitudes and `t_c` (Ouyang 2018).

## Pairing rules

- Use the dimensionless observable for `t_c`; use the scaling observable for `y_O`; use the covariance for `y_t`; cross-check `y_t` from dimensionless curvature.
- For dual models (two flow/spin representations), observables on the two sides of a duality can be combined into a difference whose crossings are exceptionally sharp (Chen 2022 `χ_diff`), but derive its scaling expansion first.

## Diagnostics

- Consistency of `t_c` across ≥ 2 independent observables.
- For each observable, the fitted correction amplitudes should be small relative to the leading term; flag observables whose corrections dominate.
- Effective exponents from different observables agree within combined errors once corrections are included.

## Failure modes

- Relying on one observable for everything.
- Choosing an observable with a vanishing `a_1` for crossings.
- Mixing observables from different sectors/universality classes in one fit.

## Implementation guidance

- `fss/io.py`: data schema tags each column with an observable class so scripts can route automatically.
- `scripts/inspect_data.py`: report per-observable type, `L` coverage, and a quick estimate of `a_1` for each dimensionless observable.
