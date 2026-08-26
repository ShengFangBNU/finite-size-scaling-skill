# Scaling-Field Framework

Model-independent starting point for every FSS analysis. Everything else in this skill is a specialization of this structure.

## Central assumption

Near a critical fixed point the singular free-energy density obeys, for a system of linear size `L`,

```
f(t, h, u_1, u_2, ...; L) = f_a + L^{-d} f_s(t L^{y_t}, h L^{y_h}, u_1 L^{y_1}, u_2 L^{y_2}, ...),
```

where `t`, `h`, `u_i` are the relevant/irrelevant scaling fields, `y_t, y_h > 0` are relevant RG exponents and `y_i < 0` irrelevant ones, and `f_a` is the analytic (background) part. By differentiating with respect to the physical fields one obtains the FSS of every observable.

An observable `O` that scales with `L^{y_O}` then has the general form

```
O = O_reg + L^{y_O} (ln L)^{ŷ_O} F(t L^{y_t} (ln L)^{ŷ_t}, u_g L^{y_g}, u_1 L^{y_1}, ...),
```

with `O_reg` the regular/background part and `(ln L)^{...}` optional multiplicative logarithms for marginal operators. This is a conceptual scaffold, not a fit form: only fit as much of it as the data and the physics justify.

## Rules of use

- Keep the fields, not the raw variables, as the arguments of `F`: map the physical control (e.g. `p - p_c`, `β - β_c`) to `t` up to a nonuniversal constant.
- Every exponent refers to a scaling field; several observables can share one exponent with different, observable-dependent amplitudes.
- The analytic part `O_reg` must be modeled explicitly (see `analytic-backgrounds.md`); folding it into the power law biases exponents.
- Distinguish exponent sets by universality class and by observable sector (thermal vs magnetic vs geometric vs geometric-of-clusters).
- Degenerate irrelevant exponents (`y_i = y_j`) force logarithmic corrections `(ln L)L^{y_i}` (see `corrections-to-scaling.md`); do not assume such logs are always from an upper critical dimension.

## Specializations in this skill

- ordinary continuous: `x = t L^{y_t}`, `y_t = 1/ν` → `ordinary-continuous-fss.md`.
- BKT / essential singularity: `x = t [ln(L/L_0)]^{1/σ}` → `bkt-fss.md`.
- Logarithmically corrected / marginal: `(ln L + c)^{ŷ}` → `logarithmic-fss.md`.
- Multi-field / crossover: `F(t L^{y_t}, g L^{y_g})` → `crossover-fss.md`.
- Upper critical dimension / dangerous irrelevant variables → `special-cases/`.

## Diagnostics

- If a dimensionless ratio collapses only after subtracting a background, suspect `O_reg`.
- If effective exponents drift smoothly with `L`, suspect crossover or missing corrections, not automatically a new exponent.
- If a fit needs an unphysical exponent (e.g. leading thermal exponent `> d`), the ansatz is incomplete — re-derive the scaling form rather than accepting the value.

## Failure modes

- Fitting the full `F` on data that justify only the first two terms.
- Silent conversion between `L` and `V` or between fields.
- Treating the exponent of an amplitude as the exponent of the observable.

## Implementation guidance

- `fss/models.py` implements the specializations, parameterized by `y_t`, `y_i`, amplitudes, and background flags; never a hard-coded model class.
- Every fit records the ansatz string and the mapping `t ↔` control variable.
