# Logarithmic finite-size scaling of the 4D self-avoiding walk

> Distilled from: Fang, Hou, Wang, Deng, "Logarithmic finite-size scaling of the four-dimensional self-avoiding walk," arXiv:2103.04340.
> This paper is the reference for **logarithmic FSS** (upper critical dimension, d_c = 4 for SAW). For Milestone 2 it contributes the discipline of *staged fitting with fixed exponents* and the *parameter-degeneracy warning*; its scaling forms are NOT part of the ordinary-FSS core.

## Physical problem
- Self-avoiding walk (SAW) in d = 4, the upper critical dimension, where mean-field exponents (ν = 1/2, η = 0) hold but are dressed by multiplicative logarithmic corrections. Determine the connective constant z_c and the logarithmic exponents ŷ_t, ŷ_h = 1/4.

## Main observables
- Dimensionless connectivity per monomer ξ_u/L (ξ_u the end-to-end distance), used as the dimensionless "R-like" observable whose crossings locate z_c.
- R_e² (mean square end-to-end distance) and R_g² (radius of gyration): O(L) ~ L² (ln L)^{2ŷ_h} — the "susceptibility-like" observables from which ŷ_h is extracted.
- Specific-heat-like quantity C ~ (ln L)^{2ŷ_t} (their Eq. ~27 area): pure log power law with NO power-law prefactor.

## Scaling ansatz
- Key form (their Eq. 19):
  ξ_u / [L (ln L + c_0)^{1/4}] = Σ_{k=0..m} a_k (z − z_c)^k [ L^{y_t} (ln L + c'_0)^{ŷ_t} ]^k
    + b_1 L^{y_1} + b_2 L^{y_2} + c_1 (z − z_c) L^{y_t + y_1} (ln L + c'_0)^{ŷ_t}.
  Here y_t = 1/ν = 2, ŷ_t = 1/4, y_1 = −1, y_2 = −2.
- The scaling argument itself carries a **shifted log**: x = (z − z_c) L^{y_t} (ln L + c'_0)^{ŷ_t}.
- Susceptibility-like: χ ~ L² (ln L)^{2ŷ_h}; the multiplicative log exponent doubles because χ ~ R² (second moment).
- Amplitude fit: amplitude A of the L^{y_A} (ln L)^{ŷ_A} term fitted with corrections b_1 L^{y_1} + b_2 L^{y_2}.

## Correction terms
- Power corrections b_1 L^{−1} + b_2 L^{−2}.
- Mixed term c_1 (z − z_c) L^{y_t + y_1} (ln L + c'_0)^{ŷ_t} — the same mixed t·L^{y_t+y_i} structure as in the FK-Ising note.
- Shifted logs (ln L + c_0) and (ln L + c'_0) with nonuniversal constants c_0, c'_0 — these are the analytic part of the log-correction factors.

## Fitting strategy
- Staged fitting is essential: if y_t and ŷ_t are both left free the fit is unstable. Fix y_t at its theoretical value (2) and fit ŷ_t; equivalently fix ŷ_t and fit y_t. Never fit both from the same data set without an independent constraint.
- c_0 and c'_0 cannot both be left free: they are strongly degenerate with the amplitudes (a_1 vs c_0 have almost the same effect on the fit). Fix one (set c'_0 = 0 or fix c_0) and fit the other.
- Use ξ_u/L (dimensionless) for z_c, and the power observables (R_e², R_g²) for the log exponents — a division of labor by observable type, exactly as the skill recommends.
- L_min selection: same rule — smallest L_min with χ²/DF ≈ 1; prefer the smallest such L_min.
- Final value and systematic error from comparing estimates over various reasonable ansätze (different m, different inclusion of b₂, c₀ vs c′₀ choices).

## Stability tests
- Vary L_min; check χ²/DF ≈ 1 and stability of the estimate.
- Vary polynomial order m of the a_k series.
- Compare fits with y_t free vs fixed, with each shifted-log constant free vs fixed — the spread over these defines systematic error.

## Error estimation
- Statistical from least-squares; **systematic dominates** and comes from the ansatz-comparison spread. They quote the systematic error separately and fold it in.

## Important methodological lessons
1. Parameter degeneracy is real and quantitative: c_0 (or c′_0) and a_1 produce almost identical changes in the fit. A fit that appears excellent can be hiding a degenerate direction — report the correlation and pin one of the degenerate parameters.
2. Leave-at-most-one of each competing pair free: (y_t, ŷ_t) and (c_0, c′_0). Staged fitting (fix the RG prediction, fit the rest) is the standard remedy.
3. Multiplicative logarithmic factors come with shifted logs (ln L + c_0); setting c_0 = 0 silently is a common source of systematic bias.
4. When the exponent of a log-power is requested, use the appropriate second-moment observable (χ, R²), not the dimensionless ratio.

## Failure modes
- Simultaneously fitting y_t and ŷ_t, or both shifted-log constants, → unstable, degenerate, misleading error bars.
- Extracting a log exponent from a dimensionless ratio that has had the leading power divided out with a wrong c_0 → the residual carries a spurious ln-L trend.

## Generalizable rules for the FSS skill
- Encode a **degeneracy warning**: flag fits where two parameters have near-identical sensitivity (high correlation) and suggest pinning one.
- The staged "fit free → pin consistent-with-zero / theoretically-known → refit" workflow is the skill's default for the correction-bearing fits.
- Keep the `log_corrected_fit` module out of Milestone 2; the ordinary core needs only the *diagnostic* for unphysical or degenerate fits.
