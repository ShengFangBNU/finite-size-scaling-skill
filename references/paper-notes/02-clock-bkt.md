# 2D q-state clock model and BKT transitions

> Distilled from: Chen, Hou, Fang, Deng, "Monte Carlo study of duality and the BKT phase transitions of the two-dimensional q-state clock model in flow representations," arXiv:2205.02642.
> This paper is the **BKT / logarithmic-FSS** reference. For Milestone 2 it is read only for its methodological discipline (staged fitting, L_min selection, systematic-error-by-ansatz-comparison); its scaling forms are NOT to be implemented until the ordinary-FSS core passes validation.

## Physical problem
- 2D q-state clock model, q = 5–9, with two BKT transitions at β_c1 and β_c2 (QLRO phase between them). Worm algorithm in high-T and low-T flow representations.
- Goals: precise β_c1, β_c2; anomalous dimensions η(β_c1) = 1/4 and η(β_c2) = 4/q²; duality between the two critical points; self-dual point β_sd.

## Main observables
- Susceptibility-like quantities χ_h (high-T flow, equals the spin susceptibility) and χ_l (low-T flow), sampled as worm returning times; χ_diff = χ_h − χ_l; correlation length ξ = ∫r g(r)dr / ∫g(r)dr from defect separation.

## Scaling ansatz
- BKT scaling is *exponentially divergent* correlation length: ξ ~ exp(b/√t) with t = (β_c1−β)/β_c1. Verified by semilog plots vs b/√t and a data collapse of ξ/L vs b√t/(ln(L/L₀))².
- Susceptibility scaling (their Eq. 31):
  χ(β,L) = L^{7/4} (ln L + C_1)^{1/8} [ a_0 + Σ_{k=1..3} a_k ε^k (ln L + C_2)^{2k} + d_1 L^{y_1} + d_2 L^{y_2} + n_0 ε + n_1 ε² (ln L + C_2)² ].
  Note the **multiplicative logarithmic factor** (ln L + C_1)^{1/8} with a shifted log and a nonuniversal constant C_1; ε = β − β_c.
- Critical power law at β_c: χ(L)/L² = L^{−η} (a_0 + d_1 L^{y_1}) (their Eq. 32) to extract η.
- Self-dual point: χ_diff = a_0 + a_1 (β − β_sd) L^{y_d} (Eq. 33); correct form is a_0 + ε L^{2−η}(a_1 ln L + a_2) (Eq. 35).

## Correction terms
- Additive corrections d_1 L^{y_1} + d_2 L^{y_2} with y_2 < y_1 < 0; in practice they set y_1 = −1 and y_2 = −2.
- Shifted logarithms (ln L + C_1) and (ln L + C_2) — C_i are nonuniversal constants that must be fit or fixed; they are the hallmark of logarithmic FSS.
- n_0 ε term: asymmetry of the scaling function; n_1 ε²(ln L + C_2)²: nonlinearity of the RG-invariant field.

## Fitting strategy
- Staged fitting: first leave everything free; identify parameters consistent with 0 (C_2, n_1, d_2 for χ_h); set them to 0 and re-fit; then also try leaving a single correction (e.g. n_0) free. Estimates from the different stages are compared.
- Fix correction exponents (y_1 = −1, y_2 = −2) rather than fitting them.
- Fix the shifted-log constant C_1 (or fit and fix its stable value) since it competes with amplitudes.
- For χ_l at β_c2, it is preferred over χ_h because it has weaker finite-size corrections — observable selection again.
- Final error is dominated by the spread over reasonable ansätze.

## Stability tests
- L_min scan with the same rule: smallest L_min with χ²/DF ≈ 1; no drop by vastly more than one unit per DF when L_min is increased.
- Same ansatz at several L_min (e.g. L_m = 24/32/48) to show stability of β_c.
- Leave parameters free vs set-to-zero comparison is itself a stability test.

## Error estimation
- One-sigma margins from the fits; systematic error obtained by comparing estimates from various reasonable fitting Ansätze (stated almost verbatim in the paper). This "systematic from ansatz comparison" is the skill's core error philosophy.

## Important methodological lessons
1. BKT/logarithmic scaling demands shifted logs (ln L + C); C is nonuniversal and must be fit or fixed, not set to 0 silently.
2. Staged fitting — free first, pin consistent-with-zero parameters, refit — is a robust workflow for models with many competing correction terms.
3. When a leading exponent comes out *unphysical* (here y_d ≈ 2.1 > d = 2 for χ_diff), treat it as a red flag that the ansatz is incomplete, not as a real result. The correct ansatz contained an L^{2−η} ln L term.

## Failure modes
- Using intersections of R_1-like observables (with zero linear amplitude) to locate p_c — covered in the equivalent-neighbor note; here the analogous pitfall is relying on a single observable.
- Overfitting: too many correction terms with too few L values → parameters consistent with 0 → pin them.
- Comparing across methods (MC vs TN) requires a stringent consistency criterion (they use |β₁−β₂| > 3σ₁ + 3σ₂), not overlap-of-errorbars.

## Generalizable rules for the FSS skill
- Implement `log_corrected_fit` only in the later milestone; but already encode the *diagnostic* that an effective exponent larger than d (or otherwise unphysical) flags a wrong ansatz.
- The generic machinery (dimensionless-observable fits, L_min scans, χ²/DF, systematic error from ansatz comparison) is identical to ordinary FSS; only the scaling forms differ.
