# Equivalent-neighbor percolation and crossover scaling

> Distilled from: Ouyang, Deng, Blöte, "Equivalent-neighbor percolation models in two dimensions: crossover between mean-field and short-range behavior," arXiv:1808.05812.
> This paper is the reference for **crossover scaling** (MF → short-range) and for **logarithmic corrections from degenerate irrelevant exponents** (the b_l ln(L) L^{y_i} terms). It also demonstrates *fixing universal constants in the fit* — a discipline Milestone 2's `correction_fit` should support.

## Physical problem
- 2D bond percolation with a finite interaction range: each site couples to z equivalent neighbors (z = 4, 8, 12, 20, 28, 36, 48, 60, 224, 1224, 4016), plus the complete-graph (mean-field) limit z → ∞.
- Questions: do all finite-z models belong to the short-range (SR) percolation universality class? What is the crossover exponent y_r (≈ 2/3) describing a finite-range perturbation at the MF fixed point? Is the crossover continuous (no tricritical point)?
- Answer: uniform crossover; finite z ⇒ SR universality (y_t = 3/4, y_h = 91/48); MF limit is an unstable fixed point with y_t = 2/3, y_h = 4/3; crossover exponent φ with φ^{−1} = y_t/y_r ≈ 1.006(7).

## Main observables
- Wrapping probabilities R_b (both directions), R_e (either direction), R_1 (one direction); Binder ratio Q (Eq. 38); density of the largest cluster c_l.
- For the complete graph, only Q (and c_l) are meaningful; use the integrated probability p_i = p·L².

## Scaling ansatz
- Free energy FSS (their Eq. 12): f(t,h,u₁,u₂,L) = f_a + L^{−d} f_s(L^{y_t}t, L^{y_h}h, L^{y_1}u₁, L^{y_2}u₂, 1). Analytic background f_a is explicit.
- Binder ratio (their Eq. 27):
  Q(t,L) = Q + Σ_k q_k t^k L^{k y_t} + b_i L^{y_i} + b_l ln(L) L^{y_i} + b_3 L^{d−2y_h} + b_4 L^{y_t−2y_h} + …
- Wrapping probabilities (their Eq. 28): R_w(t,L) = R_w + Σ_k a_k t^k L^{k y_t} + b_i L^{y_i} + b_l ln(L) L^{y_i} + …
- Concrete percolation fit form (their Eq. 39/40): corrections b_i L^{−2} + b_l ln(L) L^{−2} + b_3 L^{−3} (y_i = −2 is degenerate with the second thermal exponent, forcing the log term); c_l(p,L) = L^{y_h−2}[ Σ_k g_k(p−p_c)^k L^{k y_t} + corrections ].

## Correction terms
- Leading irrelevant terms b_i L^{y_i} and b_l ln(L) L^{y_i}: for 2D percolation y_i = −2 is degenerate, so the correction is L^{−2} together with a ln(L) L^{−2} term. This is the paper's central technical feature.
- Analytic-background terms for Q: b_3 L^{d−2y_h} (= L^{−43/24} for 2D percolation since y_h = 91/48) and b_4 L^{y_t−2y_h} (= L^{−5/48}), the latter arising because the temperature field contains a ρH² term.
- Subleading L^{−3} term. In the SR class the amplitudes are approximately proportional to z² (distance from the SR fixed point); they change sign between z = 4 and z = 8.

## Fitting strategy
- Locate p_c from R_b and R_e (whose linear amplitudes a₁ ≠ 0). R_1 is **useless for locating p_c**: self-duality forces a₁ = 0, so its curves do not intersect cleanly. Then **fix p_c** at the wrapping-derived value in all subsequent fits (Binder, c_l).
- Fix the universal constants R_b, R_e, R_1 and Q at their known SR values (R_b = 0.351642855, R_e = 0.690473725, R_1 = 0.169415435; Q = 0.87057(2)) when testing universality and when measuring correction amplitudes; fix y_t = 3/4, y_h = 91/48 at their exact values.
- For the complete graph, fix p_{i,c} = 1 and the MF exponents y_t = −2/3, y_h = −4/3 (as irrelevant corrections), fit Q.
- Narrow p-intervals with few p-dependent terms; χ²/DF as reliability gate.

## Stability tests
- L_min scan (they tabulate L_min and χ²/DF for every model/observable).
- Universality test: R and Q estimates from fits with universal constants fixed must agree across z; y_t and y_h free fits must return the SR values.
- Crossover demonstration: effective exponents y_{t,eff}(L) (Eq. 45) and y_{h,eff} plotted vs L^{y_i}; they must converge to 3/4 (SR) or 2/3 (MF); data collapse via rescaled length L_r = L/b(z) (Eq. 43).

## Error estimation
- One-sigma margins from least-squares assuming ansatz validity; quoted errors "can be a few times larger after taking into account the errors in p_c and Q" — i.e. **propagate the uncertainty of pinned/fixed quantities** into the final error. The skill's identifiability/reporting rules must flag this.
- Systematic from ansatz variation (whether a log term is included, L_min, window).

## Important methodological lessons
1. An observable with zero linear amplitude (R_1, a₁ = 0) gives no useful crossings — always check the leading amplitude before using crossings to locate the critical point. Self-duality can force a₁ = 0.
2. Degenerate irrelevant exponents (y_i = −2 degeneracy) force *logarithmic* correction terms b_l ln(L) L^{y_i}; ignoring them biases amplitudes and exponents.
3. Fixing universal constants (R, Q, y_t, y_h) at known exact values converts the fit into a measurement of nonuniversal amplitudes and p_c — a powerful and legitimate strategy when universality is already established.
4. Effective-exponent plots vs L^{y_i} or 1/ln L are the cleanest way to exhibit crossover; convergence of y_{eff} to the target value is the validation.
5. The analytic background for Q includes b_3 L^{d−2y_h} and b_4 L^{y_t−2y_h} terms from the magnetic-field-squared dependence of the temperature field — these are non-negligible and must be modeled.

## Failure modes
- Using R_1 to locate p_c → biased/failed crossings.
- Omitting the ln(L)L^{−2} correction for a degenerate y_i → systematic amplitude errors.
- Forgetting the b_3, b_4 analytic-background terms in Q fits → biased Q and exponents.
- Comparing finite-z models with the MF model as if the same scaling form applied → the MF limit has different exponents (2/3, 4/3), and Q values are different (0.381 vs 0.8706).

## Generalizable rules for the FSS skill
- The correction library must include `ln(L)·L^{y_i}` terms selectable alongside `L^{y_i}` terms — degenerate irrelevant exponents are common.
- Crossing analysis must **report the fitted leading slope/amplitude** of each observable so a zero-amplitude observable is caught automatically.
- Fit functions must support **pinning universal constants** and **propagating their uncertainties** into final errors.
- The `crossover_fit` module is a later milestone; for now, encode the *effective-exponent* diagnostic (y_eff vs L^{y_i}) as a way to detect crossover in `effective_exponent.py`.
