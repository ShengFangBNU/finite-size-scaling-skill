# Geometric properties of the FK representation of the Ising model

> Distilled from: Hou, Fang, Deng, "Geometric properties of the Fortuin-Kasteleyn representation of the 3D Ising model," arXiv:1811.03358.
> This paper is the most complete example of the **correction-to-scaling + analytic-background** machinery of ordinary continuous FSS, including the mixed t·L^{y_t+y_i} term. Directly relevant to Milestone 2.

## Physical problem
- FK (random-cluster) clusters of the 3D Ising model: their geometric observables (wrapping probabilities, cluster-size distribution moments, FK correlation length) belong to the 3D Ising universality class but have their own exponents/amplitudes.
- Goals: high-precision estimates of geometric exponents, universal wrapping probabilities, and the leading irrelevant exponent of 3D Ising, y_i ≈ −0.83.

## Main observables
- Dimensionless geometric observables: wrapping probabilities and moment ratios — locate the critical coupling K_c via intersections and fits.
- Susceptibility-like geometric quantities (FK cluster size, specific-heat-like C_e) — extract exponents via power-law fits.
- Correlation-length-like quantity ρ = ξ/L (their Eq. 8): ρ = ρ_0 + L^{y_t−d}(a + b L^{y_1}) — note the **analytic background ρ_0** plus the scaling part.

## Scaling ansatz
- Main fit form (their Eq. 5):
  O(ε, L) = O_c + Σ_k q_k ε^k L^{k y_t} + c_1 ε L^{y_t + y_i} + b_1 L^{y_i} + b_2 L^{−2} + b_3 L^{−3}.
  This is the full ordinary-FSS ansatz: thermal Taylor terms, a **mixed term** c_1 ε L^{y_t+y_i}, and two/three irrelevant corrections.
- Amplitude form (their Eq. 6): A(L) = L^{y_A}(a_0 + b_1 L^{y_1} + b_2 L^{y_2}) — fit the amplitude of the leading power, then the corrections to the amplitude.
- Specific-heat-like with analytic background (their Eq. A3): C_e = c_0 + L^{y_A}(a_0 + b_1 L^{y_1} + b_2 L^{y_2}) — an **analytic constant plus the scaled part**.

## Correction terms
- y_i (leading irrelevant exponent of 3D Ising) ≈ −0.83, from theory/RG; b_1 its amplitude.
- Analytic/subleading terms: b_2 L^{−2}, b_3 L^{−3}.
- Mixed term c_1 ε L^{y_t+y_i}: the cross product of the relevant field with the leading irrelevant field — essential when the fit window extends away from K_c.
- Analytic background c_0 for specific-heat-like quantities and ρ_0 for ξ/L-like quantities.

## Fitting strategy
- Fit with the correction exponent y_i free first to check consistency with theory, then **fix y_i** at the theoretical value (here −0.83) to reduce parameter count and stabilize the fit. This is the "fit free → then pin" workflow.
- Fit K_c simultaneously with amplitudes; also fix K_c from a dimensionless-observable fit and refit amplitudes at fixed K_c.
- Include the mixed term when data extend off-critical; drop it for a narrow window near K_c and check stability.

## Stability tests
- L_min scan with the standard χ²/DF ≈ 1 criterion.
- K_c from different observables (wrapping vs moment ratios vs ρ) must be consistent.
- y_i free vs fixed comparison.
- Window-in-ε comparison: with and without the mixed term c_1 ε L^{y_t+y_i}.

## Error estimation
- Statistical from χ² fits; systematic from ansatz variation (correction content, window, L_min). Reported final errors include the systematic spread.
- Where a parameter is pinned to a theoretical value, the error reflects only the fitted parameters; the pinning uncertainty is a separate systematic.

## Important methodological lessons
1. The mixed term c_1 ε L^{y_t + y_i} is not optional once the fit window includes more than the immediate critical neighborhood. Omitting it biases K_c and y_t.
2. Analytic backgrounds (c_0, ρ_0) must be modeled explicitly for specific-heat-like and ratio-like observables; absorbing them into the power law biases the exponent.
3. Fixing y_i at its theoretical value (−0.83 for 3D Ising) is standard practice and dramatically improves stability; a free-y_i fit is used only as a consistency check.
4. Exponents extracted from the *amplitude* of a leading power (Eq. 6) must themselves be corrected by the same irrelevant terms — the skill's "exponent must be extracted with its corrections" rule.

## Failure modes
- Dropping the mixed term off-critical → biased K_c.
- Treating a quantity with an analytic background (C_e, ξ/L) as a pure power law → wrong exponent.
- Leaving y_i free in a small-L window → unstable, spuriously large errors or false convergence.

## Generalizable rules for the FSS skill
- The `models.py` ansatz library must include: the pure thermal series, the mixed t·L^{y_t+y_i} term, multiple irrelevant terms (L^{y_i}, L^{−2}, L^{−3}), and analytic backgrounds (c_0 + L^{y_A}(…)).
- Fits should support "pin a parameter to a fixed value" and report pinned vs free parameters distinctly.
- Mandate a y_i-free/fixed consistency check in the fitting workflow.
