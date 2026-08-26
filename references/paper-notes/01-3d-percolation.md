# Bond and site percolation in three dimensions

> Distilled from: Wang, Hou, Huang, Deng, "Bond and site percolation in three dimensions," PRL 111, 240601 (2013) / arXiv:1302.0421.
> This is the cleanest demonstration of the **ordinary continuous FSS** machinery on which Milestone 2 of the skill is built.

## Physical problem
- Bond and site percolation on the simple-cubic lattice in d = 3. The parameter is the bond/site occupation probability p; the control field is ε = p − p_c.
- Goals: (i) high-precision p_c; (ii) universal dimensionless numbers (wrapping probabilities R(x), R(a); moment ratios Q1, Q2); (iii) the thermal RG exponent y_t = 1/ν.

## Main observables
- Wrapping probabilities R(x) (wrap along x only), R(a) (wrap along any direction) — dimensionless, so they scale to universal constants at criticality.
- Moment ratios Q1 = ⟨m²⟩²/⟨m⁴⟩, Q2 = ⟨m²⟩³/⟨m⁶⟩ (m = cluster density) — dimensionless, universal at p_c.
- Derivative/covariance estimators g = cov(R, N_b) = p(1−p) dR/dp, where N_b is the number of occupied bonds. g is a susceptibility-like quantity: at criticality g ~ L^{y_t}. This is the covariance estimator the skill's `derivative_scaling` module must implement.

## Scaling ansatz
- Main fit form (their Eq. 4):
  O(ε, L) = O_c + Σ_{k=1..2} q_k ε^k L^{k y_t} + b_1 L^{y_i} + b_2 L^{−2}.
  The linear term in ε gives the "crossing" behavior; the quadratic term controls the curvature of the intersection region.
- For the amplitude A of the covariance/density observable (their Eq. 5):
  A(L) = L^{y_A} (a_0 + b_1 L^{−1.2} + b_2 L^{−2}), i.e. the leading power extracted from a power-law fit, with y_i = −1.2(2) the leading irrelevant exponent of 3D percolation.
- y_t is obtained from the covariance observable: g = cov(R,N_b) ~ L^{y_t} at p_c.

## Correction terms
- Two leading corrections retained: L^{y_i} with y_i ≈ −1.2 and an analytic/analytic-extension L^{−2} term (b_2).
- Amplitudes are observable-dependent: for R(x) and R(a) the b_1 amplitude is consistent with zero, so those observables have *weaker* corrections than the moment ratios Q1, Q2. This is a quantitative example of the skill's guidance to choose observables whose leading correction amplitude vanishes.
- y_i was also treated as a fit parameter in one variant; leaving it free was unstable, so it was fixed.

## Fitting strategy
- p_c from least-squares fits of the full ansatz (not just from intersections), on a window of L and of p around the would-be p_c.
- Preferred route for y_t: fit the covariance g = p(1−p)dR/dp at p_c directly (power law in L) rather than fitting R in a neighborhood of p_c and differentiating the fit. This is a concrete, transferable lesson: **estimate exponents from the derivative observable, not by differentiating a fitted curve**.
- Also fit p_c and y_t simultaneously in multi-parameter fits; check consistency.

## Stability tests
- L_min scan: for each candidate ansatz, raise the smallest system size included until χ²/DF ~ O(1). Prefer the smallest L_min for which χ²/DF is acceptable and for which further raising L_min does not lower χ² by much more than one unit per degree of freedom.
- Ansatz comparison: fits with and without the b_2 term, with y_i free vs fixed, with different p-window widths.

## Error estimation
- Statistical: standard χ² least-squares one-sigma margins from the covariance matrix.
- Systematic: spread of estimates across reasonable ansätze (different L_min, different correction content) dominates the final error. Reported errors combine the two.

## Important methodological lessons
1. Estimate exponents from covariance/derivative observables g = cov(R,N_b) = p(1−p)dR/dp, which are linear-response estimators and are much better behaved than differentiating fits.
2. Dimensionless observables (R(x), R(a), Q1, Q2) are the workhorses: universal at p_c, so p_c can be located precisely.
3. Pick observables with vanishing leading correction amplitudes when possible; the paper shows R(x)/R(a) are cleaner than Q1/Q2.
4. The L_min selection rule "smallest L_min with χ²/DF ≈ 1, no further drop > ~1/DF" is applied exactly as the skill prescribes.

## Failure modes
- Leaving the irrelevant exponent y_i free made fits unstable → fix it to a known/theoretical value (here ≈ −1.2).
- Forcing every observable to have the same correction amplitude is wrong: amplitudes depend on the observable, so the same correction exponent can enter with different (even zero) amplitude.

## Generalizable rules for the FSS skill
- Provide `derivative_scaling` (covariance estimator g = cov(O, N_b)) as a first-class tool.
- Implement the generic ansatz O(ε,L) = O_c + Σ_k q_k ε^k L^{k y_t} + Σ_i b_i L^{y_i} with per-observable amplitudes and with L^{−2}-type analytic terms.
- Mandate L_min stability scans with the χ²/DF criterion as the primary model-selection gate, and systematic error from ansatz comparison.
