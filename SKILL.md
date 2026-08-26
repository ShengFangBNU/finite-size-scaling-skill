---
name: finite-size-scaling
description: General finite-size scaling analysis for continuous transitions, BKT transitions, logarithmically corrected scaling, crossover phenomena, geometric criticality, distribution scaling, percolation, Ising/Potts/O(N), FK clusters, self-avoiding walks, and related critical systems. Use for critical-point estimation, scaling dimensions, exponent fitting, dimensionless crossings, derivative scaling, corrections to scaling, logarithmic corrections, crossover analysis, distribution collapse, finite-size stability tests, uncertainty quantification, model comparison, and simulation design.
---

# Finite-Size Scaling

Use finite-size scaling as a hypothesis-testing framework, not merely as a data-collapse procedure.

## Core workflow

1. Inspect the data and identify the finite-size variable (`L`, `V`, or another size variable) without silently converting between them.
2. Classify the observable: dimensionless, thermodynamic, derivative, geometric, distributional, or composite.
3. Identify plausible transition/RG structures:
   - ordinary continuous;
   - BKT / essential singularity;
   - marginal or logarithmically corrected;
   - crossover / multi-field;
   - special cases such as upper-critical-dimensional or dangerous-irrelevant-variable scaling.
4. Build the simplest physically justified scaling ansatz.
5. Estimate the critical point independently whenever possible.
6. Fit critical exponents or scaling dimensions.
7. Perform mandatory `L_min` stability analysis.
8. Perform fitting-window stability analysis for near-critical fits.
9. Test corrections to scaling and analytic backgrounds.
10. Check parameter identifiability and correlations.
11. Propagate uncertainty, using bootstrap/jackknife/covariance information when available.
12. Compare competing physically justified ansätze.
13. Use data collapse only as a consistency check.
14. Report statistical and systematic uncertainties separately when possible.
15. If ambiguity remains, recommend the next simulations that would best discriminate competing hypotheses.

## General scaling-field viewpoint

Start from a generic structure such as

\[
O =
O_{\rm reg}
+
L^{y_O}
(\ln L)^{\hat y_O}
\mathcal F(
u_t L^{y_t}(\ln L)^{\hat y_t},
u_g L^{y_g},
u_1 L^{y_1},
u_2 L^{y_2},
\ldots
).
\]

This is a conceptual framework. Do not fit the full expression unless the data justify it.

### Ordinary continuous transition

If \(\xi\sim |t|^{-\nu}\), use

\[
x=u_t L^{y_t}, \qquad y_t=1/\nu.
\]

### BKT / essential singularity

If

\[
\xi \sim \exp(b t^{-\sigma}),
\]

use an essential-singularity finite-size variable such as

\[
x=t[\ln(L/L_0)]^{1/\sigma}.
\]

Do not force a finite \(\nu\).

### Logarithmically corrected scaling

Allow

\[
O(L)\sim L^{y_O}(\ln L+c)^{\hat y_O}
\]

and logarithmically modified thermal fields when justified.

### Crossover scaling

If an additional parameter \(g\) moves the system between regimes/fixed points, test

\[
O(t,g,L)=L^{y_O}F(tL^{y_t},gL^{y_g}).
\]

Do not interpret drifting effective exponents as a new universality class before testing crossover.

## Observable classes

### Dimensionless observables

Examples: Binder ratios, wrapping probabilities, correlation-length ratios.

Near an ordinary critical point, use an expansion such as

\[
R(t,L)=R^*+\sum_k a_k t^k L^{k y_t}+\text{corrections}.
\]

### Derivative observables

For a dimensionless ratio,

\[
\partial_t R|_{t_c}\sim L^{y_t}.
\]

Prefer covariance-based derivatives when raw Monte Carlo information permits them. If the control parameter is a bond/edge probability \(p\) and the bond count \(N_b\) is sampled, the linear-response estimator is

\[
g=\mathrm{cov}(R,N_b)=p(1-p)\,\partial_p R,
\]

which obeys \(g|_{p_c}\sim L^{y_t}\). Extract exponents from the derivative observable directly; never obtain \(y_t\) by differentiating a fitted \(R(p)\) curve.

### Crossings

Intersections of \(R(t,L)\) curves for different \(L\) estimate the critical point only when the linear amplitude \(a_1\neq 0\). Some observables have a vanishing linear amplitude (e.g. forced by a self-duality), which makes their crossings useless for locating \(t_c\). Check the fitted leading slope/amplitude before relying on crossings, and cross-validate with at least one other observable.

### Scaling observables

At criticality,

\[
O(L)\sim L^{y_O}.
\]

### Geometric observables

Examples include largest-cluster mass, cluster radius, shortest path, backbone mass, hull length, loop length.

### Distribution observables

Examples:

\[
P(X,L)=L^{-y_X}\widetilde P(X/L^{y_X}),
\]

\[
n(s,L)=s^{-\tau}\widetilde n(s/L^{d_F}).
\]

Do not infer a power law from a visually straight log-log segment alone.

### Composite observables

For differences, ratios, duality combinations, or cancellations of observables, derive the expected scaling expansion before assigning a new exponent.

## Correction library

Consider only physically justified corrections:

- irrelevant fields \(L^{y_i}\), \(y_i<0\);
- multiple irrelevant fields;
- analytic corrections such as \(L^{-1}\), \(L^{-2}\);
- analytic backgrounds;
- mixed thermal-irrelevant terms such as \(tL^{y_t+y_i}\);
- multiplicative logarithms;
- terms such as \(L^{y_i}\ln L\) when scaling-field degeneracy or theory justifies them.

A universal correction exponent and its observable-dependent amplitude are different concepts.

## Mandatory fitting discipline

For every important fit:

- report parameter estimates and uncertainties;
- report \(\chi^2\), degrees of freedom, \(\chi^2/\mathrm{dof}\), and p-value when appropriate;
- record `L_min`, `L_max`, and the control-parameter fitting window;
- inspect residuals;
- scan over `L_min`;
- scan over fitting window where relevant;
- inspect parameter correlation / identifiability;
- compare alternative justified ansätze.

Do not prefer a fit solely because it has the smallest nominal error bar.

The standard `L_min` selection rule is: choose the smallest `L_min` for which \(\chi^2/\mathrm{dof}\approx 1\), and such that increasing `L_min` further does not lower \(\chi^2\) by much more than about one unit per degree of freedom. Prefer the smallest such `L_min` (maximal data retention) unless a robustness argument requires a larger cutoff.

## Parameter-identifiability guard

If parameters are strongly correlated or multiple parameter combinations give nearly identical fit quality:

- reduce the model;
- fix theoretically well-established nuisance parameters and rerun;
- perform profile scans;
- report the degeneracy explicitly.

Use a staged fitting protocol for correction-heavy fits: fit with all candidate parameters free, identify parameters consistent with zero or with their theoretical values, pin them, and refit; compare estimates across stages. Watch for specific known degeneracies — e.g. a shifted-log constant \((\ln L+c)\) competing with an amplitude, or two exponents competing with each other; pin one member of the pair and report the correlation. When universal constants are pinned at known values, propagate the uncertainty of the pinned values into the final error.

Never report artificial precision from an overparameterized fit.

## Model comparison

Possible competing scenarios include:

- pure power law;
- ordinary correction-to-scaling;
- logarithmically corrected scaling;
- BKT;
- multi-field crossover.

Use goodness of fit, residuals, stability, identifiability, information criteria where meaningful, and physical consistency.

If available system sizes do not distinguish the scenarios, say so.

## Fit red flags

- An effective exponent outside its physically allowed range (e.g. a leading thermal exponent larger than the space dimension) usually means the ansatz is incomplete, not that the exponent is novel.
- Strong degeneracy between a shifted-log constant \((\ln L+c)\) and an amplitude, or between two exponents, means the data cannot separate the two effects; pin one and report the correlation.
- A fit requiring many more parameters than informative data points is overparameterized; reduce and compare.
- Sign changes or strong drift of effective exponents with \(L\) may signal crossover, a missing correction, or the wrong fixed point — investigate before concluding.

## Data collapse

Use collapse only after independent or constrained estimates of the critical parameters exist.

Quantify collapse quality whenever possible.

Never infer an exponent solely from visual overlap.

## Reporting

Every final report must contain:

- data and observable classification;
- explicit scaling hypotheses;
- critical-point method;
- exponent/scaling-dimension estimates;
- corrections tested;
- `L_min` and fitting-window stability;
- parameter-identifiability assessment;
- statistical uncertainty;
- systematic uncertainty;
- competing-hypothesis comparison;
- collapse diagnostics;
- cross-observable consistency;
- remaining ambiguity;
- recommended additional simulations.

Classify conclusions as:

- **Robust**
- **Suggestive**
- **Unresolved**

## Reference routing

Read only the reference modules needed for the problem:

- `references/scaling-field-framework.md`
- `references/ordinary-continuous-fss.md`
- `references/dimensionless-crossings.md`
- `references/corrections-to-scaling.md`
- `references/analytic-backgrounds.md`
- `references/logarithmic-fss.md`
- `references/bkt-fss.md`
- `references/crossover-fss.md`
- `references/geometric-fss.md`
- `references/distribution-fss.md`
- `references/observable-selection.md`
- `references/fitting-and-systematics.md`
- `references/experimental-design.md`
- `references/special-cases/`

## Non-negotiable rules

- Never infer critical exponents solely from visual collapse.
- Never choose an ansatz because it agrees with an expected exponent.
- Never hide poor fit quality.
- Never interpret crossover immediately as a new universality class.
- Never assume every logarithm comes from an upper critical dimension.
- Never assume one correction term dominates every observable.
- Never treat strongly correlated fit parameters as independently determined.
- Never locate a critical point from crossings of an observable whose linear amplitude is zero.
- Never report an unphysical effective exponent without first questioning the ansatz.
