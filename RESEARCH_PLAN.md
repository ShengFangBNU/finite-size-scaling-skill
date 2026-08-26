# Research Plan: Build a General-Purpose Finite-Size Scaling Agent Skill

## 0. Mission

Build a reusable, research-grade `finite-size-scaling` Agent Skill that can be used by local coding/research agents (Codex, Claude Code, OpenCode, or equivalent) to analyze finite-size numerical data in statistical physics and critical phenomena.

The skill must be general. It must not be organized around one specific model or around high-dimensional systems only.

The intended scope includes:

- ordinary continuous phase transitions;
- percolation;
- Ising, Potts, and O(N) models;
- Fortuin–Kasteleyn cluster observables;
- self-avoiding walks;
- geometric critical phenomena;
- BKT transitions;
- logarithmically corrected scaling;
- upper-critical-dimensional special cases;
- crossover phenomena;
- multi-scaling-field problems;
- cluster-size and probability-distribution scaling.

The finished system should behave as an **RG-aware finite-size-scaling research agent**, not as a data-collapse plotting utility.

---

## 1. Core scientific principles

The implementation must enforce the following principles.

### 1.1 FSS is hypothesis testing

The system should compare physically motivated scaling hypotheses instead of optimizing a collapse first and interpreting it afterwards.

### 1.2 Separate physics from numerics

The Agent is responsible for:

- classifying the observable;
- identifying plausible scaling structures;
- choosing candidate ansätze;
- interpreting results.

Numerical code is responsible for:

- optimization;
- covariance estimation;
- bootstrap/jackknife;
- goodness-of-fit;
- stability scans;
- residuals;
- information criteria;
- collapse metrics.

### 1.3 Data collapse is secondary

Critical points and exponents should, whenever possible, be independently estimated before collapse.

### 1.4 Stability is mandatory

Every asymptotic quantity must be tested against:

- minimum size \(L_{\min}\);
- fitting window;
- correction terms;
- polynomial order;
- fixed/free nuisance parameters.

### 1.5 Statistical and systematic uncertainty are distinct

One covariance matrix from one fit is not an uncertainty analysis.

### 1.6 Ambiguity is an allowed result

If available sizes cannot distinguish, for example,

\[
L^{-\omega}
\]

from

\[
(\ln L)^{-q},
\]

the system should report the ambiguity and recommend new simulations.

---

# 2. Primary methodological literature

The first implementation stage must distill the following papers as methodology sources.

## Paper 1 — 3D percolation

**Bond and site percolation in three dimensions**

Extract:

- wrapping/dimensionless observables;
- critical-point fits;
- near-critical polynomial expansions;
- irrelevant and analytic corrections;
- derivative/covariance estimators;
- \(L_{\min}\) testing;
- goodness-of-fit practice.

## Paper 2 — BKT / clock model

**Monte Carlo study of duality and the Berezinskii-Kosterlitz-Thouless phase transitions of the two-dimensional q-state clock model in flow representations**

Extract:

- essential singularity;
- logarithmic BKT scaling fields;
- multiplicative logarithms;
- duality/composite observables;
- misleading effective power laws.

## Paper 3 — 4D SAW

**Logarithmic finite-size scaling of the self-avoiding walk at four dimensions**

Extract:

- marginal scaling;
- multiplicative logarithms;
- shifted logarithms;
- identifiability and parameter degeneracy;
- constrained vs free fits;
- clean-observable selection;
- inability of finite sizes to distinguish log vs power corrections.

## Paper 4 — FK geometry

**Geometric properties of the Fortuin-Kasteleyn representation of the Ising model**

Extract:

- geometric observables;
- fractal dimensions;
- cluster-size distributions;
- full distribution collapse;
- hyperscaling tests;
- mixed thermal–irrelevant corrections;
- analytic backgrounds.

## Paper 5 — crossover

**Equivalent-neighbor percolation models in two dimensions: Crossover between mean-field and short-range behavior**

Extract:

- general scaling-field formulation;
- multi-variable scaling;
- crossover fields and crossover lengths;
- effective-exponent flow;
- multiple irrelevant fields;
- logarithmic corrections from degenerate fields;
- distinguishing crossover from new universality.

---

# 3. Deliverables

The Agent must create the following repository structure.

```text
finite-size-scaling/
├── SKILL.md
├── README.md
├── pyproject.toml
├── RESEARCH_PLAN.md
├── IMPLEMENTATION_SUMMARY.md
├── SELF_REVIEW.md
│
├── references/
│   ├── scaling-field-framework.md
│   ├── ordinary-continuous-fss.md
│   ├── dimensionless-crossings.md
│   ├── corrections-to-scaling.md
│   ├── analytic-backgrounds.md
│   ├── logarithmic-fss.md
│   ├── bkt-fss.md
│   ├── crossover-fss.md
│   ├── geometric-fss.md
│   ├── distribution-fss.md
│   ├── observable-selection.md
│   ├── fitting-and-systematics.md
│   ├── experimental-design.md
│   ├── special-cases/
│   │   ├── upper-critical-dimension.md
│   │   ├── dangerous-irrelevant-variables.md
│   │   └── boundary-condition-effects.md
│   └── paper-notes/
│       ├── 01-3d-percolation.md
│       ├── 02-clock-bkt.md
│       ├── 03-4d-saw-log.md
│       ├── 04-fk-geometry.md
│       └── 05-equivalent-neighbor-crossover.md
│
├── fss/
│   ├── __init__.py
│   ├── io.py
│   ├── models.py
│   ├── fitting.py
│   ├── statistics.py
│   ├── diagnostics.py
│   ├── collapse.py
│   └── plotting.py
│
├── scripts/
│   ├── inspect_data.py
│   ├── crossing.py
│   ├── derivative_scaling.py
│   ├── critical_power_fit.py
│   ├── correction_fit.py
│   ├── log_corrected_fit.py
│   ├── bkt_fit.py
│   ├── crossover_fit.py
│   ├── effective_exponent.py
│   ├── distribution_collapse.py
│   ├── stability_scan.py
│   ├── bootstrap.py
│   └── model_compare.py
│
├── examples/
└── tests/
```

The Agent may improve this layout if it documents the reason.

---

# 4. Phase A — Repository and environment audit

## Tasks

1. Recursively inspect the working directory.
2. Search for:
   - local copies of the five papers;
   - existing FSS scripts;
   - existing analysis pipelines;
   - `AGENTS.md`;
   - `CLAUDE.md`;
   - `SKILL.md`;
   - environment files.
3. Record Python version and installed scientific libraries.
4. Do not modify existing research code at this stage.
5. Create a short implementation log.

## Acceptance criterion

The Agent understands the local environment before implementation starts.

---

# 5. Phase B — Literature distillation

Create one methodology note per paper.

Each note must contain:

1. physical problem;
2. observables;
3. scaling ansatz;
4. correction terms;
5. fitting workflow;
6. stability tests;
7. uncertainty treatment;
8. key methodological innovations;
9. failure modes;
10. reusable rules for the general skill.

Do not write only a conventional paper summary.

The goal is to answer:

> What should a future FSS Agent learn from this paper?

## Acceptance criterion

Every rule introduced into `SKILL.md` can be traced to either:

- general scaling theory;
- one or more distilled methodological lessons;
- explicit numerical-statistics best practice.

---

# 6. Phase C — General theory architecture

Organize the knowledge base using three orthogonal axes.

## 6.1 Transition / RG structure

Support:

- ordinary continuous;
- BKT / essential singularity;
- marginal/log-corrected;
- crossover/multi-field;
- special upper-critical-dimensional cases;
- dangerous irrelevant variables where relevant.

## 6.2 Observable structure

Support:

- dimensionless observables;
- order parameters;
- susceptibilities;
- derivatives;
- geometric quantities;
- fractal observables;
- probability distributions;
- cluster-size distributions;
- composite observables.

## 6.3 Correction structure

Support:

- irrelevant-field corrections;
- multiple irrelevant fields;
- analytic corrections;
- analytic backgrounds;
- mixed thermal–irrelevant terms;
- multiplicative logarithms;
- \(L^{y_i}\ln L\)-type terms;
- observable-dependent amplitudes.

## Acceptance criterion

The architecture is not organized simply as “2D / 3D / high-d”.

---

# 7. Phase D — Canonical data model

Design a flexible table-based schema.

Suggested canonical meanings:

```text
L
V
control
observable
error
sample_id
```

Allow user mappings such as:

```text
T -> control
beta -> control
p -> control
r -> control
N -> V
```

Important rule:

> Never silently replace \(V\) by \(L\), or assume \(V=L^d\), unless geometry and dimension are explicitly known.

Implement:

```text
fss/io.py
scripts/inspect_data.py
```

The data inspector should report:

- candidate size columns;
- parameter columns;
- observables;
- uncertainties;
- missing values;
- number of sizes;
- size range;
- parameter coverage;
- whether raw samples appear to be present.

---

# 8. Phase E — Ordinary continuous-transition core

This is the first numerical milestone.

Implement:

## 8.1 Critical power law

\[
O(L)=aL^y
\]

with optional:

\[
O(L)=aL^y(1+bL^{-\omega}).
\]

## 8.2 Dimensionless observable fits

Support a near-critical form such as

\[
R(t,L)=R_c+a_1tL^{y_t}+a_2t^2L^{2y_t}+\cdots.
\]

Allow irrelevant corrections.

## 8.3 Pairwise crossings

Estimate crossings for size pairs \((L,sL)\).

Support extrapolation:

\[
t_\times(L,sL)=t_c+aL^{-\lambda}+\cdots.
\]

## 8.4 Derivative scaling

Support:

\[
\partial_t R|_{t_c}\sim L^{y_t}.
\]

Prepare an interface for covariance-derived derivatives if raw MC data are supplied.

## Acceptance criterion

Synthetic ordinary-transition data recover known parameters within uncertainty.

---

# 9. Phase F — Corrections and backgrounds

Create reusable model components rather than separate hard-coded fit scripts.

Support:

\[
L^{-\omega},
\]

multiple irrelevant terms,

\[
L^{-1},L^{-2},
\]

regular backgrounds,

\[
O(L)=c_0+L^y(a_0+\cdots),
\]

and mixed terms such as

\[
tL^{y_t-\omega}.
\]

The code must allow physically justified combinations without producing an uncontrolled combinatorial model generator.

## Acceptance criterion

The Agent can compare a pure power fit against leading-correction fits and detect when small sizes bias the result.

---

# 10. Phase G — Mandatory stability engine

Implement `stability_scan.py`.

At minimum scan:

\[
L_{\min}.
\]

For near-critical fits also scan fitting-window width.

Each scan must record:

- model;
- \(L_{\min}\);
- \(L_{\max}\);
- parameter window;
- estimates;
- uncertainties;
- \(\chi^2/\mathrm{dof}\);
- p-value;
- number of points;
- warnings.

Generate stability plots.

## Agent interpretation rule

Prefer an asymptotic plateau supported by acceptable fit quality.

Do not select the fit with the smallest error bar automatically.

## Acceptance criterion

Synthetic data with strong small-\(L\) corrections visibly converge toward the known asymptotic exponent when \(L_{\min}\) increases.

---

# 11. Phase H — Statistical diagnostics

Implement common fitting results containing:

- parameter values;
- standard errors;
- covariance matrix;
- correlation matrix;
- \(\chi^2\);
- degrees of freedom;
- reduced \(\chi^2\);
- p-value;
- AIC;
- AICc where meaningful;
- BIC;
- residuals.

## Identifiability guard

Warn when:

- covariance is ill-conditioned;
- absolute parameter correlations are extremely high;
- optimizer solutions depend strongly on initialization;
- parameters are practically unconstrained.

Add optional profile scans.

## Acceptance criterion

The software refuses to present a strongly degenerate multi-parameter fit as a high-precision result.

---

# 12. Phase I — Bootstrap and correlated uncertainty

Implement:

- parametric or nonparametric bootstrap where appropriate;
- raw-sample bootstrap;
- jackknife where useful;
- uncertainty propagation from \(t_c\) to exponent fits.

If covariance matrices for multiple data points/observables are available, allow correlated fitting.

Never silently treat known correlated data as independent.

---

# 13. Phase J — Logarithmic FSS

Implement a general log-correction module.

Support at least:

\[
O(L)=aL^y(\ln L+c)^{\hat y}
\]

and selected correction extensions.

Fit in stages:

1. fixed \(y\), free \(\hat y\);
2. fixed \(\hat y\), free \(y\);
3. optionally free both only if identifiable;
4. test sensitivity to \(c\) or \(L_0\).

Generate profile/correlation diagnostics.

## Key scientific test

Construct synthetic data showing that over a restricted size range:

\[
L^{-\omega}
\]

can mimic a logarithmic correction.

The framework must be capable of returning:

> unresolved between power correction and logarithmic correction.

---

# 14. Phase K — BKT module

Implement essential-singularity scaling.

General form:

\[
\xi\sim\exp(bt^{-\sigma}).
\]

Support a finite-size variable such as

\[
x=t[\ln(L/L_0)]^{1/\sigma}.
\]

Allow \(\sigma=1/2\) as a constrained BKT hypothesis.

Do not enforce it by default.

Support multiplicative logarithms in observables where justified.

## Validation

Generate synthetic BKT-like data and demonstrate that naive ordinary power-law FSS gives unstable or misleading effective exponents.

---

# 15. Phase L — Crossover module

Implement effective-exponent diagnostics and two-field scaling support.

Conceptual form:

\[
O(t,g,L)=L^{y_O}F(tL^{y_t},gL^{y_g}).
\]

At minimum provide:

- local/pairwise effective exponents;
- effective-exponent flow versus \(L\);
- crossover variable \(gL^{y_g}\);
- crossover collapse;
- approximate crossover length.

## Validation

Synthetic data should transition between two apparent scaling regimes.

A naive one-regime fit should show drift.

The crossover analysis should detect that drift rather than declaring a continuously varying exponent.

---

# 16. Phase M — Geometric FSS

Support critical geometric observables such as:

\[
C_1\sim L^{d_F},
\]

cluster radius,

shortest-path length,

backbone mass,

hulls,

loops,

or user-defined geometric quantities.

Do not assume standard thermodynamic exponents apply.

Support direct tests of fractal dimensions and relations among geometric exponents.

---

# 17. Phase N — Distribution FSS

Implement:

\[
P(X,L)=L^{-y_X}\widetilde P(X/L^{y_X})
\]

and

\[
n(s,L)=s^{-\tau}\widetilde n(s/L^{d_F}).
\]

Include:

- normalization diagnostics;
- cutoff-scale estimation;
- moment scaling;
- collapse;
- tail diagnostics;
- hyperscaling checks.

Important:

> Do not estimate \(\tau\) solely from a manually selected straight region on a log-log plot.

Use statistically defensible distribution-fitting methods where possible.

---

# 18. Phase O — Model/ansatz comparison

Create a comparison framework for physically motivated candidate models.

Typical comparisons:

- pure power vs correction-to-scaling;
- power correction vs logarithmic correction;
- ordinary continuous vs BKT;
- asymptotic exponent vs crossover scenario.

Use:

- fit quality;
- residual structure;
- AIC/AICc/BIC where appropriate;
- parameter stability;
- identifiability;
- theoretical consistency.

The system must permit:

> current sizes do not discriminate.

---

# 19. Phase P — Data collapse engine

Implement collapse as a final validation layer.

Support:

## Ordinary

\[
x=tL^{y_t}.
\]

## Log corrected

\[
x=tL^{y_t}(\ln L)^{\hat y_t}.
\]

## BKT

\[
x=t[\ln(L/L_0)]^{1/\sigma}.
\]

## Crossover

Two-variable diagnostics based on:

\[
(tL^{y_t},gL^{y_g}).
\]

Quantify collapse quality.

Never make visual overlap the primary fit criterion unless the method itself is explicitly defined and numerically optimized.

---

# 20. Phase Q — Experimental-design advisor

This is a required research feature.

When competing hypotheses remain unresolved, recommend the simulations most useful for discrimination.

Possible actions:

- increase \(L_{\max}\);
- add intermediate sizes;
- simulate closer to \(t_c\);
- improve statistics at crossings;
- measure a cleaner dimensionless observable;
- measure derivative/covariance observables;
- add raw samples;
- add interaction-range/coupling values;
- choose sizes where competing fits diverge maximally.

The output must explain why each new simulation helps.

---

# 21. Phase R — Synthetic validation suite

Create controlled datasets with known answers.

Required cases:

1. ordinary transition with irrelevant correction;
2. pure power law with small-size bias;
3. multiplicative logarithm;
4. log-vs-power ambiguity;
5. BKT-like essential scaling;
6. crossover between two regimes;
7. geometric fractal scaling;
8. distribution FSS.

All generators should use fixed random seeds.

All tests must verify recovery within expected uncertainty, not exact floating-point equality.

---

# 22. Phase S — Real benchmark cases

After synthetic tests pass, add compact benchmark data from known models or existing local research outputs.

Suggested coverage:

- standard continuous transition;
- percolation dimensionless crossing;
- geometric observable;
- log-corrected example;
- crossover example.

Do not overwrite original research data.

Copy only small representative datasets into a regression-test directory if permitted.

Compare new results with previous analyses.

Investigate discrepancies instead of forcing agreement.

---

# 23. Phase T — Standard output format

A full analysis should produce:

```text
analysis/
├── tables/
├── figures/
├── fits/
├── diagnostics/
└── FSS_REPORT.md
```

`FSS_REPORT.md` should contain:

1. data description;
2. observable classification;
3. candidate scaling hypotheses;
4. critical-point analysis;
5. exponent/scaling-dimension fits;
6. corrections;
7. \(L_{\min}\) stability;
8. window stability;
9. parameter correlations;
10. statistical uncertainty;
11. systematic uncertainty;
12. model comparison;
13. collapse diagnostics;
14. cross-observable consistency;
15. unresolved issues;
16. proposed new simulations;
17. final conclusion.

Conclusion labels:

- **Robust**
- **Suggestive**
- **Unresolved**

---

# 24. Phase U — CLI and API

First expose reliable low-level scripts.

Examples:

```bash
python scripts/inspect_data.py data.csv
```

```bash
python scripts/critical_power_fit.py \
  data.csv \
  --size L \
  --observable chi \
  --error chi_err
```

```bash
python scripts/stability_scan.py ...
```

Only after the low-level routines are stable, add a higher-level interface such as:

```bash
fss analyze config.yaml
```

Do not prematurely build a large framework.

---

# 25. Phase V — Scientific adversarial review

After implementation, act as a skeptical referee.

Create `SELF_REVIEW.md`.

Check specifically for:

1. visual data-collapse bias;
2. underestimated uncertainty;
3. missing finite-size corrections;
4. overparameterized models;
5. parameter degeneracy;
6. inappropriate BKT fits;
7. unjustified logarithmic terms;
8. crossover mistaken for new universality;
9. misuse of hyperscaling;
10. missing analytic backgrounds;
11. ignored correlations;
12. failure to propagate \(t_c\) uncertainty;
13. overinterpretation of effective exponents;
14. weak \(L_{\min}\) testing;
15. model selection based only on information criteria.

Classify issues by severity.

Fix all high-severity issues that can reasonably be fixed.

---

# 26. Execution strategy for the local Agent

Do not attempt the entire project in one giant implementation pass.

Use the following milestones.

## Milestone 1 — Scientific design

Complete:

- environment audit;
- five paper notes;
- theory references;
- final `SKILL.md`;
- data schema.

No advanced numerical implementation yet.

### Exit criterion

The architecture is scientifically coherent and the skill instructions do not overfit to one model class.

---

## Milestone 2 — Core ordinary FSS

Implement:

- I/O;
- power-law fits;
- corrections;
- crossings;
- derivative scaling;
- common fit-result structure;
- \(L_{\min}\) scan;
- window scan;
- residual diagnostics.

### Exit criterion

All ordinary synthetic tests pass.

---

## Milestone 3 — Statistical robustness

Implement:

- bootstrap;
- correlated fits;
- identifiability warnings;
- profile scans;
- model comparison.

### Exit criterion

The software correctly flags deliberately overparameterized synthetic examples.

---

## Milestone 4 — Nonstandard scaling

Implement:

- logarithmic FSS;
- BKT;
- crossover.

### Exit criterion

The system distinguishes clearly distinguishable synthetic cases and reports ambiguity in deliberately ambiguous cases.

---

## Milestone 5 — Geometry and distributions

Implement:

- geometric FSS;
- distribution scaling;
- moment/hyperscaling checks.

### Exit criterion

Known synthetic \(d_F\) and \(\tau\) are recovered within uncertainty.

---

## Milestone 6 — Research-agent integration

Implement:

- standard report;
- experimental-design advisor;
- examples;
- high-level CLI if justified.

### Exit criterion

A local Agent can analyze a new dataset using only the skill instructions plus the data.

---

## Milestone 7 — Real-world validation and self-review

Run representative real datasets.

Create:

- `IMPLEMENTATION_SUMMARY.md`;
- `SELF_REVIEW.md`.

Fix all high-severity scientific problems.

---

# 27. Definition of Done

The project is not complete until:

- `SKILL.md` is concise and usable;
- detailed theory is in references;
- ordinary continuous FSS works;
- corrections and backgrounds work;
- \(L_{\min}\) and fitting-window scans work;
- residual diagnostics work;
- uncertainty is quantified;
- identifiability warnings work;
- log FSS works;
- BKT is treated separately;
- crossover is supported;
- geometric FSS works;
- distribution FSS works;
- model comparison can return "unresolved";
- synthetic tests pass;
- at least several real or benchmark datasets have been tried;
- the skill can recommend additional simulations;
- `pytest` passes;
- `SELF_REVIEW.md` has been completed.

---

# 28. Prompt for the local Agent

Use the following as the project-launch instruction:

> You are building the general-purpose finite-size-scaling Agent Skill described in `RESEARCH_PLAN.md`. Work autonomously and scientifically. Begin by auditing the repository and reading the five primary methodological papers if local PDFs are available. Do not jump directly to data-collapse code. Implement the project milestone by milestone. At the end of each milestone, run tests, update `IMPLEMENTATION_SUMMARY.md`, and commit only logically coherent changes if version control is available. Never fabricate paper-specific methodology. If a scientific ambiguity cannot be resolved from available data or references, document it instead of inventing a rule. Proceed until the current milestone's exit criteria are satisfied before moving on.

Recommended first run:

> Execute Milestones 1 and 2 only. Do not start BKT/log/crossover modules until the ordinary-FSS core, stability engine, and synthetic tests pass.

Recommended second run:

> Continue with Milestones 3 and 4. Focus on parameter identifiability, uncertainty propagation, logarithmic FSS, BKT, and crossover. Add adversarial synthetic tests designed to fool naive power-law fitting.

Recommended third run:

> Continue with Milestones 5–7. Add geometric/distribution FSS, research-report generation, experimental-design advice, representative real-data regression tests, and the full scientific self-review.

---

# 29. V2 extensions

Do not implement these until V1 is reliable:

- Bayesian hierarchical FSS;
- Gaussian-process scaling functions;
- automated symbolic generation of RG expansions;
- joint multi-observable global fits;
- multi-level covariance models;
- corrections for autocorrelated Monte Carlo time series;
- automated reweighting interfaces;
- histogram reweighting;
- neural/data-driven scaling-function inference;
- active-learning simulation design;
- automated universality-class comparison across models;
- dynamic FSS;
- anisotropic FSS;
- quantum critical finite-size / finite-temperature scaling.

These are future extensions, not V1 requirements.
