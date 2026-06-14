# RCuredOrdinal

Two-stage estimation for an **ordinal outcome** regressed on a **censored covariate** in the presence of a **cured fraction**.

Some individuals never experience the event (they are *cured*, so the covariate is undefined for them); among those who will, the event time is observed only if it occurred before a cross-sectional measurement time. Complete-case analysis discards everyone without an observed event; a naive analysis pools the cured with the not-yet-observed uncured and is biased. This code fits a model that separates the two and estimates covariate effects for every subgroup — including the uncured whose event has not yet been observed.

Reference application: the Adult Perthes' study (International Perthes' Disease Study Group), relating patient-reported physical functioning (SF-36 Role Physical) to age-at-total-hip-replacement.

> **Status: research code, not an installable package.** Despite the `R/` layout and roxygen comments, there is no `DESCRIPTION` or `NAMESPACE`, so `install.packages()` / `devtools::install_github()` will not work. Load the code with `source()` as shown below. See [Known issues](#known-issues).

---

## Method

Estimation proceeds in two stages.

**Stage 1 — mixture cure model**, fitted by EM (`R/cure_model.R`):

- *incidence*: logistic regression for the probability of being uncured;
- *latency*: Weibull proportional-hazards regression for the event time among the uncured.

The E-step produces, for each censored individual, a posterior probability of being uncured. These become the weights passed to Stage 2.

**Stage 2 — ordinal pseudo-likelihood** (`R/ordinal_pseudo_likelihood.R`), maximized by `optim(method = "L-BFGS-B")` and seeded by naive and complete-case fits (`R/initialization.R`). Individuals with an observed event contribute directly; censored individuals contribute a weighted mixture of a cured component and an uncured-not-yet-observed component.

**Variance** (`R/variance.estimation.R`): a stacked estimating-equations sandwich that propagates first-stage uncertainty into the second-stage standard errors. The ordinary inverse-Hessian understates it.

Two ordinal link families are supported: `"PO"` (proportional odds / cumulative logit) and `"ACAT"` (adjacent categories).

---

## Requirements

R, plus:

```r
install.packages(c(
  "MASS", "VGAM", "survival", "eha", "smcure", "dplyr",
  "numDeriv", "matrixStats", "forcats", "emmeans", "ggplot2"
))
```

The plotting scripts additionally use `scales`, `tidyr`, `gridExtra`, and `resample`.

---

## Setup

Clone the repo and open `RCuredOrdinal.Rproj`, or set the working directory to the repository root. **All paths in the scripts are relative to the repository root**, so run from there.

```r
source("analysis/_setup.R")   # loads packages, sets the seed, sources R/
```

---

## Quick start — simulation

```r
source("analysis/main_simulation.R")
```

This generates data from known parameters and fits the two-stage estimator alongside the naive and complete-case comparators. Edit the top of the script to change the truth, `n`, `replications`, or the outcome link.

To fit a single dataset directly:

```r
source("analysis/_setup.R")

dat <- gen_demo_data(n = 2000, k = 3, par = true.params, outcome.model = "PO")

fit <- ordcure(
  cureform  = delta ~ Z1 + Z2,        # incidence:  uncured yes/no
  survform  = Tau   ~ Z1 + Z2,        # latency:    event time | uncured
  formula.a = V ~ R + Z1 + Z2,        # cured
  formula.b = V ~ R,                  # uncured, event not yet observed
  formula.c = V ~ R + Tau + R:Tau,    # uncured, event observed
  formula.e = V ~ Z1 + Z2,            # effects shared across uncured subgroups
  Tau   = "Tau",                      # observed time = min(event, censoring)
  R     = "R",                        # measurement time
  delta = "delta",                    # 1 if the event was observed
  data  = dat,
  outcome.model = "PO",
  var = TRUE,                         # sandwich variance
  verbose = TRUE
)

fit$par.list     # estimates by block
fit$variance     # sandwich variance (when var = TRUE)
fit$first.stage  # cure-model fit and weights
fit$naive
fit$cc
```

The response must be an **ordered factor**. Covariates entering the outcome blocks must be numeric or 0/1 dummies — `.build_design()` in `R/outcome_probabilities.R` forms interaction columns by multiplication and does not apply factor contrasts.

---

## Parameter blocks

The code uses `a` / `e` / `b` / `c` where the thesis uses Greek letters:

| Code | Symbol | Subgroup | Formula argument |
|------|--------|----------|------------------|
| `a` | α | cured | `formula.a` |
| `b` | β | uncured, event **not yet** observed | `formula.b` |
| `c` | γ | uncured, event **observed** | `formula.c` |
| `e` | η | effect shared across both uncured subgroups | `formula.e` |

A covariate belongs to `formula.e` **or** to `formula.b`/`formula.c`, never both — otherwise the effects are not separately identified. Note also that the event time may appear in `formula.c` but not in `formula.b`: for the not-yet-observed group it has not happened, and the observed-time column equals the measurement time there.

Under `"PO"` with `k > 2`, the category thresholds are optimized as (first threshold, successive gaps) so monotonicity becomes a non-negativity bound; `ordcure()` converts them back to cumulative thresholds before returning.

---

## Perthes' data analysis

The dataset is **not** in this repository. Create a `config.R` at the repository root (it is gitignored):

```r
PERTHES_DATA_PATH <- "/full/path/to/data.txt"
```

Then:

```r
source("R/run_Perthes_Data.R")
```

The script expects a tab-separated file with columns `age.THR`, `age`, `THR`, `sex`, `unilateral`, `dysplasia`, `diag_age6_7`, `diag_age8_11`, `diag_age11`, `chPtrt.surgery`, `chPtrt.activity.restrict`, and `sf36_RP`. It caps `age.THR` at `age`, jitters the ties backwards within the same year under `set.seed(17)`, and collapses the three middle SF-36 categories to give K = 3.

---

## Repository layout

```
R/
  ordcure.R                    main entry point; validation and orchestration
  cure_model.R                 stage 1: Weibull mixture cure model by EM
  initialization.R             naive and complete-case fits; stage-2 starting values
  ordinal_pseudo_likelihood.R  stage 2: negative pseudo log-likelihood
  outcome_probabilities.R      category probabilities under PO and ACAT
  variance.estimation.R        stacked estimating-equations sandwich
  gen_data.R                   synthetic data generation
  run_simulation.R             replication loop
  simulation_utils.R           summaries, coverage indicators
  run_Perthes_Data.R           Perthes' analysis driver
  EstPerthesData.R             legacy — see Known issues
analysis/
  _setup.R                     packages, seed, sources R/
  main_simulation.R            simulation study driver
plots/
  Plots.R                      legacy — see Known issues
  EstPerthes-plots.R           legacy — see Known issues
  *.png                        outcome-probability figures
```

---

## Known issues

- **Not an installable package.** No `DESCRIPTION`, `NAMESPACE`, `man/`, or tests. The roxygen blocks and `@export` tags are not currently used to generate anything.
- **Stale scripts.** `R/EstPerthesData.R`, `plots/Plots.R`, and `plots/EstPerthes-plots.R` `source()` files that no longer exist in the repository (`gen.probs.R`, `psudeo.likelihood.R`, `gen.data.R`, `run.sim.R`, and `analysis.R`, the last of which is gitignored). They will not run as-is and predate the current `R/` layout. The `.png` files in `plots/` are their output.
- **`build_summary_tables()`** is commented out in `analysis/main_simulation.R` with a note that it needs fixing.
- **The `R` argument to `ordcure()`** is validated but not used in the function body; the measurement time enters through the outcome formulas instead.
- **First-stage Hessian.** `variance.estimation.R` computes both a numeric Hessian (`numDeriv::hessian`, used) and an analytic one (`M2`, computed but not wired in). See the note at the top of that file to switch.

---

## Citation

Ziv, S. (2026). *Ordinal Outcome Regression with Censored Covariates in the Presence of a Cured Fraction.* MSc thesis, Department of Statistics, University of Haifa. Advisor: Bella Vakulenko-Lagun. Joint work with M. B. Millis and H. K. W. Kim (International Perthes' Study Group).

With thanks to the International Perthes' Study Group and the study participants. Supported by the Israel Science Foundation (grant 1219/23).
