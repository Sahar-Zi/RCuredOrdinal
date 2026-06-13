# ============================================================
# ordcure(): two-stage estimation for an ordinal outcome with a
# cured fraction and a censored covariate.
#
# Stage 1 (cure_model.R)        : Weibull mixture cure model by EM.
# Stage 2 (pseudo_likelihood.R) : ordinal pseudo-likelihood, seeded
#                                 by naive/CC fits (initialization.R).
# Variance (variance.R)         : optional stacked-equations sandwich.
# ============================================================

# ---- input validation -------------------------------------
#' @noRd
validate_ordcure_inputs <- function(
    data, survform, cureform,
    formula.a, formula.e, formula.b, formula.c,
    Tau, R, delta, outcome.model
) {
  
  ## 1. Outcome model
  if (!outcome.model %in% c("PO", "ACAT"))
    stop("outcome.model must be one of: 'PO', 'ACAT'")
  
  ## 2. Required variables present in data
  all_forms <- list(survform, cureform, formula.a, formula.e, formula.b, formula.c)
  form_vars <- unique(unlist(lapply(all_forms, all.vars)))
  required_vars <- unique(c(form_vars, Tau, R, delta))
  
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0)
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
  
  ## 3. Response consistency and type
  resp <- all.vars(formula.a)[1]
  
  same_resp <- vapply(
    list(formula.e, formula.b, formula.c),
    function(f) all.vars(f)[1] == resp,
    logical(1)
  )
  if (!all(same_resp))
    stop("All outcome formulas (a, e, b, c) must share the same response variable")
  
  if (!is.ordered(data[[resp]]))
    stop("Outcome variable must be an ordered factor")
  
  ## 4. Forbidden variable placement
  outcome_vars <- unique(unlist(lapply(
    list(formula.a, formula.e, formula.b, formula.c), all.vars
  )))
  
  if (delta %in% outcome_vars)
    stop("delta cannot appear in outcome model formulas")
  
  if (Tau %in% all.vars(formula.a))
    stop("Tau cannot appear in formula.a")
  
  ## 5. Basic semantic checks
  if (!is.numeric(data[[Tau]]) || any(data[[Tau]] < 0, na.rm = TRUE))
    stop("Tau must be a non-negative numeric variable")
  
  if (!is.logical(data[[delta]]) &&
      !all(na.omit(unique(data[[delta]])) %in% c(0, 1)))
    stop("delta must be logical or binary (0/1)")
  
  invisible(TRUE)
}

# ---- complete-case formula --------------------------------
# Combines the RHS of formula.e and formula.c into the single
# complete-case formula used to fit the delta = 1 rows.
#' @noRd
build_cc_formula <- function(formula.e, formula.c, move_interactions_last = TRUE) {
  
  if (!inherits(formula.e, "formula") || !inherits(formula.c, "formula"))
    stop("formula.e and formula.c must be formula objects")
  
  lhs.e <- all.vars(formula.e)[1]
  lhs.c <- all.vars(formula.c)[1]
  
  if (lhs.e != lhs.c)
    stop("formula.e and formula.c must have the same response")
  
  terms.e <- attr(terms(formula.e), "term.labels")
  terms.c <- attr(terms(formula.c), "term.labels")
  rhs <- unique(c(terms.e, terms.c))
  
  ## keep interactions after main effects for stable column ordering
  if (move_interactions_last) {
    is_inter <- grepl("[:*]|I\\(", rhs)
    rhs <- c(rhs[!is_inter], rhs[is_inter])
  }
  
  if (length(rhs) == 0) {
    as.formula(paste(lhs.e, "~ 1"))
  } else {
    as.formula(paste(lhs.e, "~", paste(rhs, collapse = " + ")))
  }
}

#' Fit an ordinal cure model with a censored covariate
#'
#' Fits an ordinal regression model with a cured fraction and a
#' parametric (Weibull) survival component, allowing different covariate
#' effects before and after the event time. Estimation is two-stage: a
#' mixture cure model by EM, followed by an ordinal pseudo-likelihood.
#'
#' @param survform Latency (survival) formula, e.g. `Tau ~ Z1 + Z2`.
#' @param cureform Incidence (cure) formula, e.g. `delta ~ Z1 + Z2`.
#' @param formula.a Naive outcome formula (fitted on the censored rows).
#' @param formula.e Shared-effect outcome formula.
#' @param formula.b Outcome formula for censored ("beta") effects.
#' @param formula.c Outcome formula for observed ("gamma") effects.
#' @param Tau Character name of the observed-time column.
#' @param R Character name of the censored covariate.
#' @param delta Character name of the event-indicator column.
#' @param data A data frame containing all model variables.
#' @param outcome.model `"PO"` (proportional odds) or `"ACAT"`
#'   (adjacent category).
#' @param var Logical; if `TRUE`, also compute the sandwich variance.
#' @param verbose Logical; if `TRUE`, print stage-by-stage progress.
#'
#' @return An object of class `"ordcure"`: a list with `par.list`
#'   (estimated parameters by block), the stage-2 `opt` object, the
#'   `first.stage` fit, the `naive` and `cc` summaries, and, when
#'   `var = TRUE`, a `variance` component.
#'
#' @export
ordcure <- function(
    survform, cureform,
    formula.a, formula.e, formula.b, formula.c,
    Tau, R, delta, data,
    outcome.model = c("PO", "ACAT"),
    var = FALSE,
    verbose = FALSE
) {
  
  outcome.model <- match.arg(outcome.model)
  
  ## helper: stage-by-stage progress, silent unless verbose
  log_step <- function(...) if (verbose) message(...)
  
  log_step("validation")
  validate_ordcure_inputs(
    data, survform, cureform,
    formula.a, formula.e, formula.b, formula.c,
    Tau, R, delta, outcome.model
  )
  
  ## ---- combine CC formula ----
  formula.cc <- build_cc_formula(formula.e, formula.c)
  
  ## ---- stage 1: mixture cure model ----
  log_step("first stage (cure model)")
  stage1 <- fit_cure_weibull_em(
    data, survform, cureform, Tau, delta,
    tol = 1e-8, maxit = 200, verbose = verbose
  )
  
  ## ---- naive + complete-case fits ----
  log_step("naive and complete-case fits")
  nc <- fit_naive_cc_models(formula.a, formula.cc, data, delta, outcome.model)
  
  ## ---- stage-2 starting values ----
  k <- length(levels(model.extract(model.frame(formula.cc, data), "response")))
  
  log_step("build initial parameters")
  init <- build_initial_parameters(
    nc,
    list(e = formula.e, b = formula.b, c = formula.c),
    outcome.model, k
  )
  
  ## ---- stage 2: ordinal pseudo-likelihood ----
  log_step("stage-2 optimization")
  opt <- optim(
    par           = init$par,
    fn            = negloglik_stage2,
    weights       = stage1$weights,
    k             = k,
    lengths       = init$lengths,
    delta         = delta,
    response      = all.vars(formula.a)[1],
    outcome.model = outcome.model,
    data          = data,
    method        = "L-BFGS-B",
    lower         = init$lower
  )
  
  ## ---- assemble parameters by block ----
  ## First-stage par order: cure coefficients, then latency coefficients
  ## plus Weibull shape/scale. Index by design-matrix width so this does
  ## not depend on the cure model having an intercept.
  n_cure     <- ncol(model.matrix(cureform, data))
  stage1_par <- stage1$par
  
  theta <- unpack_stage2_parameters(opt$par, init$lengths)
  
  par.list <- list(
    d   = stage1_par[seq_len(n_cure)],
    Tau = stage1_par[(n_cure + 1):length(stage1_par)],
    a   = theta$a,
    e   = theta$e,
    b   = theta$b,
    c   = theta$c
  )
  
  ## undo the PO threshold reparameterization (gaps -> cumulative)
  if (outcome.model == "PO" && k > 2) {
    par.list$a[seq_len(k - 1)] <- cumsum(par.list$a[seq_len(k - 1)])
    par.list$b[seq_len(k - 1)] <- cumsum(par.list$b[seq_len(k - 1)])
    par.list$c[seq_len(k - 1)] <- cumsum(par.list$c[seq_len(k - 1)])
  }
  
  out <- list(
    par.list    = par.list,
    opt         = opt,
    first.stage = stage1,
    naive       = nc$naive,
    cc          = nc$cc
  )
  
  ## ---- optional sandwich variance ----
  if (var) {
    log_step("variance estimation")
    out$variance <- variance_est(
      est           = par.list,
      outcome.model = outcome.model,
      data          = data,
      survform      = survform,
      cureform      = cureform,
      Tau           = Tau,
      delta         = delta,
      V             = all.vars(formula.a)[1],
      lengths       = init$lengths
    )
  }
  
  class(out) <- "ordcure"
  out
}