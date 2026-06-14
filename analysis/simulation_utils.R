############################################################
## Simulation utilities
############################################################

# Validate the inputs to run_simulation(). Internal to the analysis scripts.
validate_sim_inputs <- function(seed, n, k, replications, par,
                                outcome.model, survform, cureform,
                                formula.a, formula.b, formula.c, formula.e,
                                Tau, R, delta, var, alpha) {
  
  ## simulation control
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1 || seed < 0))
    stop("seed must be a non-negative numeric scalar or NULL")
  
  if (!is.numeric(replications) || replications <= 0 ||
      replications != as.integer(replications))
    stop("replications must be a positive integer")
  
  if (!is.numeric(n) || n <= 0 || n != as.integer(n))
    stop("n must be a positive integer")
  
  if (n < 30)
    warning("n is very small; ordinal cure models may be unstable")
  
  if (!is.numeric(k) || length(k) != 1 || k < 2 || k != as.integer(k))
    stop("k must be an integer >= 2")
  
  ## model choice
  if (!outcome.model %in% c("PO", "ACAT"))
    stop("outcome.model must be either 'PO' or 'ACAT'")
  
  ## true parameters
  if (!is.list(par))
    stop("par must be a list")
  
  if (any(!is.finite(unlist(par))))
    stop("par contains non-finite values")
  
  if (!all(vapply(par, is.numeric, logical(1))))
    stop("all elements of par must be numeric vectors")
  
  required_blocks <- c("a", "b", "c", "d", "e", "Tau", "Tau.shape", "Tau.scale")
  missing_blocks  <- setdiff(required_blocks, names(par))
  if (length(missing_blocks) > 0)
    stop("par is missing components: ", paste(missing_blocks, collapse = ", "))
  
  ## formulas
  formulas <- list(survform, cureform, formula.a, formula.b, formula.c, formula.e)
  if (!all(vapply(formulas, inherits, logical(1), what = "formula")))
    stop("All model formulas must be formula objects")
  
  ## estimation options
  if (!is.logical(var) || length(var) != 1)
    stop("var must be TRUE or FALSE")
  
  if (var && (!is.numeric(alpha) || alpha <= 0 || alpha >= 1))
    stop("alpha must be in (0,1) when var = TRUE")
  
  ## variable names
  if (!all(vapply(list(Tau, R, delta),
                  function(x) is.character(x) && length(x) == 1, logical(1))))
    stop("Tau, R, delta must be single character strings")
  
  invisible(TRUE)
}

# Coverage indicator: TRUE when the true value lies in the Wald CI.
CI_indicator <- function(est, se, alpha, true)
  abs(est - true) <= qnorm(1 - alpha / 2) * se

# Per-parameter summary across replications: average estimate, bias,
# empirical SD, and (when supplied) mean SE and empirical coverage.
summary_sim <- function(est, true, se = NULL, ci = NULL,
                        se.inv = NULL, ci.inv = NULL) {
  
  avg  <- colMeans(est)
  bias <- true - avg
  SD   <- colSds(est)
  
  out <- data.frame(
    true = true, avg = avg, bias = bias, SD = SD,
    row.names = names(true),
    check.names = FALSE
  )
  
  if (!is.null(se)) {
    out[["SE-sandwich"]] <- colMeans(se)
    out[["CP-sandwich"]] <- colMeans(ci)
  }
  
  if (!is.null(se.inv)) {
    out[["SE-inv"]] <- colMeans(se.inv)
    out[["CP-inv"]] <- colMeans(ci.inv)
  }
  
  out
}

# --------------------------------------------------
# Combined summary table with relative efficiency
# --------------------------------------------------
# Aligns the proposed estimator's per-parameter SEs against the naive
# and complete-case (CC) estimators and reports relative efficiency
#   RE = Var(competitor) / Var(proposed),
# so RE > 1 favors the proposed estimator. The naive estimator only
# covers the alpha (V.alpha.*) block; the CC estimator covers the
# gamma/eta (V.gamma.* / V.eta.*) blocks; first-stage and beta
# parameters have no competitor and are left as NA.
#
# Requires a simulation run with var = TRUE.
build_summary_tables <- function(sim) {
  
  est   <- sim$summary.table.est
  naive <- sim$summary.table.naive
  cc    <- sim$summary.table.cc
  
  if (!"SE-sandwich" %in% colnames(est))
    stop("build_summary_tables() needs a run with var = TRUE")
  
  ## bare parameter key (strip the leading block prefix: a. b. c. e. d. Tau.)
  key    <- sub("^[A-Za-z]+\\.", "", rownames(est))
  se_est <- est[["SE-sandwich"]]
  
  ## look up a competitor table by key, aligned to est's rows
  lookup <- function(tab, eligible) {
    se <- cp <- rep(NA_real_, nrow(est))
    hit <- eligible & key %in% rownames(tab)
    se[hit] <- tab[key[hit], "SE-sandwich"]
    cp[hit] <- tab[key[hit], "CP-sandwich"]
    list(se = se, cp = cp)
  }
  
  nv <- lookup(naive, grepl("^V\\.alpha\\.", key))
  cv <- lookup(cc,    grepl("^V\\.(gamma|eta)\\.", key))
  
  base_cols <- intersect(
    c("true", "avg", "bias", "SD", "SE-sandwich", "SE-inv", "CP-sandwich", "CP-inv"),
    colnames(est)
  )
  out <- est[, base_cols, drop = FALSE]
  
  out[["SE-naive"]] <- nv$se
  out[["CP-naive"]] <- nv$cp
  out[["RE-naive"]] <- (nv$se / se_est)^2
  out[["SE-cc"]]    <- cv$se
  out[["CP-cc"]]    <- cv$cp
  out[["RE-cc"]]    <- (cv$se / se_est)^2
  
  ## split first-stage from second-stage parameters for readability
  first_stage <- grepl("^(d|Tau)\\.", rownames(out))
  
  list(
    first_stage  = out[first_stage, , drop = FALSE],
    second_stage = out[!first_stage, , drop = FALSE]
  )
}