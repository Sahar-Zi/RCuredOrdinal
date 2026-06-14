############################################################
## Simulation driver
## Ordinal Cure Model with Censored Covariates
############################################################

#' Run a Monte Carlo simulation of the ordinal cure model
#'
#' For each replication, generates a data set with gen_demo_data() and
#' fits it with ordcure(). Failed replications are recorded and skipped
#' (their seeds are reported); summaries are computed over the successful
#' ones. Per replication the seed is `seed + i - 1`, so runs are
#' reproducible and each replication is independent of failures elsewhere.
#'
#' @param seed Base seed; replication i uses `seed + i - 1`.
#' @param n Sample size.
#' @param k Number of ordinal categories.
#' @param par List of true parameters (see gen_demo_data()).
#' @param outcome.model "PO" or "ACAT".
#' @param survform,cureform,formula.a,formula.e,formula.b,formula.c Model formulas.
#' @param Tau,R,delta Character column names.
#' @param var Logical; compute the sandwich variance (needed for SE/CP/RE).
#' @param replications Number of replications.
#' @param alpha Nominal level for coverage.
#' @param verbose Passed to ordcure().
#'
#' @return A list of raw estimates, optional SE/CI matrices, and the
#'   per-estimator summary tables.
run_simulation <- function(
    seed, n, k, par,
    outcome.model = c("PO", "ACAT"),
    survform, cureform,
    formula.a, formula.e, formula.b, formula.c,
    Tau, R, delta,
    var = FALSE, replications = 1, alpha = 0.05,
    verbose = FALSE
) {
  
  outcome.model <- match.arg(outcome.model)
  
  validate_sim_inputs(
    seed, n, k, replications, par, outcome.model,
    survform, cureform, formula.a, formula.b, formula.c, formula.e,
    Tau, R, delta, var, alpha
  )
  
  ## ---- storage --------------------------------------------
  p <- length(unlist(par))
  
  est.val   <- matrix(NA_real_, replications, p)
  naive.val <- matrix(NA_real_, replications, k + 2)
  cc.val    <- matrix(NA_real_, replications, k + 4)
  
  est.se   <- est.se.i <- matrix(NA_real_, replications, p)
  naive.se <- matrix(NA_real_, replications, k + 2)
  cc.se    <- matrix(NA_real_, replications, k + 4)
  ci.est   <- ci.est.i <- matrix(NA_real_, replications, p)
  ci.naive <- matrix(NA_real_, replications, k + 2)
  ci.cc    <- matrix(NA_real_, replications, k + 4)
  
  seeds   <- integer(replications)
  success <- logical(replications)
  
  ## true values for the naive / CC estimators
  true.a  <- unlist(par$a)
  true.cc <- c(unlist(par$c)[1:(k - 1)], unlist(par$e), unlist(par$c)[k:length(par$c)])
  
  ## name the columns so the result is self-describing (used by plot())
  colnames(est.val) <- colnames(est.se) <- colnames(est.se.i) <-
    colnames(ci.est) <- colnames(ci.est.i) <- names(unlist(par))
  colnames(naive.val) <- colnames(naive.se) <- colnames(ci.naive) <- names(true.a)
  colnames(cc.val)    <- colnames(cc.se)    <- colnames(ci.cc)    <- names(true.cc)
  
  ## ---- replication loop -----------------------------------
  for (i in seq_len(replications)) {
    
    seed_i   <- seed + i - 1
    seeds[i] <- seed_i
    set.seed(seed_i)
    
    dat <- gen_demo_data(n = n, k = k, par = par, outcome.model = outcome.model)
    
    fit <- tryCatch(
      ordcure(
        survform = survform, cureform = cureform,
        formula.a = formula.a, formula.e = formula.e,
        formula.b = formula.b, formula.c = formula.c,
        Tau = Tau, R = R, delta = delta, data = dat,
        outcome.model = outcome.model, var = var, verbose = verbose
      ),
      error = function(e) {
        message("replication ", i, " (seed ", seed_i, ") failed: ",
                conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) next
    
    est.val[i, ]   <- unlist(fit$par.list)
    naive.val[i, ] <- unlist(fit$naive[, "Estimate"])
    cc.val[i, ]    <- unlist(fit$cc[, "Estimate"])
    
    if (var) {
      est.se[i, ]   <- sqrt(diag(fit$variance$stacked.v.est))
      est.se.i[i, ] <- sqrt(diag(-fit$variance$G.tilde.inv))
      naive.se[i, ] <- unlist(fit$naive[, "Std. Error"])
      cc.se[i, ]    <- unlist(fit$cc[, "Std. Error"])
      
      ci.est[i, ]   <- CI_indicator(est.val[i, ], est.se[i, ],   alpha, unlist(par))
      ci.est.i[i, ] <- CI_indicator(est.val[i, ], est.se.i[i, ], alpha, unlist(par))
      ci.naive[i, ] <- CI_indicator(naive.val[i, ], naive.se[i, ], alpha, true.a)
      ci.cc[i, ]    <- CI_indicator(cc.val[i, ], cc.se[i, ],     alpha, true.cc)
    }
    
    success[i] <- TRUE
  }
  
  ## ---- failure accounting ---------------------------------
  n_ok <- sum(success)
  if (n_ok == 0) stop("all replications failed")
  if (n_ok < replications)
    warning(replications - n_ok, " of ", replications, " replications failed")
  
  keep <- success
  sub  <- function(M) M[keep, , drop = FALSE]
  
  ## ---- raw output -----------------------------------------
  out <- list(
    seeds          = seeds,
    success        = success,
    n_replications = replications,
    n_success      = n_ok,
    n              = n,
    outcome.model  = outcome.model,
    true.params    = par,
    est.values       = sub(est.val),
    naive.est.values = sub(naive.val),
    cc.est.values    = sub(cc.val)
  )
  
  if (var) {
    out$est.se.values     <- sub(est.se)
    out$est.se.inv.values <- sub(est.se.i)
    out$CI.est            <- sub(ci.est)
    out$CI.est.inv        <- sub(ci.est.i)
    out$naive.se          <- sub(naive.se)
    out$cc.se             <- sub(cc.se)
    out$CI.naive          <- sub(ci.naive)
    out$CI.cc             <- sub(ci.cc)
  }
  
  ## ---- summaries ------------------------------------------
  out$summary.table.est <- summary_sim(
    est    = sub(est.val), true = unlist(par),
    se     = if (var) sub(est.se)   else NULL,
    ci     = if (var) sub(ci.est)   else NULL,
    se.inv = if (var) sub(est.se.i) else NULL,
    ci.inv = if (var) sub(ci.est.i) else NULL
  )
  
  out$summary.table.naive <- summary_sim(
    est = sub(naive.val), true = true.a,
    se  = if (var) sub(naive.se) else NULL,
    ci  = if (var) sub(ci.naive) else NULL
  )
  
  out$summary.table.cc <- summary_sim(
    est = sub(cc.val), true = true.cc,
    se  = if (var) sub(cc.se) else NULL,
    ci  = if (var) sub(ci.cc) else NULL
  )
  
  class(out) <- "ordcure_sim"
  out
}