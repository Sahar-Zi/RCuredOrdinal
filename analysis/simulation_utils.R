############################################################
## Simulation utilities
## Ordinal Cure Model with Censored Covariates
############################################################

validate_sim_inputs <- function(seed, n, k, replication, par,
                                outcome.model, survform, cureform,
                                formula.a, formula.b, formula.c, formula.e,
                                Tau, R, delta, var, alpha) {
  
  ## Simulation control
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1 || seed < 0))
    stop("seed must be a non-negative numeric scalar or NULL")
  
  if (!is.numeric(replication) || replication <= 0 ||
      replication != as.integer(replication))
    stop("replication must be a positive integer")
  
  if (!is.numeric(n) || n <= 0 || n != as.integer(n))
    stop("n must be a positive integer")
  
  if (n < 30)
    warning("n is very small; ordinal cure models may be unstable")
  
  if (!is.numeric(k) || length(k) != 1 ||
      k < 2 || k != as.integer(k))
    stop("k must be an integer >= 2")
  
  ## Model choice
  if (!outcome.model %in% c("PO", "ACAT"))
    stop("outcome model must be either 'PO' or 'ACAT'")
  
  ## True parameters
  if (!is.list(par))
    stop("par must be a list")
  
  if (any(!is.finite(unlist(par))))
    stop("par contains non-finite values")
  
  if (!all(sapply(par, is.numeric)))
    stop("all elements of par must be numeric vectors")
  
  required_blocks <- c("a", "b", "c", "d", "e", "Tau", "shape", "scale")
  
  missing_blocks <- setdiff(required_blocks, names(par))
  if (length(missing_blocks) > 0)
    stop("par is missing components: ",
         paste(missing_blocks, collapse = ", "))
  
  ## Formulas
  formulas <- list(survform, cureform, formula.a,
                   formula.b, formula.c, formula.e)
  
  if (!all(sapply(formulas, inherits, "formula")))
    stop("All model formulas must be formula objects")
  
  ## Estimation options
  if (!is.logical(var) || length(var) != 1)
    stop("var must be TRUE or FALSE")
  
  if (var && (!is.numeric(alpha) || alpha <= 0 || alpha >= 1))
    stop("alpha must be in (0,1) when var = TRUE")
  
  ## Variable names
  if (!all(sapply(list(Tau, R, delta),
                  function(x) is.character(x) && length(x) == 1)))
    stop("Tau, R, delta must be single character strings")
  
  invisible(TRUE)
}


CI_indicator <- function(est, se, alpha, true)
  abs(est - true) <= qnorm(1 - alpha / 2) * se

summary_sim <- function(
    est,
    true,
    se = NULL,
    ci = NULL,
    se.inv = NULL,
    ci.inv = NULL
) {
  
  avg  <- colMeans(est)
  bias <- true - avg
  SD   <- colSds(est)
  
  out <- data.frame(
    true = true,
    avg = avg,
    bias = bias,
    SD = SD
  )
  
  if (!is.null(se)) {
    out$`SE-sandwich` <- colMeans(se)
    out$`CP-sandwich` <- colMeans(ci)
  }
  
  if (!is.null(se.inv)) {
    out$`SE-inv` <- colMeans(se.inv)
    out$`CP-inv` <- colMeans(ci.inv)
  }
  
  out
}

# --------------------------------------------------
# Build aligned summary tables for simulation output
# --------------------------------------------------

build_summary_tables <- function(sim) {
  
  est   <- sim$summary.table.est
  naive <- sim$summary.table.naive
  cc    <- sim$summary.table.cc
  
  # ---- add parameter column ----
  est$param   <- rownames(est)
  naive$param <- rownames(naive)
  cc$param    <- rownames(cc)
  
  # ---- helper: map eta -> gamma ----
  map_eta_to_gamma <- function(p) {
    sub("^e\\.V\\.eta\\.", "V.gamma.", p)
  }
  
  cc_lookup <- cc
  rownames(cc_lookup) <- cc_lookup$param
  
  # ---- initialize CC columns in est ----
  est$RE     <- NA
  
  for (i in seq_len(nrow(est))) {
    
    p_est <- est$param[i]
    
    if (grepl("^e\\.V\\.eta\\.", p_est)) {
      
      p_cc <- map_eta_to_gamma(p_est)
      
      if (p_cc %in% rownames(cc_lookup)) {
        
        est[p, "avg"] <- cc[p_cc, "avg"]
        est[p, "SE"]  <- cc[p_cc, "SE"]
        est[p, "CP"]  <- cc[p_cc, "CP"]
        
        # Relative efficiency: Var(CC) / Var(EST)
        est[p, "RE"] <- (cc[p_cc, "SE"]^2) /
          (est[p, "SE.sandwich"]^2)
      }
    }
  }
  
  # ---- naive alignment (row names already match) ----
  naive_out <- naive
  colnames(naive_out) <- paste0(colnames(naive_out), "_naive")
  
  # ---- split first-stage vs second-stage ----
  first_stage_idx <- grepl("^(d|Tau|shape|scale)", rownames(est))
  
  est_first  <- est[first_stage_idx, ]
  est_second <- est[!first_stage_idx, ]
  
  # ---- final column order ----
  est_cols <- c(
    "true", "avg", "bias", "SD",
    "SE-sandwich", "SE-inv",
    "CP-sandwich", "CP-inv", "RE"
  )
  
  est_first  <- est_first[, est_cols]
  est_second <- est_second[, est_cols]
  
  return(list(
    est_first_stage  = est_first,
    est_second_stage = est_second,
    naive            = naive_out
  ))
}
