# ============================================================
# Stage-2 initialization.
#
# Fits the "naive" ordinal model (delta = 0 rows) and the
# "complete-case" model (delta = 1 rows), then maps their
# coefficients onto the stage-2 (alpha / eta / beta / gamma)
# parameterization to build starting values and lower bounds for
# the pseudo-likelihood optimizer.
#
# Internal (not exported); called by ordcure().
# ============================================================

# Standard errors from a fit's vcov, guarding against non-positive
# variances (which signal a non-PD information matrix). Returns NA in
# that case rather than silently reporting zero.
#' @noRd
.safe_se <- function(fit) {
  v <- diag(vcov(fit))
  if (any(v < 0, na.rm = TRUE)) {
    warning("non-positive variances in fit; standard errors set to NA")
    return(rep(NA_real_, length(v)))
  }
  sqrt(v)
}

#' Fit naive and complete-case ordinal models
#'
#' @param formula.a Naive-model formula (fitted on `delta == 0`).
#' @param formula.cc Complete-case formula (fitted on `delta == 1`).
#' @param data A data frame.
#' @param delta Character name of the event-indicator column.
#' @param outcome.model `"PO"` or `"ACAT"`.
#'
#' @return A list with `naive` and `cc`, each a matrix with columns
#'   `Estimate` and `Std. Error`.
#' @noRd
fit_naive_cc_models <- function(
    formula.a, formula.cc, data, delta, outcome.model
) {

  if (!inherits(formula.a, "formula") || !inherits(formula.cc, "formula"))
    stop("formula.a and formula.cc must be formula objects")

  if (!outcome.model %in% c("PO", "ACAT"))
    stop("outcome.model must be 'PO' or 'ACAT'")

  ## ---- split data by event indicator ----
  data0 <- data[data[[delta]] == 0, , drop = FALSE]   # naive
  data1 <- data[data[[delta]] == 1, , drop = FALSE]   # complete case

  ## ---- outcome info ----
  k <- length(levels(data[[all.vars(formula.cc)[1]]]))

  vars.a  <- attr(terms(formula.a),  "term.labels")
  vars.cc <- attr(terms(formula.cc), "term.labels")

  ## ---- fit models ----
  if (outcome.model == "PO" && k > 2) {

    naive.fit <- MASS::polr(formula.a,  data = data0, Hess = TRUE)
    cc.fit    <- MASS::polr(formula.cc, data = data1, Hess = TRUE)

    ## polr orders parameters as (coefficients, thresholds); the stage-2
    ## convention is (thresholds, coefficients), so reorder accordingly.
    est.naive <- unname(c(naive.fit$zeta, naive.fit$coefficients))
    est.cc    <- unname(c(cc.fit$zeta,   cc.fit$coefficients))

    std       <- .safe_se(naive.fit)
    se.naive  <- unname(c(tail(std, k - 1), head(std, length(std) - k + 1)))

    std       <- .safe_se(cc.fit)
    se.cc     <- unname(c(tail(std, k - 1), head(std, length(std) - k + 1)))

  } else {

    naive.fit <- VGAM::vglm(formula.a,  VGAM::acat(reverse = FALSE, parallel = TRUE), data = data0)
    cc.fit    <- VGAM::vglm(formula.cc, VGAM::acat(reverse = FALSE, parallel = TRUE), data = data1)

    est.naive <- unname(coef(naive.fit))
    est.cc    <- unname(coef(cc.fit))

    se.naive  <- unname(.safe_se(naive.fit))
    se.cc     <- unname(.safe_se(cc.fit))
  }

  names(est.naive) <- names(se.naive) <-
    c(paste0("V.alpha.", seq_len(k - 1)),
      paste0("V.alpha.", vars.a))

  names(est.cc) <- names(se.cc) <-
    c(paste0("V.gamma.", seq_len(k - 1)),
      paste0("V.gamma.", vars.cc))

  summary.n <- cbind(est.naive, se.naive)
  summary.c <- cbind(est.cc,   se.cc)
  colnames(summary.n) <- colnames(summary.c) <- c("Estimate", "Std. Error")

  list(
    naive = summary.n,
    cc    = summary.c
  )
}

# Pull the coefficients matching a formula's RHS terms out of the
# complete-case estimates and relabel them with a target prefix
# (eta / beta / gamma).
#
# NOTE: the `Tau:R` -> `R:Tau` swap normalizes the one interaction
# ordering that terms() and the fitted names disagree on. This is a
# known fragile spot; revisit if more interaction terms are added.
#' @noRd
extract_params <- function(estimates, formula, include_intercept = FALSE, prefix) {

  terms <- attr(terms(formula), "term.labels")
  terms[terms == "Tau:R"] <- "R:Tau"

  pat <- paste0("^V\\.gamma\\.(", paste(terms, collapse = "|"), ")$")
  matched <- estimates[grep(pat, names(estimates))]
  matched <- matched[match(terms, sub("^V\\.gamma\\.", "", names(matched)))]

  if (include_intercept) {
    intercepts <- estimates[grep("^V\\.gamma\\.[0-9]+$", names(estimates))]
    matched <- c(intercepts, matched)
  }

  names(matched) <- sub("^V\\.gamma", paste0("V.", prefix), names(matched))
  matched
}

#' Build stage-2 starting values, block lengths, and lower bounds
#'
#' @param nc Output of `fit_naive_cc_models()`.
#' @param formulas List with elements `e`, `b`, `c`.
#' @param outcome.model `"PO"` or `"ACAT"`.
#' @param k Number of ordinal levels.
#'
#' @return A list with `par` (starting values), `lengths` (block sizes),
#'   and `lower` (per-parameter lower bounds for `optim`).
#' @noRd
build_initial_parameters <- function(nc, formulas, outcome.model, k) {

  a.par <- nc$naive[, "Estimate"]
  e.par <- extract_params(nc$cc[, "Estimate"], formulas$e, FALSE, "eta")
  b.par <- extract_params(nc$cc[, "Estimate"], formulas$b, TRUE,  "beta")
  c.par <- extract_params(nc$cc[, "Estimate"], formulas$c, TRUE,  "gamma")

  lengths <- c(length(a.par), length(e.par), length(b.par), length(c.par))

  ## ---- PO reparameterization: store thresholds as (first, gaps) ----
  if (outcome.model == "PO" && k > 2) {
    a.par[seq_len(k - 1)] <- c(a.par[1], diff(a.par[seq_len(k - 1)]))
    b.par[seq_len(k - 1)] <- c(b.par[1], diff(b.par[seq_len(k - 1)]))
    c.par[seq_len(k - 1)] <- c(c.par[1], diff(c.par[seq_len(k - 1)]))
  }

  ## ---- combine parameters ----
  par <- c(a.par, e.par, b.par, c.par)

  ## ---- lower bounds (gaps between thresholds must be non-negative) ----
  lower.b <- rep(-Inf, length(par))

  if (outcome.model == "PO" && k > 2) {
    idx.a <- seq_len(k - 1)
    idx.b <- lengths[1] + lengths[2] + seq_len(k - 1)
    idx.c <- lengths[1] + lengths[2] + lengths[3] + seq_len(k - 1)

    lower.b[idx.a] <- c(-Inf, rep(0, k - 2))
    lower.b[idx.b] <- c(-Inf, rep(0, k - 2))
    lower.b[idx.c] <- c(-Inf, rep(0, k - 2))
  }

  list(
    par     = par,
    lengths = lengths,
    lower   = lower.b
  )
}
