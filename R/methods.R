# ============================================================
# S3 methods for objects of class "ordcure"
#
# Give ordcure() an lm()-like interface:
#   print, coef, vcov, confint, summary (+print), predict
# ============================================================

#' @export
print.ordcure <- function(x, digits = 4, ...) {
  cat("Ordinal cure model (ordcure)\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }
  cat("\nOutcome model: ", x$outcome.model,
      "   |   n = ", x$nobs, "\n\n", sep = "")
  cat("Coefficients:\n")
  print(round(unlist(x$par.list), digits))
  invisible(x)
}

#' Extract coefficients
#'
#' @param object An "ordcure" fit.
#' @param block Optional block name ("d", "Tau", "a", "e", "b", "c") to
#'   return just that parameter group.
#' @param ... Ignored.
#' @export
coef.ordcure <- function(object, block = NULL, ...) {
  if (!is.null(block)) {
    if (!block %in% names(object$par.list))
      stop("block must be one of: ", paste(names(object$par.list), collapse = ", "))
    return(object$par.list[[block]])
  }
  unlist(object$par.list)
}

#' Sandwich variance-covariance matrix
#'
#' Rows/columns are labelled to match coef(); requires a fit with var = TRUE.
#' @export
vcov.ordcure <- function(object, ...) {
  if (is.null(object$variance))
    stop("no variance available; refit with ordcure(..., var = TRUE)")
  V  <- object$variance$stacked.v.est
  nm <- names(unlist(object$par.list))   # same order as the sandwich blocks
  dimnames(V) <- list(nm, nm)
  V
}

#' Wald confidence intervals
#' @export
confint.ordcure <- function(object, parm, level = 0.95, ...) {
  est <- coef(object)
  se  <- sqrt(diag(vcov(object)))
  z   <- qnorm(1 - (1 - level) / 2)
  ci  <- cbind(est - z * se, est + z * se)
  a   <- (1 - level) / 2
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' Summarize an ordcure fit
#'
#' Builds a coefficient table (estimate, SE, z, two-sided Wald p-value).
#' @export
summary.ordcure <- function(object, ...) {
  
  est <- unlist(object$par.list)
  
  if (!is.null(object$variance)) {
    se <- sqrt(diag(object$variance$stacked.v.est))   # same order as est
    z  <- est / se
    p  <- 2 * pnorm(abs(z), lower.tail = FALSE)        # two-sided Wald
    tab <- cbind(Estimate = est, `Std. Error` = se,
                 `z value` = z, `Pr(>|z|)` = p)
  } else {
    tab <- cbind(Estimate = est)
  }
  rownames(tab) <- names(est)
  
  structure(
    list(call          = object$call,
         outcome.model = object$outcome.model,
         coefficients  = tab,
         logLik        = -object$opt$value,
         nobs          = object$nobs,
         has.var       = !is.null(object$variance)),
    class = "summary.ordcure"
  )
}

#' @export
print.summary.ordcure <- function(x, digits = 4, signif.stars = TRUE, ...) {
  cat("Ordinal cure model (ordcure)\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }
  cat("\nOutcome model: ", x$outcome.model,
      "   |   n = ", x$nobs,
      "   |   logLik = ", round(x$logLik, 2), "\n\n", sep = "")
  
  if (x$has.var) {
    stats::printCoefmat(x$coefficients, digits = digits,
                        signif.stars = signif.stars,
                        P.values = TRUE, has.Pvalue = TRUE)
  } else {
    cat("Coefficients (no variance; refit with var = TRUE for SE/z/p):\n")
    print(round(x$coefficients, digits))
  }
  invisible(x)
}

#' Predicted category probabilities
#'
#' @param object An "ordcure" fit.
#' @param newdata Data frame of covariates. Must contain the event-indicator
#'   column (used to choose the beta/gamma effects within the uncured
#'   component) and any covariates referenced by the model formulas.
#' @param type "probs" returns the n x k probability matrix; "class" returns
#'   the most probable category as a factor.
#' @param component "uncured" returns the observed-data component (beta for
#'   delta = 0, gamma for delta = 1); "cured" returns the baseline (alpha)
#'   component.
#' @param delta Name of the column in newdata that selects the uncured
#'   sub-model (beta where delta = 0, gamma where delta = 1). Defaults to the
#'   column the model was fit with; override it when the prediction regime
#'   differs (e.g. a plotting grid keyed on an "R > Tau" indicator).
#' @param se.fit If TRUE, also return delta-method standard errors and a
#'   Wald confidence band for each predicted probability (requires a fit with
#'   var = TRUE). These SEs are computed numerically against
#'   compute_outcome_probs and combined with vcov(object), so they adapt
#'   automatically to k, the formulas, and the link.
#' @param level Confidence level for the band when se.fit = TRUE.
#' @param ... Ignored.
#' @return If se.fit = FALSE: a matrix (type = "probs") or factor
#'   (type = "class"). If se.fit = TRUE: a list with fit, se.fit, lower, upper
#'   (each an n x k matrix).
#' @export
predict.ordcure <- function(object, newdata,
                            type = c("probs", "class"),
                            component = c("uncured", "cured"),
                            delta = object$vars$delta,
                            se.fit = FALSE, level = 0.95, ...) {
  
  type      <- match.arg(type)
  component <- match.arg(component)
  
  if (missing(newdata) || is.null(newdata))
    stop("newdata is required")
  
  if (!delta %in% names(newdata))
    stop("newdata must contain the event-indicator column '", delta, "'")
  
  ## predicted probabilities for one parameter vector (returns n x k)
  .probs_for <- function(par.list) {
    pr <- compute_outcome_probs(
      par = par.list, data = newdata, delta = delta,
      k = object$k, outcome.model = object$outcome.model)
    if (component == "cured") pr$D0 else pr$D1
  }
  
  P <- .probs_for(object$par.list)
  colnames(P) <- object$ylevels
  
  if (type == "class") {
    cls <- object$ylevels[max.col(P, ties.method = "first")]
    return(factor(cls, levels = object$ylevels, ordered = TRUE))
  }
  
  if (!se.fit) return(P)
  
  ## ---- delta-method SEs (numeric, model-agnostic) ----
  if (is.null(object$variance))
    stop("se.fit = TRUE requires a fit with var = TRUE")
  
  par0 <- unlist(object$par.list)
  V    <- object$variance$stacked.v.est        # ordered as unlist(par.list)
  
  ## flatten n x k -> length n*k vector of a perturbed parameter vector
  flat <- function(v) as.vector(.probs_for(utils::relist(v, object$par.list)))
  J  <- .num_jacobian(flat, par0)              # (n*k) x p
  se <- sqrt(pmax(0, rowSums((J %*% V) * J)))  # diag(J V J') without forming it
  SE <- matrix(se, nrow = nrow(P), dimnames = dimnames(P))
  
  z     <- qnorm(1 - (1 - level) / 2)
  lower <- pmax(0, P - z * SE)
  upper <- pmin(1, P + z * SE)
  list(fit = P, se.fit = SE, lower = lower, upper = upper)
}

## Central-difference Jacobian (base R; no numDeriv dependency for prediction).
## f: R^p -> R^m ; returns m x p matrix.
.num_jacobian <- function(f, x, rel = 1e-6) {
  f0 <- f(x)
  m  <- length(f0); p <- length(x)
  J  <- matrix(0, m, p)
  for (j in seq_len(p)) {
    dx <- rel * max(abs(x[j]), 1)
    xp <- x; xp[j] <- xp[j] + dx
    xm <- x; xm[j] <- xm[j] - dx
    J[, j] <- (f(xp) - f(xm)) / (2 * dx)
  }
  J
}