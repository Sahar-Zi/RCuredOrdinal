# ============================================================
# Stage 1: mixture cure model (logistic incidence + Weibull PH
# latency), fitted by EM. Produces the posterior "uncured" weights
# that seed the stage-2 ordinal pseudo-likelihood.
#
# Internal (not exported); called by ordcure().
# ============================================================

#' Fit the first-stage Weibull mixture cure model by EM
#'
#' @param data A data frame.
#' @param survform Latency (survival) formula, e.g. `Tau ~ Z1 + Z2`.
#' @param cureform Incidence (cure) formula, e.g. `delta ~ Z1 + Z2`.
#' @param Tau Character name of the (observed) time column.
#' @param delta Character name of the event-indicator column.
#' @param tol Convergence tolerance on the parameter L2 step.
#' @param maxit Maximum number of EM iterations.
#' @param verbose If `TRUE`, print per-iteration progress.
#'
#' @return A list with the parameter vector `par` (named: cure
#'   coefficients, latency coefficients, Weibull `shape`, `scale`),
#'   posterior weights `weights` (`w0`, `w1`), the log-likelihood
#'   `loglik`, and the iteration count `iter`.
#' @noRd
fit_cure_weibull_em <- function(
    data, survform, cureform, Tau, delta,
    tol = 1e-8, maxit = 200, verbose = FALSE
) {

  n  <- nrow(data)
  i0 <- data[[delta]] == 0

  ## ---- design matrices (cached) ----
  Xc <- model.matrix(cureform, data)                       # incidence
  Xt <- model.matrix(update(survform, . ~ . - 1), data)    # latency (no intercept)

  ## ---- initial values ----
  d.CC <- data[data[[delta]] == 1, , drop = FALSE]

  surv_init <- as.formula(
    paste0("Surv(time=", Tau, ", event=rep(1, nrow(d.CC))) ~ ",
           paste(colnames(Xt), collapse = " + "))
  )

  fit_lat  <- eha::phreg(surv_init, data = d.CC)
  fit_cure <- glm(cureform, family = binomial, data = data)

  ## eha::phreg stores, as the last two coefficients, log(scale) and
  ## log(shape); recover the Weibull shape (a) and scale (k) below.
  b    <- coef(fit_cure)
  beta <- fit_lat$coef[1:(length(fit_lat$coef) - 2)]
  a    <- exp(fit_lat$coef[length(fit_lat$coef)])        # Weibull shape
  k    <- exp(-fit_lat$coef[length(fit_lat$coef) - 1])   # Weibull scale

  par <- c(b, beta, a, k)

  ## ---- log-likelihoods (negative, weighted) ----
  loglik_cure <- function(b, w) {
    eta <- Xc %*% b
    -sum(w * eta - log1p(exp(eta)))
  }

  loglik_lat <- function(theta, w) {
    beta <- theta[1:(length(theta) - 2)]
    a    <- theta[length(theta) - 1]
    k    <- theta[length(theta)]

    eta <- Xt %*% beta
    t   <- data[[Tau]]

    ll <- data[[delta]] * (eta + log(a) + (a - 1) * log(t) - a * log(k)) -
      exp(eta + a * log(t) - a * log(k))
    -sum(w * ll)
  }

  ## ---- EM loop ----
  loglik_prev <- Inf

  for (iter in seq_len(maxit)) {

    ## E-step: posterior probability of being uncured for censored rows.
    pi_cure <- plogis(Xc %*% b)
    surv    <- exp(-((data[[Tau]] / k)^a) * exp(Xt %*% beta))

    w <- rep(1, n)
    w[i0] <- (pi_cure[i0] * surv[i0]) /
      (pi_cure[i0] * surv[i0] + 1 - pi_cure[i0])

    ## M-step: weighted incidence and latency fits.
    fit1 <- optim(b, loglik_cure, w = w, method = "BFGS")
    fit2 <- optim(
      c(beta, a, k), loglik_lat, w = w,
      method = "L-BFGS-B",
      lower = c(rep(-Inf, ncol(Xt)), 1e-6, 1e-6)
    )

    b    <- fit1$par
    beta <- fit2$par[1:(length(fit2$par) - 2)]
    a    <- fit2$par[length(fit2$par) - 1]
    k    <- fit2$par[length(fit2$par)]

    loglik <- fit1$value + fit2$value
    diff.L <- abs(loglik_prev - loglik)
    diff   <- sqrt(sum((par - c(b, beta, a, k))^2))
    par    <- c(b, beta, a, k)

    if (verbose)
      message("iter=", iter, " log.L=", -loglik,
              " diff=", diff, " diff.L=", diff.L)

    if (diff < tol) break
    loglik_prev <- loglik
  }

  ## ---- posterior weights ----
  ## Recompute at the converged parameters so the returned weights are
  ## consistent with the final fit (the loop's last E-step used the
  ## previous M-step's values).
  pi_cure <- as.vector(plogis(Xc %*% b))
  surv    <- as.vector(exp(-((data[[Tau]] / k)^a) * exp(Xt %*% beta)))

  wdf <- data.frame(
    w0 = ifelse(i0, 1 - pi_cure, 1),
    w1 = ifelse(i0, pi_cure * surv, 1)
  )

  ## ---- parameter names ----
  cure_terms <- colnames(Xc)
  cure_terms[cure_terms == "(Intercept)"] <- "0"
  cure_names <- paste0("b.", cure_terms)

  lat_names <- paste0("zeta.", colnames(Xt))

  names(par) <- c(cure_names, lat_names, "shape", "scale")

  list(
    par     = par,
    weights = wdf,
    loglik  = -loglik,
    iter    = iter
  )
}
