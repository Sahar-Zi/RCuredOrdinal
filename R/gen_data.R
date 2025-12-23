#' Generate synthetic data for ordinal cure models
#'
#' Simulates covariates, cure status, survival times, and ordinal outcomes
#' under the ordinal cure model with Weibull latency.
#'
#' @param n Integer. Sample size.
#' @param k Integer. Number of ordinal categories.
#' @param par List. Model parameters (same structure as ordcure()).
#' @param outcome.model Character. "PO" or "ACAT".
#'
#' @return A data.frame containing simulated data.
#' @export
gen_demo_data <- function(n, k, par, outcome.model = c("PO", "ACAT")) {
  
  outcome.model <- match.arg(outcome.model)
  
  ## ============================================================
  ## 1. Generate baseline covariates
  ## ============================================================
  z0 <- 1
  u <- runif(n, 0, 1)
  z1 <- runif(n, 0, 1)
  z2 <- (z1 < 0.5 & u > 2/3) + (z1 > 0.5 & u < 2/3)
  Z <- cbind(z0,z1,z2)  # intercept included
  
  ## ============================================================
  ## 2. Cure indicator D
  ## ============================================================
  D <- rbinom(n, 1, plogis(Z %*% par$d))
  
  ## ============================================================
  ## 3. Latent survival time Tau
  ## ============================================================
  shape <- par$Tau.shape
  scale <- par$Tau.scale
  
  linpred_T <- Z %*% c(0, par$Tau)  # no intercept
  scale_i  <- exp(-linpred_T / shape) / scale
  
  Tau_latent <- ifelse(
    D == 1,
    rweibull(n, shape = shape, scale = scale_i),
    Inf
  )
  
  ## ============================================================
  ## 4. Generate censoring / inspection time R
  ## ============================================================
  u <- runif(n)
  
  gen_age <- function(u) {
    if (u <= 0.8) {
      1.25 * u + 0.5
    } else {
      uniroot(
        function(r) r^2 - 4 * r + 2.75 + 1.25 * u,
        c(1.5, 2)
      )$root
    }
  }
  
  rz1 <- sapply(u, gen_age) # age distribution for Z2=1
  rz0 <- runif(n, 0.3, 2) # age distribution for Z2=0
  R <- ifelse(z2, rz1, rz0)
  
  ## ============================================================
  ## 5. Observed time and event indicator
  ## ============================================================
  Tau_obs <- pmin(Tau_latent, R)
  delta   <- as.integer(Tau_latent < R & D == 1)
  
  data <- data.frame(
    Z1 = z1,
    Z2 = z2,
    D  = D,
    Tau = Tau_obs,
    R  = R,
    delta = delta
  )
  
  ## ============================================================
  ## 6. Generate ordinal outcome V
  ## ============================================================
  probs <- compute_outcome_probs(
    par = par,
    data = data,
    delta = "delta",
    k = k,
    outcome.model = outcome.model
  )
  
  ## cumulative probabilities
  Uv <- runif(n)
  cum_probs <- lapply(probs, function(df){t(apply(df, 1, cumsum))})
  cum_probs_by_D <- cum_probs$D0
  cum_probs_by_D[D==1,] <- cum_probs$D1[D==1,]
  
  V <- sapply(1:n, function(i) {
    which(cum_probs_by_D[i,] > Uv[i])[1]  # Find first column index where v[i,] > u[i]
  })

  data$V <- ordered(V)
  
  data
}

