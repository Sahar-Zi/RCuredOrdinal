############################################################
## Simulation utilities
## Ordinal Cure Model with Censored Covariates
############################################################

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