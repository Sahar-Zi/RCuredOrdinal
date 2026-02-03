# ============================================================
# Stage-2 pseudo log-likelihood for ordinal cure model
# ============================================================
# Computes the negative pseudo log-likelihood used in the
# second-stage optimization of ordcure().
#
# The likelihood combines:
#  - observed outcomes (delta = 1)
#  - mixture contributions for censored observations (delta = 0)
#
# Supports PO and ACAT ordinal outcome models.
# ============================================================

expand_po_thresholds <- function(par, lengths, k) {
  if (k <= 2) return(par)
  
  idx <- cumsum(c(0, lengths))
  
  for (block in c(1, 3, 4)) {  # a, b, c blocks
    start <- idx[block] + 1
    end   <- start + k - 2
    par[start:end] <- cumsum(par[start:end])
  }
  
  par
}

unpack_stage2_parameters <- function(par, lengths) {
  
  idx <- cumsum(c(0, lengths))
  
  list(
    a = par[(idx[1] + 1):idx[2]],
    e = if (lengths[2] > 0) par[(idx[2] + 1):idx[3]] else NULL,
    b = par[(idx[3] + 1):idx[4]],
    c = par[(idx[4] + 1):idx[5]]
  )
}

negloglik_stage2 <- function(
    par,
    lengths,
    k,
    delta,
    response,
    data,
    weights,
    outcome.model = c("PO", "ACAT")
) {
  
  outcome.model <- match.arg(outcome.model)
  
  ## ---- expand PO thresholds if needed ----
  par_expanded <- par
  if (outcome.model == "PO" && k > 2) {
    par_expanded <- expand_po_thresholds(par_expanded, lengths, k)
  }
  
  ## ---- unpack parameters ----
  theta <- unpack_stage2_parameters(par_expanded, lengths)

  ## ---- compute outcome probabilities ----
  probs <- compute_outcome_probs(
    par = theta,
    data = data,
    delta = delta,
    k = k,
    outcome.model = outcome.model
  )
  
  ## ---- extract observed category probabilities ----
  y <- data[[response]]
  idx <- cbind(seq_len(nrow(data)), y)

  p_D0 <- probs$D0[idx]  # P(Y | delta = 0)
  p_D1 <- probs$D1[idx]  # P(Y | delta = 1)
  
  ## ---- log-likelihood contributions ----
  ll_obs <- sum(log(p_D1[data[[delta]] == 1]))
  
  ll_mix <- sum(
    log(
      weights$w1 * p_D1 + weights$w0 * p_D0
    )[data[[delta]] == 0]
  )
  
  ## ---- return negative log-likelihood ----
  -(ll_obs + ll_mix)
}
