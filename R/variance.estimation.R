# ============================================================
# Sandwich (stacked estimating-equations) variance for ordcure().
#
# Accounts for first-stage (cure model) uncertainty when estimating
# the variance of the second-stage (ordinal) parameters. Internal
# (not exported); called by ordcure() when var = TRUE.
#
# NOTE on the first-stage Hessian: two implementations are available.
#   M  - numeric, via numDeriv::hessian (used below).
#   M2 - analytic, via obs.log.likelihood.hess() (computed but not
#        wired in). To switch, replace `M` with `M2` in the `psi`
#        line and the `cov.mat2[1:p, 1:p]` assignment.
# ============================================================

#' @noRd
variance_est <- function(est, outcome.model = c("PO", "ACAT"),
                         data, survform, cureform, Tau, delta, V, lengths) {
  
  outcome.model <- match.arg(outcome.model)
  n <- nrow(data)
  
  #### 1. Parameter vectors --------------------------------------------------
  
  gamma_est <- c(est$d, est$Tau)              # first-stage params
  theta_est <- c(est$a, est$e, est$b, est$c)  # second-stage params
  full_est  <- c(theta_est, gamma_est)
  
  mm.b    <- model.matrix(cureform, data)
  mm.beta <- model.matrix(update(survform, . ~ . - 1), data)
  
  p     <- length(gamma_est)
  q     <- length(theta_est)
  omega <- matrix(0, p + q, p + q)
  e2    <- matrix(0, q, q)
  
  #### 2. First-stage analytic gradients m (vectorized) ----------------------
  
  obs.log.likelihood.grads.vec <- function(par, mm.D, mm.T, Tau, delta) {
    
    kD <- ncol(mm.D)
    kT <- ncol(mm.T)
    
    b    <- par[1:kD]
    zeta <- par[kD + (1:kT)]
    a    <- par[length(par) - 1]
    k    <- par[length(par)]
    
    Z.b <- mm.D %*% b
    pi  <- as.vector(exp(Z.b) / (1 + exp(Z.b)))
    
    Z.z   <- mm.T %*% zeta
    log.s <- -((Tau / k)^a * exp(Z.z))
    
    grad.b    <- matrix(0, n, kD)
    grad.zeta <- matrix(0, n, kT)
    grad.a    <- numeric(n)
    grad.k    <- numeric(n)
    
    # delta = 1
    idx1 <- which(delta == 1)
    if (length(idx1) > 0) {
      grad.b[idx1, ]   <- (1 - pi[idx1]) * mm.D[idx1, , drop = FALSE]
      grad.a[idx1]     <- (1 / a) + log(Tau[idx1]) - log(k) +
        (log(Tau[idx1] / k) * log.s[idx1])
      grad.k[idx1]     <- -(a / k) * (1 + log.s[idx1])
      grad.zeta[idx1, ] <- (1 + log.s[idx1]) * mm.T[idx1, , drop = FALSE]
    }
    
    # delta = 0
    idx0 <- which(delta == 0)
    if (length(idx0) > 0) {
      num  <- 1 - exp(log.s[idx0])
      den  <- 1 - pi[idx0] + pi[idx0] * exp(log.s[idx0])
      frac <- num / den
      
      grad.b[idx0, ]   <- frac * pi[idx0] * (1 - pi[idx0]) * mm.D[idx0, , drop = FALSE]
      grad.a[idx0]     <- (pi[idx0] / den) * (exp(log.s[idx0]) *
                                                log(Tau[idx0] / k) * log.s[idx0])
      grad.k[idx0]     <- (pi[idx0] / den) * (-exp(log.s[idx0]) * (a / k) * log.s[idx0])
      grad.zeta[idx0, ] <- (pi[idx0] / den) * (exp(log.s[idx0]) * log.s[idx0]) *
        mm.T[idx0, , drop = FALSE]
    }
    
    cbind(grad.b, grad.zeta, grad.a, grad.k)
  }
  
  m <- t(obs.log.likelihood.grads.vec(
    par   = gamma_est,
    mm.D  = mm.b,
    mm.T  = mm.beta,
    Tau   = data[[Tau]],
    delta = data[[delta]]
  ))
  
  #### 3. Vectorized pseudo-loglikelihood (per-obs, for g) --------------------
  
  vectorized_pseudo_loglik_g <- function(par, data, outcome.model,
                                         mm.D, mm.T, lengths,
                                         k.levels, delta, Tau, V, Hessian = FALSE) {
    
    ## parameter unpacking
    par.list <- list(
      a = par[1:(lengths[1])],
      e = par[lengths[1] + (1:(lengths[2]))],
      b = par[lengths[1] + lengths[2] + (1:(lengths[3]))],
      c = par[lengths[1] + lengths[2] + lengths[3] + (1:(lengths[4]))]
    )
    if (lengths[2] == 0) par.list$e <- NULL
    
    b    <- par[sum(lengths) + (1:ncol(mm.D))]
    zeta <- par[sum(lengths) + ncol(mm.D) + (1:ncol(mm.T))]
    a    <- par[length(par) - 1]
    k    <- par[length(par)]
    
    v.probs <- compute_outcome_probs(
      par = par.list, data = data,
      delta = delta, k = k.levels, outcome.model = outcome.model
    )
    
    idx          <- cbind(seq_len(nrow(data)), data[[V]])
    p.alpha      <- v.probs$D0[idx]
    p.beta.gamma <- v.probs$D1[idx]
    
    Z.b    <- as.vector(mm.D %*% b)
    log_pi <- Z.b - log1p(exp(Z.b))   # log(pi)
    pi     <- exp(log_pi)
    
    Tau_i <- data[[Tau]]
    Z.z   <- as.vector(mm.T %*% zeta)
    log_s <- -((Tau_i / k)^a * exp(Z.z))
    
    del <- data[[delta]]
    out <- numeric(nrow(data))
    
    out[del == 1] <- log(p.beta.gamma[del == 1])
    out[del == 0] <- log(
      pi[del == 0] * exp(log_s[del == 0]) * p.beta.gamma[del == 0] +
        (1 - pi[del == 0]) * p.alpha[del == 0]
    )
    
    if (Hessian) return(as.numeric(sum(out)))
    out
  }
  
  #### 4. g (per-obs gradient) and G (joint Hessian) -------------------------
  
  k.levels <- length(unique(data[[V]]))
  
  g <- t(numDeriv::jacobian(
    vectorized_pseudo_loglik_g, x = full_est, data = data,
    outcome.model = outcome.model, mm.D = mm.b, mm.T = mm.beta,
    lengths = lengths, k.levels = k.levels,
    delta = delta, Tau = Tau, V = V
  ))[1:q, ]
  
  G <- numDeriv::hessian(
    vectorized_pseudo_loglik_g, x = full_est, data = data,
    outcome.model = outcome.model, mm.D = mm.b, mm.T = mm.beta,
    lengths = lengths, k.levels = k.levels,
    delta = delta, Tau = Tau, V = V, Hessian = TRUE
  )
  
  #### 5. First-stage log-likelihood (for numeric Hessian M) -----------------
  
  obs.log.likelihood <- function(par, mm.D, mm.T, Tau, delta) {
    b    <- par[1:ncol(mm.D)]
    zeta <- par[ncol(mm.D) + (1:ncol(mm.T))]
    a    <- par[length(par) - 1]
    k    <- par[length(par)]
    
    Z.b    <- (b %*% t(mm.D))
    log_pi <- Z.b - log1p(exp(Z.b))   # log(pi)
    pi     <- exp(log_pi)
    
    Z.z   <- (zeta %*% t(mm.T))
    log_h <- log(a / k) + (a - 1) * log(Tau / k) + Z.z
    log_s <- -(Tau / k)^a * exp(Z.z)
    
    sum(delta * (log_pi + log_h + log_s) +
          (1 - delta) * log1p(-pi + pi * exp(log_s)))
  }
  
  M <- numDeriv::hessian(
    obs.log.likelihood, x = gamma_est, mm.D = mm.b, mm.T = mm.beta,
    Tau = data[[Tau]], delta = data[[delta]]
  )
  
  #### 6. First-stage analytic Hessian M2 (alternative to M) ------------------
  
  obs.log.likelihood.hess <- function(par, mm.D, mm.T, Tau, delta) {
    
    kD   <- ncol(mm.D)
    kT   <- ncol(mm.T)
    kpar <- kD + kT + 2
    
    b    <- par[1:kD]
    zeta <- par[kD + (1:kT)]
    a    <- par[length(par) - 1]
    k    <- par[length(par)]
    
    Z.b   <- mm.D %*% b
    pi    <- as.vector(exp(Z.b) / (1 + exp(Z.b)))
    
    Z.z   <- mm.T %*% zeta
    log.s <- -((Tau / k)^a * exp(Z.z))
    s     <- exp(log.s)
    
    H <- matrix(0, kpar, kpar)
    
    idx_b <- 1:kD
    idx_z <- kD + (1:kT)
    idx_a <- kD + kT + 1
    idx_k <- kD + kT + 2
    
    # ---- delta = 1 ----
    idx1 <- which(delta == 1)
    if (length(idx1) > 0) {
      
      pi1        <- pi[idx1]
      ls1        <- log.s[idx1]            # log(s)
      log_tau_k1 <- log(Tau[idx1] / k)     # log(tau/k)
      xD1        <- mm.D[idx1, , drop = FALSE]
      xT1        <- mm.T[idx1, , drop = FALSE]
      
      # b-b
      w_bb <- -pi1 * (1 - pi1)
      H[idx_b, idx_b] <- H[idx_b, idx_b] + t(xD1) %*% (w_bb * xD1)
      
      # zeta-zeta
      w_zz <- ls1
      H[idx_z, idx_z] <- H[idx_z, idx_z] + t(xT1) %*% (w_zz * xT1)
      
      # a-a
      H[idx_a, idx_a] <- H[idx_a, idx_a] + sum(-1 / a^2 + log_tau_k1^2 * ls1)
      
      # k-k
      H[idx_k, idx_k] <- H[idx_k, idx_k] + sum((a / k^2) * (1 + ls1) - (a / k)^2 * ls1)
      
      # zeta-a
      w_za <- log_tau_k1 * ls1
      za   <- colSums(w_za * xT1)
      H[idx_z, idx_a] <- H[idx_z, idx_a] + za
      H[idx_a, idx_z] <- H[idx_a, idx_z] + za
      
      # zeta-k
      w_zk <- -(a / k) * ls1
      zk   <- colSums(w_zk * xT1)
      H[idx_z, idx_k] <- H[idx_z, idx_k] + zk
      H[idx_k, idx_z] <- H[idx_k, idx_z] + zk
      
      # a-k
      H[idx_a, idx_k] <- H[idx_a, idx_k] + sum(-(1 / k) * (1 + ls1) - (a / k) * log_tau_k1 * ls1)
      H[idx_k, idx_a] <- H[idx_k, idx_a] + sum(-(1 / k) * (1 + ls1) - (a / k) * log_tau_k1 * ls1)
    }
    
    # ---- delta = 0 ----
    idx0 <- which(delta == 0)
    if (length(idx0) > 0) {
      
      pi0        <- pi[idx0]
      ls0        <- log.s[idx0]
      s0         <- s[idx0]
      log_tau_k0 <- log(Tau[idx0] / k)
      xD0        <- mm.D[idx0, , drop = FALSE]
      xT0        <- mm.T[idx0, , drop = FALSE]
      
      D   <- 1 - pi0 + pi0 * s0       # denominator
      psl <- pi0 * s0 * ls0           # pi*s*log(s)
      CF  <- (ls0 + 1) - psl / D      # common factor
      
      # b-b
      w_bb <- pi0 * (1 - pi0) * ((1 - s0) * (1 - 2 * pi0) / D -
                                   pi0 * (1 - pi0) * (1 - s0)^2 / D^2)
      H[idx_b, idx_b] <- H[idx_b, idx_b] + t(xD0) %*% (w_bb * xD0)
      
      # zeta-zeta
      w_zz <- psl / D * CF
      H[idx_z, idx_z] <- H[idx_z, idx_z] + t(xT0) %*% (w_zz * xT0)
      
      # a-a
      H[idx_a, idx_a] <- H[idx_a, idx_a] + sum(psl / D * log_tau_k0^2 * CF)
      
      # k-k
      H[idx_k, idx_k] <- H[idx_k, idx_k] +
        sum(psl * (a / k)^2 / D * (-CF) + psl / D * (a / k^2))
      
      # zeta-a
      w_za <- psl * log_tau_k0 / D * CF
      za   <- colSums(w_za * xT0)
      H[idx_z, idx_a] <- H[idx_z, idx_a] + za
      H[idx_a, idx_z] <- H[idx_a, idx_z] + za
      
      # zeta-k
      w_zk <- -psl * (a / k) / D * CF
      zk   <- colSums(w_zk * xT0)
      H[idx_z, idx_k] <- H[idx_z, idx_k] + zk
      H[idx_k, idx_z] <- H[idx_k, idx_z] + zk
      
      # a-k
      ak <- sum(-psl * log_tau_k0 * (a / k) / D * CF - psl / (D * k))
      H[idx_a, idx_k] <- H[idx_a, idx_k] + ak
      H[idx_k, idx_a] <- H[idx_k, idx_a] + ak
      
      # b-zeta
      w_bz <- pi0 * (1 - pi0) * (-s0 * ls0 / D - (1 - s0) * psl / D^2)
      bz   <- t(xD0) %*% (w_bz * xT0)
      H[idx_b, idx_z] <- H[idx_b, idx_z] + bz
      H[idx_z, idx_b] <- H[idx_z, idx_b] + t(bz)
      
      # b-a
      w_ba <- -pi0 * (1 - pi0) * s0 * ls0 * log_tau_k0 / D * (1 + (1 - s0) * pi0 / D)
      ba   <- colSums(w_ba * xD0)
      H[idx_b, idx_a] <- H[idx_b, idx_a] + ba
      H[idx_a, idx_b] <- H[idx_a, idx_b] + ba
      
      # b-k
      w_bk <- pi0 * (1 - pi0) * s0 * ls0 * (a / k) / D * (1 + (1 - s0) * pi0 / D)
      bk   <- colSums(w_bk * xD0)
      H[idx_b, idx_k] <- H[idx_b, idx_k] + bk
      H[idx_k, idx_b] <- H[idx_k, idx_b] + bk
    }
    
    H
  }
  
  M2 <- t(obs.log.likelihood.hess(
    par   = gamma_est,
    mm.D  = mm.b,
    mm.T  = mm.beta,
    Tau   = data[[Tau]],
    delta = data[[delta]]
  ))
  
  #### 7. Sandwich assembly ---------------------------------------------------
  
  stacked <- rbind(m, g)
  
  psi <- -solve(M, m)   # influence of first stage on its own params
  
  G.theta <- G[1:q, 1:q]
  G.gamma <- G[1:q, (q + 1):(q + p)]
  
  ## g.psi / e2 are components of an alternative middle-matrix formula,
  ## kept for reference; the returned variance uses `omega` below.
  g.psi <- g + G.gamma %*% psi
  
  for (i in seq_len(n)) {
    e2    <- e2 + tcrossprod(g.psi[, i])
    omega <- omega + tcrossprod(stacked[, i])
  }
  
  ## block "bread" matrix (lower-triangular: first stage does not depend
  ## on second-stage parameters)
  cov.mat2 <- matrix(0, p + q, p + q)
  cov.mat2[1:p, 1:p]                     <- M
  cov.mat2[(p + 1):(p + q), 1:p]         <- G.gamma
  cov.mat2[(p + 1):(p + q), (p + 1):(p + q)] <- G.theta
  
  cov.mat2 <- solve(cov.mat2)
  
  stacked.v.est <- cov.mat2 %*% omega %*% t(cov.mat2)
  
  rownames(stacked.v.est) <- colnames(stacked.v.est) <-
    names(c(gamma_est, theta_est))
  
  list(
    stacked.v.est = stacked.v.est,
    G.tilde.inv   = cov.mat2
  )
}