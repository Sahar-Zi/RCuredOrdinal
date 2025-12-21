variance.est <- function(est, outcome.model = c("PO", "ACAT"), data, survform, cureform, Tau, delta, V, lengths) {
  n  <- nrow(data)
  
  #### 1. Build parameter vectors ---------------------------------------------
  
  gamma_est <- c(est$d, est$Tau)
  
  theta_est <- c(est$a, est$e, est$b, est$c)
  
  full_est  <- c(theta_est, gamma_est)
  
  mm.b <- model.matrix(cureform, data)
  mm.beta <- model.matrix(update(survform, . ~ . - 1), data)

  ## 5. Allocate matrices ----
  p  <- length(gamma_est)
  q  <- length(theta_est)
  m  <- matrix(NA_real_, p, n)
  g  <- matrix(NA_real_, q, n)
  omega <- matrix(0, p+q, p+q)
  e2 <- matrix(0,  q, q)
  
  #### 2. First-stage analytic gradients mᵢ (vectorized) -----------------------
  obs.log.likelihood.grads.vec <- function(par, mm.D, mm.T, Tau, delta){
    
    kD <- ncol(mm.D)
    kT <- ncol(mm.T)
    
    b    <- par[1:kD]
    zeta <- par[kD + (1:kT)]
    a    <- par[length(par)-1]
    k    <- par[length(par)]
    
    Z.b <- mm.D %*% b
    pi  <- as.vector( exp(Z.b) / (1 + exp(Z.b)) )
    
    Z.z  <- mm.T %*% zeta
    log.s <- -( (Tau / k)^a * exp(Z.z) )
    
    grad.b <- matrix(0, n, kD)
    grad.zeta <- matrix(0, n, kT)
    grad.a <- numeric(n)
    grad.k <- numeric(n)
    
    # delta = 1
    idx1 <- which(delta == 1)
    if (length(idx1) > 0) {
      grad.b[idx1, ] <- (1 - pi[idx1]) * mm.D[idx1, , drop = FALSE]
      grad.a[idx1]   <- (1/a) + log(Tau[idx1]) - log(k) +
        (log(Tau[idx1]/k) * log.s[idx1])
      grad.k[idx1]   <- -(a/k) * (1 + log.s[idx1])
      grad.zeta[idx1,] <- (1 + log.s[idx1]) * mm.T[idx1, , drop = FALSE]
    }
    
    # delta = 0
    idx0 <- which(delta == 0)
    if (length(idx0) > 0) {
      num  <- 1 - exp(log.s[idx0])
      den  <- 1 - pi[idx0] + pi[idx0] * exp(log.s[idx0])
      frac <- num / den
      
      grad.b[idx0, ] <- frac * pi[idx0] * (1 - pi[idx0]) * mm.D[idx0, , drop = FALSE]
      grad.a[idx0]   <- (pi[idx0]/den) * (exp(log.s[idx0]) *
                                            log(Tau[idx0]/k) * log.s[idx0])
      grad.k[idx0]   <- (pi[idx0]/den) * (-exp(log.s[idx0]) * (a/k) * log.s[idx0])
      grad.zeta[idx0,] <- (pi[idx0]/den) * (exp(log.s[idx0]) * log.s[idx0]) *
        mm.T[idx0, , drop = FALSE]
    }
    
    cbind(grad.b, grad.zeta, grad.a, grad.k)
  }

  m <- t(obs.log.likelihood.grads.vec(
    par = gamma_est,
    mm.D = mm.b,
    mm.T = mm.beta,
    Tau  = data[[Tau]],
    delta = data[[delta]]))

  #### 3. Build vectorized pseudo-loglikelihood for g ------------------------
  
  vectorized_pseudo_loglik_g <- function(par, data, outcome.model,
                                         mm.D, mm.T, lengths,
                                         k.levels, delta, Tau, V, Hessian = F) {
    ## parameter unpacking ------------
    par.list <- list(
      a = par[1:(lengths[1])],
      e = par[lengths[1] + (1:(lengths[2]))],
      b = par[lengths[1] + lengths[2] + (1:(lengths[3]))],
      c = par[lengths[1] + lengths[2] + lengths[3] + (1:(lengths[4]))])
    if (lengths[2]==0) par.list$e <- NULL
    
    b <- par[sum(lengths) + (1:ncol(mm.D))]
    zeta <- par[sum(lengths) + ncol(mm.D) + (1:ncol(mm.T))]
    a <- par[length(par)-1]
    k <- par[length(par)]
    
    v.probs <- estimate.outcome.probabilities(
      par  = par.list,
      data = data,  # single row to df
      delta = delta, k = k.levels, outcome.model = outcome.model)
    
    idx <- cbind(seq_len(nrow(data)), data[,V])
    p.alpha      <- v.probs$D0[idx]
    p.beta.gamma <- v.probs$D1[idx]

    Z.b   <- as.vector(mm.D %*% b)
    log_pi <- Z.b - log1p(exp(Z.b)) #log(pi)
    pi <- exp(log_pi)
    
    Tau_i <- data[,Tau]
    Z.z   <- as.vector(mm.T %*% zeta)

    log_s <- -((Tau_i/k)^a * exp(Z.z))

    del <- data[,delta]
    out <- numeric(nrow(data))

    out[del == 1] <- log(p.beta.gamma[del == 1])

    out[del == 0] <- log(
      pi[del == 0] * exp(log_s[del == 0]) * p.beta.gamma[del == 0] +
        (1 - pi[del == 0]) * p.alpha[del == 0]
    )
    
    if(Hessian) return(as.numeric(sum(out)))
    return(out)
  }

  #### 4. Compute ALL g with one jacobian call -------------------------------
  g <- t(jacobian(vectorized_pseudo_loglik_g, x = full_est, data = data, outcome.model = outcome.model,
                  mm.D = mm.b, mm.T = mm.beta, lengths = lengths, k.levels = n_distinct(data[,V]),
                  delta = delta, Tau = Tau, V = V))[1:q,]

  G <- hessian(vectorized_pseudo_loglik_g, x = full_est, data = data, outcome.model = outcome.model,
               mm.D = mm.b, mm.T = mm.beta, lengths = lengths, k.levels = n_distinct(data[[V]]), 
               delta = delta, Tau = Tau, V = V, Hessian = T)

  ## First‑stage log‑likelihood Gradients for Hessian
  obs.log.likelihood <- function(par, mm.D, mm.T, Tau, delta){
    b <- par[1:ncol(mm.D)]
    zeta <- par[ncol(mm.D)+(1:ncol(mm.T))]
    a <- par[length(par)-1]
    k <- par[length(par)]

    Z.b <- (b %*% t(mm.D))
    log_pi <- Z.b - log1p(exp(Z.b))  # log(pi)
    pi <- exp(log_pi)

    Z.z <- (zeta %*% t(mm.T))
    log_h <- log(a / k) + (a - 1) * log(Tau / k) + Z.z
    log_s <- - (Tau / k)^a * exp(Z.z)

    return(sum(delta * (log_pi + log_h + log_s) +
                 (1 - delta) * log1p(-pi + pi * exp(log_s))))
  }

  # ## 3. Pseudo log‑likelihood h (whole data) ----
  # pseudo.log.like.h <- function(par, data, outcome.model, mm.D, mm.T, lengths, k.levels, delta, Tau, V) {
  # 
  #   ## parameter unpacking ------------
  #   par.list <- list(
  #     a = par[1:(lengths[1])],
  #     b = par[lengths[1] + (1:(lengths[2]))],
  #     c = par[lengths[1] + lengths[2] + (1:(lengths[3]))])
  # 
  #   if(length(lengths)==4) par.list$e <- par[(lengths[1]+lengths[2]+lengths[3]+1):sum(lengths)]
  # 
  #   b <- par[sum(lengths) + (1:ncol(mm.D))]
  #   zeta <- par[sum(lengths) + ncol(mm.D) + (1:ncol(mm.T))]
  #   a <- par[length(par)-1]
  #   k <- par[length(par)]
  # 
  #   v.probs <- estimate.outcome.probabilities(
  #     par  = par.list, data = data, delta = delta, k = k.levels, outcome.model = outcome.model)
  # 
  #   p.alpha      <- v.probs$D0[cbind(seq_len(nrow(data)), data[[V]])]
  #   p.beta.gamma <- v.probs$D1[cbind(seq_len(nrow(data)), data[[V]])]
  # 
  #   Z.b <- (b %*% t(mm.D))
  #   log_pi <- Z.b - log1p(exp(Z.b))  # log(pi)
  #   pi <- exp(log_pi)
  # 
  #   Z.z <- (zeta %*% t(mm.T))
  #   log_h <- log(a / k) + (a - 1) * log(data[[Tau]] / k) + Z.z
  #   log_s <- - (data[[Tau]] / k)^a * exp(Z.z)
  # 
  #   del <- data[[delta]]
  #   sum(log(p.beta.gamma[del == 1])) +
  #     sum(log(pi[del == 0] * exp(log_s[del == 0]) * p.beta.gamma[del == 0] +
  #               (1 - pi[del == 0]) * p.alpha[del == 0]))
  # }

  ## 6. Loop over rows (still row‑wise because grad() needs it) ----
  # for (i in seq_len(n)) {
  #   if (i %% 10 == 0) cat("Row", i, "of", n, "\n")
  #   m[, i] <- obs.log.likelihood.grads(par = gamma_est, mm.D = mm.b[i,], mm.T = mm.beta[i,],
  #                                      Tau = data[[Tau]][i], delta = data[[delta]][i])
  # 
  #   g[, i] <- grad(pseudo.log.like.g, x = full_est, data = data[i, ], outcome.model = outcome.model,
  #                  mm.D = mm.b[i,], mm.T = mm.beta[i,], lengths = lengths, k.levels = n_distinct(data[[V]]), delta = delta, Tau = Tau, V = V)[1:q]
  #   }
  
  stacked <- rbind(m,g)
  
  ## 7. Hessians (whole data) ----
  M <- hessian(obs.log.likelihood, x = gamma_est, mm.D = mm.b, mm.T = mm.beta,
               Tau = data[[Tau]], delta = data[[delta]])

  psi  <- -solve(M, m)
  v <- tcrossprod(psi)

  # G <- hessian(pseudo.log.like.h, x = full_est, data = data, outcome.model = outcome.model,
               # mm.D = mm.b, mm.T = mm.beta, lengths = lengths, k.levels =n_distinct(data[[V]]), delta = delta, Tau = Tau, V = V)
  G.theta  <- G[1:q, 1:q]
  G.gamma  <- G[1:q, (q + 1):(q + p)]

  ## 8. Sandwich variance formulas ----

  #v <- solve(M, e1 / n) %*% t(solve(M))

  g.psi  <- g + G.gamma %*% psi
  
  for (i in seq_len(n)) {
    e2 <- e2 + tcrossprod(g.psi[, i])
    omega <- omega + tcrossprod(stacked[, i])
  }

  #V <- solve(G.theta, e2) %*% t(solve(G.theta))
  
  cov.mat2 <- matrix(0, p+q, p+q)
  cov.mat2[1:p,1:p] <- M
  cov.mat2[(p+1):(p+q),1:p] <- G.gamma
  cov.mat2[(p+1):(p+q),(p+1):(p+q)] <- G.theta
  #cov.mat2[1:q,1:q] <- G.theta
  #cov.mat2[1:q,(q+1):(q+p)] <- G.gamma
  #cov.mat2[(q+1):(q+p),(q+1):(q+p)] <- M
  cov.mat2 <- solve(cov.mat2)
  
  stacked.v.est <- cov.mat2 %*% omega %*% t(cov.mat2)
  
  rownames(stacked.v.est) <- colnames(stacked.v.est) <- names(c(gamma_est, theta_est))
  
  return(list(stacked.v.est = stacked.v.est,
              G.tilde.inv = cov.mat2))
}
