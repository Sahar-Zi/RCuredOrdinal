gen.data <- function(n, p.z2 = 0.5, D.params, T.params, 
                     V.alpha0, V.beta0, V.gamma0, V.alphaZ, V.gammaZ,
                     V.alphaR, V.gammaR, V.gammaT, V.gammaTR, V.model="PO"){
  
  start.time <- Sys.time()
  
  
  if(length(D.params) != 3){stop("D Creation Error: params and variables are not the same length")}
  if(length(T.params) != 5){stop("T Creation Error: params and variables are not the same length")}
  if(any(T.params[1:2] <= 0)){stop("T Creation Error: scale or shape less than or equal to zero")}
  
  # Generation of Z
  z0 <- 1
  u <- runif(n, 0, 1)
  z1 <- runif(n, 0, 1)
  z2 <- (z1 < 0.5 & u > p.z2) + (z1 > 0.5 & u < p.z2)
  Z <- cbind(z0,z1,z2)
  
  # Generation of D
  d <- rbinom(n, 1, plogis(Z %*% D.params))
  
  # Generation of T
  scale.T <- par$scale
  shape.T <- par$shape
  scale.tau <- exp(-(Z %*% c(0,par$Tau))/shape.T)/scale.T
  tau <- ifelse(d==1, rweibull(n, shape = shape.T, scale = scale.tau), Inf)
  
  # Generation of R
  u <- runif(n, 0, 1)
  gen.age <- function(u) ifelse(u <= 0.8, 
                                  1.25 * u + 0.5, 
                                  uniroot(function(r, u) r^2 - 4 * r + 2.75 + 1.25 * u, c(1.5, 2), u = u)$root)
  rz1 <- sapply(u, gen.age) # age distribution for Z2=1
  rz0 <- runif(n, 0.3, 2) # age distribution for Z2=0
  r <- ifelse(z2, rz1, rz0)
  
  # Generation of delta
  delta <- (tau<r)*(d==1)+0
  
  simulated_data <- data.frame(Z1 = z1, Z2 = z2, D = d, Tau = pmin(tau, r), R = r, TR=ifelse(delta, tau*r, Inf), delta = delta)

  # Generate V
  if(V.model=="PO"){
    v.probs <- probs.po(alpha0 = V.alpha0, alphaZ = V.alphaZ, alphaR = V.alphaR, 
                        beta0 = V.beta0, gammaZ = V.gammaZ,
                        gamma0 = V.gamma0, gammaR = V.gammaR,
                        gammaT = V.gammaT, gammaTR = V.gammaTR,
                        rangeZ = 1:2, rangeR = 5, rangeT = 4, rangeTR = 6, dat = as.matrix(simulated_data))
  }
  if(V.model=="ACAT"){
    v.probs <- probs.acat(alpha0 = V.alpha0, alphaZ = V.alphaZ, alphaR = V.alphaR, 
                        beta0 = V.beta0, gammaZ = V.gammaZ,
                        gamma0 = V.gamma0, gammaR = V.gammaR,
                        gammaT = V.gammaT, gammaTR = V.gammaTR,
                        rangeZ = 1:2, rangeR = 5, rangeT = 4, rangeTR = 6, dat = as.matrix(simulated_data))
  }
  
  u <- runif(n, 0, 1)
  v.cummulative.probs <- lapply(v.probs, function(df){t(apply(df, 1, cumsum))})
  v.cummulative.by.D <- v.cummulative.probs$D0
  v.cummulative.by.D[d==1,] <- v.cummulative.probs$D1[d==1,]  # If d is 1, take rows from D1 else D0
  
  v <- sapply(1:n, function(i) {
    which(v.cummulative.by.D[i,] > u[i])[1]  # Find first column index where v[i,] > u[i]
  })
  
  #print(cbind(d,v.cummulative.by.D,u,v))
  
  simulated_data <- cbind(simulated_data, V=ordered(v))
  
  # condition1.1 <- sum(simulated_data$D == 0 & simulated_data$delta == 0 & simulated_data$V == 1)
  # condition1.2 <- sum(simulated_data$D == 0 & simulated_data$delta == 0 & simulated_data$V == 2)
  # condition1.3 <- sum(simulated_data$D == 0 & simulated_data$delta == 0 & simulated_data$V == 3)
  # condition2.1 <- sum(simulated_data$D == 1 & simulated_data$delta == 0 & simulated_data$V == 1)
  # condition2.2 <- sum(simulated_data$D == 1 & simulated_data$delta == 0 & simulated_data$V == 2)
  # condition2.3 <- sum(simulated_data$D == 1 & simulated_data$delta == 0 & simulated_data$V == 3)
  # condition3.1 <- sum(simulated_data$D == 1 & simulated_data$delta == 1 & simulated_data$V == 1)
  # condition3.2 <- sum(simulated_data$D == 1 & simulated_data$delta == 1 & simulated_data$V == 2)
  # condition3.3 <- sum(simulated_data$D == 1 & simulated_data$delta == 1 & simulated_data$V == 3)
  # print(c(condition1.1,condition1.2,condition1.3,
  #         condition2.1,condition2.2,condition2.3,
  #         condition3.1,condition3.2,condition3.3))
  # 
  #print(paste0("TIME - Generating Data: ",Sys.time()-start.time))
  return(simulated_data)
}

gen.demo.data <- function(n, k, par, outcome.model=c("PO","ACAT")){
  # --------------------------------------------------
  # what if par is a list of vectors:
  #
  #   a: {V.alpha1  | V.alpha2  | V.alphaZ1 | V.alphaZ2 | V.alphaR}
  #   b: {V.beta1   | V.beta2   | V.betaZ1}
  #   c: {V.gamma1  | V.gamma2  | V.gammaZ1 | V.gammaZ2 | V.gammaR | V.gammaT}
  #   d: {b0        | b_Z1      | b_Z2}
  # Tau: {Scale     | Shape     | beta0     | beta_Z1   | beta_Z2}
  # --------------------------------------------------
  # Generation of Z
  z0 <- 1
  u <- runif(n, 0, 1)
  z1 <- runif(n, 0, 1)
  z2 <- (z1 < 0.5 & u > 2/3) + (z1 > 0.5 & u < 2/3)
  Z <- cbind(z0,z1,z2)
  
  # Generation of D
  d <- rbinom(n, 1, plogis(Z %*% par$d))
  
  # Generation of T
  scale.T <- par$scale
  shape.T <- par$shape
  scale.tau <- exp(-(Z %*% c(0,par$Tau))/shape.T)/scale.T
  tau <- ifelse(d==1, rweibull(n, shape = shape.T, scale = scale.tau), Inf)
  
  # Generation of R
  u <- runif(n, 0, 1)
  gen.age <- function(u) ifelse(u <= 0.8, 
                                1.25 * u + 0.5, 
                                uniroot(function(r, u) r^2 - 4 * r + 2.75 + 1.25 * u, c(1.5, 2), u = u)$root)
  rz1 <- sapply(u, gen.age) # age distribution for Z2=1
  rz0 <- runif(n, 0.3, 2) # age distribution for Z2=0
  r <- ifelse(z2, rz1, rz0)
  
  # Generation of delta
  delta <- (tau<r)*(d==1)+0
  
  simulated_data <- data.frame(Z1 = z1, Z2 = z2, D = d, Tau = pmin(tau, r), R = r, delta = delta)
  
  v.probs <- estimate.outcome.probabilities(par=par, data=simulated_data, delta="delta", k=k, outcome.model = outcome.model)
  
  u <- runif(n, 0, 1)
  v.cummulative.probs <- lapply(v.probs, function(df){t(apply(df, 1, cumsum))})
  v.cummulative.by.D <- v.cummulative.probs$D0
  v.cummulative.by.D[d==1,] <- v.cummulative.probs$D1[d==1,]  # If d is 1, take rows from D1 else D0
  
  v <- sapply(1:n, function(i) {
    which(v.cummulative.by.D[i,] > u[i])[1]  # Find first column index where v[i,] > u[i]
  })
  
  simulated_data <- cbind(simulated_data, V=ordered(v))
}
