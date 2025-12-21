probs.po <- function(dat, alpha0, beta0, gamma0, alphaZ, gammaZ, alphaR, gammaR, gammaT, gammaTR, rangeZ, rangeR, rangeT, rangeTR){
  
  if((length(alpha0) != length(beta0)) | (length(alpha0) != length(gamma0))){
    print("Error in probs function - length of V is determined by the length of intercepts but they have different lengths")
    return()
  }

  #if((alpha0[-length(alpha0)] > alpha0[-1]) | (beta0[-length(beta0)] > beta0[-1]) | (gamma0[-length(gamma0)] > gamma0[-1])){
  #  print("Proportional odds error - intercepts must hold theta_1<...<theta_k")
  #  return()
  #}
  if(length(alphaZ) != length(gammaZ)){
    print("Error in probs function - lengths of Z coefficient is not equal")
    return()
  }
  
  dat <- as.matrix(dat)
  k <- length(alpha0)
  d <- nrow(dat)
  
  # creating the inner products
  xi.alpha <- as.numeric(dat[,rangeZ] %*% alphaZ) + (dat[,rangeR] * alphaR)
  xi.equal.effect.beta.gamma <- as.numeric(dat[,rangeZ] %*% gammaZ) + (dat[,rangeR] * gammaR)
  xi.gamma <- xi.equal.effect.beta.gamma + (dat[,rangeT] * gammaT) + (dat[,rangeTR] * gammaTR)
  
  # for parameters alpha: 
  # => D = 0, delta = 0
  base.g.D0 <- 1/(1+exp(matrix(alpha0, nrow = d, ncol = k, byrow = TRUE) + matrix(xi.alpha, nrow = d, ncol = k, byrow = FALSE)))
  probs.0 <- cbind(1,base.g.D0)-cbind(base.g.D0,0)

  # for parameters beta/gamma:
  # => D = 1, delta = 1
  base.g.D1 <- 1/(1+exp(matrix(gamma0, nrow = d, ncol = k, byrow = TRUE) + matrix(xi.gamma, nrow = d, ncol = k, byrow = FALSE)))
  
  # => D = 1, delta = 0
  condition1 <- (dat[,rangeT]==dat[,rangeR])
  base.g.D1[condition1,] <- 1/(1+exp(matrix(beta0, nrow = d, ncol = k, byrow = TRUE) + matrix(xi.equal.effect.beta.gamma, nrow = d, ncol = k, byrow = FALSE)))[condition1,]
  probs.1 <- cbind(1,base.g.D1)-cbind(base.g.D1,0)
  return(list("D0" = probs.0, "D1" = probs.1))
}

probs.acat <- function(dat, alpha0, beta0, gamma0, alphaZ, gammaZ, alphaR, gammaR, gammaT, gammaTR, rangeZ, rangeR, rangeT, rangeTR){
  if((length(alpha0) != length(beta0)) | (length(alpha0) != length(gamma0))){
    print("Error in probs function - length of V is determined by the length of intercepts but they have different lengths")
    return()
  }
  if(length(alphaZ) != length(gammaZ)){
    print("Error in probs function - lengths of Z coefficient is not equal")
    return()
  }
  
  dat <- as.matrix(dat)
  k <- length(alpha0)
  d <- nrow(dat)

  # creating the inner products
  xi.alpha <- as.numeric(dat[,rangeZ] %*% alphaZ) + (dat[,rangeR] * alphaR)
  xi.equal.effect.beta.gamma <- as.numeric(dat[,rangeZ] %*% gammaZ) + (dat[,rangeR] * gammaR)
  xi.gamma <- xi.equal.effect.beta.gamma + (dat[,rangeT] * gammaT) + (dat[,rangeTR] * gammaTR)
  
  # for parameters alpha: 
  # => D = 0, delta = 0
  g.D0 <- apply(cbind(1,exp(matrix(alpha0, nrow = d, ncol = k, byrow = TRUE) + matrix(xi.alpha, nrow = d, ncol = k, byrow = FALSE))),1,cumprod)
  base.g.D0 <- matrix(rep(colSums(g.D0), each = k+1),nrow = k+1, byrow = FALSE)
  probs.0 <- t(g.D0/base.g.D0)
  
  # for parameters beta/gamma:
  # => D = 1, delta = 1
  g.D1 <- apply(cbind(1,exp(matrix(gamma0, nrow = d, ncol = k, byrow = TRUE) + matrix(xi.gamma, nrow = d, ncol = k, byrow = FALSE))),1,cumprod)
  
  # => D = 1, delta = 0
  condition1 <- (dat[,rangeT]==dat[,rangeR])
  g.D1[,condition1] <- apply(cbind(1,exp(matrix(beta0, nrow = d, ncol = k, byrow = TRUE) + matrix(xi.equal.effect.beta.gamma, nrow = d, ncol = k, byrow = FALSE))),1,cumprod)[,condition1]
  
  base.g.D1 <- matrix(rep(colSums(g.D1), each = k+1),nrow = k+1, byrow = FALSE)
  probs.1 <- t(g.D1/base.g.D1)
  
  print(paste0("TIME - Probs: ",Sys.time()-start.time.p))
  
  return(list("D0" = probs.0, "D1" = probs.1))
}

estimate.outcome.probabilities <- function(par.list, data, delta, k, outcome.model=c("PO","ACAT")){
  
  # --------------------------------------------------
  # what if par is a list of vectors:
  #
  # a: {V.alpha1  | V.alpha2  | V.alphaZ1 | V.alphaZ2 | V.alphaR}
  # b: {V.beta1   | V.beta2   | V.betaZ1}
  # c: {V.gamma1  | V.gamma2  | V.gammaZ1 | V.gammaZ2 | V.gammaR | V.gammaT}
  # --------------------------------------------------
  
  #par.list$b <- c(par.list$b, par.list$c[!sub("V.gamma", "", names(par.list$c)) %in% sub("V.beta", "", names(par.list$b))])
  
  # print(par.list)
  a.cvars <- NULL
  e.cvars <- NULL
  b.cvars <- NULL
  c.cvars <- NULL
  mm.a <- rep(1,length(par.list$a))
  mm.e <- rep(1,length(par.list$e)+1)
  mm.b <- rep(1,length(par.list$b))
  mm.c <- rep(1,length(par.list$c))
  
  
  a.inter <- par.list$a[1:(k-1)]
  if(length(par.list$a)>=k) a.cvars <- par.list$a[k:length(par.list$a)]
  if(length(par.list$e)>0) e.cvars <- par.list$e
  b.inter <- par.list$b[1:(k-1)]
  if(length(par.list$b)>=k) b.cvars <- par.list$b[k:length(par.list$b)]
  c.inter <- par.list$c[1:(k-1)]
  if(length(par.list$c)>=k) c.cvars <- par.list$c[k:length(par.list$c)]
  
  
  if(!is.null(a.cvars)) mm.a <- model.matrix(as.formula(paste("~",paste(sub("V.alpha.","",names(a.cvars)),collapse="+"))), data=data)
  if(!is.null(e.cvars)) mm.e <- model.matrix(as.formula(paste("~ -1 +",paste(sub("V.eta.","",names(e.cvars)),collapse="+"))), data=data)
  if(!is.null(b.cvars)) mm.b <- model.matrix(as.formula(paste("~",paste(sub("V.beta.","",names(b.cvars)),collapse="+"))), data=data)
  if(!is.null(c.cvars)) mm.c <- model.matrix(as.formula(paste("~",paste(sub("V.gamma.","",names(c.cvars)),collapse="+"))), data=data)
  
  if(!is.null(a.cvars)) {
    xi.a <- mm.a %*% t(cbind(a.inter, matrix(a.cvars, nrow = k-1, ncol = length(a.cvars), byrow=TRUE)))
  } else{
    xi.a <- t(mm.a %*% t(matrix(a.inter, nrow = nrow(data), ncol = k - 1, byrow = TRUE)))
  }
  if(!is.null(b.cvars)) {
    xi.b <- mm.b %*% t(cbind(b.inter, matrix(b.cvars, nrow = k-1, ncol = length(b.cvars), byrow=TRUE))) 
  } else{
    xi.b <- matrix(b.inter, nrow = nrow(data), ncol = k - 1, byrow = TRUE)
  }
  if(!is.null(c.cvars)) {
    xi.c <- mm.c %*% t(cbind(c.inter, matrix(c.cvars, nrow = k-1, ncol = length(c.cvars), byrow=TRUE)))
  } else{
    xi.c <- matrix(c.inter, nrow = nrow(data), ncol = k - 1, byrow = TRUE)
  }
  
  if(!is.null(e.cvars)) {
    xi.b <- xi.b + matrix(mm.e %*% e.cvars, nrow = nrow(data), ncol = k - 1, byrow = FALSE)
    xi.c <- xi.c + matrix(mm.e %*% e.cvars, nrow = nrow(data), ncol = k - 1, byrow = FALSE)
  }
  
  condition1 <- (data[,delta]==0)
  
  if(outcome.model=="PO"){
    
    xi.a <- exp(-xi.a)
    xi.b <- exp(-xi.b)
    xi.c <- exp(-xi.c)

    # for parameters alpha: 
    # => D = 0, delta = 0
    base.g.D0 <- 1/(1+xi.a)
    probs.0 <- cbind(base.g.D0,1)-cbind(0,base.g.D0)
    # for parameters beta/gamma:
    # => D = 1, delta = 1
    base.g.D1 <- 1/(1+xi.c)
    
    # => D = 1, delta = 0
    base.g.D1[condition1,] <- 1/(1+xi.b)[condition1,]
    probs.1 <- cbind(base.g.D1,1)-cbind(0,base.g.D1)
  }
  if(outcome.model=="ACAT"){
    
    xi.a <- cbind(1,exp(xi.a))
    xi.b <- cbind(1,exp(xi.b))
    xi.c <- cbind(1,exp(xi.c))
    
    ## Initialize with 1's
    g.D0 <- c.xi.b <- g.D1 <- matrix(1, nrow = nrow(data), ncol = k)
    
    ## Vectorized adjacent-column products for columns 2:k
    idx <- 2:k
    g.D0[, idx]   <- xi.a[, idx] * xi.a[, idx - 1]
    c.xi.b[, idx] <- xi.b[, idx] * xi.b[, idx - 1]
    g.D1[, idx]   <- xi.c[, idx] * xi.c[, idx - 1]
    
    
    # for parameters alpha: 
    probs.0 <- g.D0/rowSums(g.D0)

    # for parameters beta/gamma:
    g.D1[condition1,] <- c.xi.b[condition1,]
    probs.1 <- g.D1/rowSums(g.D1)
  }
  
  return(list("D0" = unname(probs.0), "D1" = unname(probs.1)))
}
