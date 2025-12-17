naive.analysis <- function(dat, SE=TRUE, V.model="PO"){
  # Naive analysis is not dynamic - I need to generalize the ranges
  
  # Check input validity
  if (!"delta" %in% colnames(dat)) {
    warning("The input data frame must contain a column named 'delta'.")
  }
  if (!all(c("V", "Z1", "Z2", "R") %in% colnames(dat))) {
    warning("The input data frame must contain columns 'V', 'Z1', 'Z2', and 'R'.")
  }
  if (!is.logical(SE)) {
    warning("The argument 'SE' must be a logical value (TRUE or FALSE).")
  }
  if (!V.model %in% c("PO", "ACAT")) {
    warning("The argument 'V.model' must be either 'PO' or 'ACAT'.")
  }
  
  # Subset data
  delta0 <- subset(dat, delta == 0, select=c("Z1","Z2","R","delta","V"))
  if (nrow(delta0) == 0) {
    warning("No rows in the data satisfy the condition delta == 0.")
  }
  
  # Helper function to fit the model
  fit_model <- function(model, formula, data){
    fit <- tryCatch({
      if (model == "PO") polr(formula, data = data, Hess = TRUE)
      else if (model == "ACAT") vglm(formula, acat(reverse = FALSE, parallel = TRUE), data = data)
    },
    error = function(e){
      warning(paste("Naive model fitting failed for", model, "error:", e$message))
    })
    return(fit)
  }
  
  # Helper function to extract estimates
  extract_estimates <- function(fit, model, data){
    if (model == "PO") est.val <- unname(c(fit$zeta, -fit$coefficients))
    else if (model == "ACAT") est.val <- unname(coef(fit))
    if(length(est.val) != nlevels(data$V)-1 + ncol(data)-2) warning(paste("Naive model fitting failed for", model, "error: Too much / Missing coefficient estimates."))
    names(est.val) <- c("V.alpha1", "V.alpha2", "V.alphaZ1", "V.alphaZ2", "V.alphaR")
    return(est.val)
  }
  
  # Helper function to compute standard errors
  compute_standard_errors <- function(fit, model, data){
    if (model == "PO") {
      std <- sqrt(diag(vcov(fit)))
      est.se <- unname(c(tail(std, 2), head(std, length(std) - 2)))
    }
    else if (model == "ACAT"){
      est.se <- unname(sqrt(diag(vcov(fit))))
    }
    if(length(est.se) != nlevels(data$V)-1 + ncol(data)-2) warning(paste("Naive model fitting failed for", model, "error: Too much / Missing coefficient estimates"))
        names(est.se) <- c("V.alpha1", "V.alpha2", "V.alphaZ1", "V.alphaZ2", "V.alphaR")
    return(est.se)
  }
  
  # Fit the model
  fit <- fit_model(V.model, V ~ Z1 + Z2 + R, delta0)
  # Extract estimates
  est.val <- extract_estimates(fit, V.model, delta0)
  # Extract SE
  if (SE) {
    est.se <- compute_standard_errors(fit, V.model, delta0)
    return(list(est = est.val, est.se = est.se))
  }
  
  return(est.val)
}

complete.case.analysis <- function(dat, SE=TRUE, V.model="PO"){
  # Complete Case analysis is not dynamic - I need to generalize the ranges
  
  # Check input validity
  if (!"delta" %in% colnames(dat)) {
    warning("The input data frame must contain a column named 'delta'.")
  }
  if (!all(c("V", "Z1", "Z2", "R", "Tau") %in% colnames(dat))) {
    warning("The input data frame must contain columns 'V', 'Z1', 'Z2', 'R' and 'T'.")
  }
  if (!is.logical(SE)) {
    warning("The argument 'SE' must be a logical value (TRUE or FALSE).")
  }
  if (!V.model %in% c("PO", "ACAT")) {
    warning("The argument 'V.model' must be either 'PO' or 'ACAT'.")
  }
  
  # Subset data
  delta1 <- subset(dat, delta == 1, select = c("Z1","Z2","R","Tau","TR","delta","V"))
  if (nrow(delta1) == 0) {
    warning("No rows in the data satisfy the condition delta == 1.")
  }
  
  # Helper function to fit the model
  fit_model <- function(model, formula, data){
    fit <- tryCatch({
      if (model == "PO") polr(formula, data = data, Hess = TRUE)
      else if (model == "ACAT") vglm(formula, acat(reverse = FALSE, parallel = TRUE), data = data)
    },
    error = function(e){
      warning(paste("Complete case model fitting failed for", model, "error:", e$message))
    })
    return(fit)
  }
  
  # Helper function to extract estimates
  extract_estimates <- function(fit, model, data){
    if (model == "PO") est.val <- unname(c(fit$zeta, -fit$coefficients))
    else if (model == "ACAT") est.val <- unname(coef(fit))
    if(length(est.val) != nlevels(data$V)-1 + ncol(data)-2) warning(paste("Complete case model fitting failed for", model, "error: Too much / Missing coefficient estimates."))
    names(est.val) <- c("V.gamma1", "V.gamma2", "V.gammaZ1", "V.gammaZ2", "V.gammaR", "V.gammaT", "V.gammaTR")
    return(est.val)
  }
  
  # Helper function to compute standard errors
  compute_standard_errors <- function(fit, model, data){
    if (model == "PO") {
      std <- sqrt(diag(vcov(fit)))
      est.se <- unname(c(tail(std, 2), head(std, length(std) - 2)))
    }
    else if (model == "ACAT"){
      est.se <- unname(sqrt(diag(vcov(fit))))
    }
    if(length(est.se) != nlevels(data$V)-1 + ncol(data)-2) warning(paste("Complete case model fitting failed for", model, "error: Too much / Missing coefficient estimates"))
    names(est.se) <- c("V.gamma1", "V.gamma2", "V.gammaZ1", "V.gammaZ2", "V.gammaR", "V.gammaT", "V.gammaTR")
    return(est.se)
  }
  
  # Fit the model
  fit <- fit_model(V.model, V ~ Z1 + Z2 + R + Tau + TR, delta1)
  # Extract estimates
  est.val <- extract_estimates(fit, V.model, delta1)
  # Extract SE
  if (SE) {
    est.se <- compute_standard_errors(fit, V.model, delta1)
    return(list(est = est.val, est.se = est.se))
  }
  
  return(est.val)
}

proposed.method.analysis <- function(dat, rangeZ, rangeT, rangeR, rangeTR, rangeDelta, rangeV, V.model="PO"){
  
  # Check input validity
  if (!"delta" %in% colnames(dat)) {
    warning("The input data frame must contain a column named 'delta'.")
  }
  if (!all(c("V", "Z1", "Z2", "R", "Tau") %in% colnames(dat))) {
    warning("The input data frame must contain columns 'V', 'Z1', 'Z2', 'R' and 'T'.")
  }
  if (!V.model %in% c("PO", "ACAT")) {
    warning("The argument 'V.model' must be either 'PO' or 'ACAT'.")
  }
  
  # Ensure input data is valid
  delta0 <- subset(dat, delta == 0)
  if (nrow(delta0) == 0) {
    warning("No rows with delta == 0 in the dataset.")
  }
  
  # fit a PH mixture cure model:
  # formula = Surv(Tau,delta)~All Zs, cureform = ~All Zs
  cd <- colnames(dat)
  c.cox <- paste0("~", paste(cd[rangeZ], collapse="+"))
  f.cox <- paste0("Surv(",cd[rangeT], ",", cd[rangeDelta],")", c.cox)
  
  #DO NOT PASS THE INTERACTION TR, for some strange reason it changes the estimation
  # Fit the PH mixture cure model (excluding problematic columns)
  sink("NUL")
  fit.cox <- tryCatch(
    smcure(
      formula = as.formula(f.cox), 
      cureform = as.formula(c.cox), 
      data=dat[c(-3,-6)], # GENERALIZE ME!!!
      model="ph", 
      Var = FALSE),
    error = function(e){
      warning(paste("smcure failure:", e$message))
      })
  sink()
  
  compute_cox_weights <- function(fit, z, delta) {
    # Compute weights for the PH mixture cure model
    Zs <- cbind(1, z)
    w0 <- 1 / (1 + exp(fit$b %*% t(Zs)))
    w1 <- 1 - w0
    ebetaX <- exp(fit$beta %*% t(z))
    s1 <- fit.cox$s[delta == 0]^ebetaX
    
    weights <- data.frame(
      w0 = rep(NA, length(delta)),
      w1 = rep(NA, length(delta))
    )
    weights$w0[delta == 0] <- w0
    weights$w1[delta == 0] <- w1 * s1
    return(weights)
  }
  
  beta.est <- c(fit.cox$b, fit.cox$beta)
  names(beta.est) <- c("D:b0", "D:bZ1", "D:bZ2", "T:bZ1", "T:bZ2")
  
  print("Cure model estimates:")
  print(beta.est)
  
  cox.weights <- compute_cox_weights(fit = fit.cox, z = delta0[, rangeZ], delta = dat[, rangeDelta])
  
  naive.init.values <- naive.analysis(dat, SE=FALSE, V.model=V.model)
  complete.case.init.values <- complete.case.analysis(dat, SE=FALSE, V.model=V.model)
  
  # Full params:
  # V.alpha0, V.alphaZ, V.alphaR
  # V.beta0, V.gammaZ, 
  # V.gamma0, V.gammaR, V.gammaT, V.gammaTR,
  #
  # in our case we have 3 outcomes and 2 Zs (include column of 1) so we have 14 parameters
  
  par.init1 <- c(
   naive.init.values, # alpha initial params
   complete.case.init.values["V.gamma1"], complete.case.init.values["V.gamma2"], # beta initial params = gamma initial params
   complete.case.init.values) # Gamma initial params
  
  #par.init2 <- c(
  #  naive.init.values, # alpha initial params
  #  naive.init.values["V.alpha1"], naive.init.values["V.alpha2"], # beta initial params = alpha initial params
  #  complete.case.init.values) # Gamma initial params
  
  names(par.init1) <- c(#names(par.init2) <- c(
    "V.alpha1", "V.alpha2", "V.alphaZ1", "V.alphaZ2", "V.alphaR",
    "V.beta1", "V.beta2", "V.gamma1", "V.gamma2",
    "V.gammaZ1", "V.gammaZ2", "V.gammaR", "V.gammaT", "V.gammaTR")
  
  lower.b <- c(-Inf, -Inf, -Inf, -Inf, -Inf,
               -Inf, -Inf, -Inf, -Inf,
               -Inf, -Inf, -Inf, -Inf, -Inf)
  
  if(V.model=="PO"){
    lower.b <- c(-Inf, 1e-8, -Inf, -Inf, -Inf,
                 -Inf, 1e-8, -Inf, -Inf,
                 -Inf, 1e-8, -Inf, -Inf, -Inf)
    par.init1["V.alpha2"] <- par.init1["V.alpha2"]-par.init1["V.alpha1"]
    par.init1["V.beta2"] <- par.init1["V.beta2"]-par.init1["V.beta1"]
    par.init1["V.gamma2"] <- par.init1["V.gamma2"]-par.init1["V.gamma1"]
  }
  sol1 <- tryCatch(
    optim(par=par.init1, fn=pseudo.neg.log.like,
          cox.weights=cox.weights, rangeZ=rangeZ, rangeT=rangeT, rangeR=rangeR,
          rangeTR=rangeTR, rangeV=rangeV, V.model=V.model, dat=dat, method="L-BFGS-B",
          lower=lower.b, control=list(trace = 1, maxit=2000)),
    error = function(e){
      warning(paste("optim failure:", e$message))
    })
  
  if(V.model=="PO"){
    sol1$par["V.alpha2"] <- sol1$par["V.alpha2"]+sol1$par["V.alpha1"]
    sol1$par["V.beta2"] <- sol1$par["V.beta2"]+sol1$par["V.beta1"]
    sol1$par["V.gamma2"] <- sol1$par["V.gamma2"]+sol1$par["V.gamma1"]
  }
  #print(sol1)
  #sol2 <- try(optim(par=par.init2, fn=neg.loglik.est.W, d=d, method="L-BFGS-B", lower=lower.b, control=list(trace=0, maxit=2000)))
  
  #if(sol1$value < sol2$value)
  return(list(cure=beta.est ,est=sol1$par, negloglik=sol1$value, solution="Solution 1"))
  #return(list(surv=est, CC=complete.case, Naive=naive.case, V=sol2$par, negloglik=sol2$value, solution="Solution 2"))
}

ordcuredfit <- function(formula.c, formula.a, equal.effect=NULL, Tau, R, delta, data, outcome.model=c("PO","ACAT")) {
  # PH mixture cure model weights and parameters function
  cox.PH.weights <- function(fit.cox, mfd0, Z, delta, data){
    # Compute weights for the PH mixture cure model
    w0 <- 1 / (1 + exp(fit.cox$b %*% t(cbind(1, mfd0[,Z]))))
    ebetaX <- exp(fit.cox$beta %*% t(mfd0[,Z]))
    s1 <- fit.cox$s[data[,delta] == 0]^ebetaX
    
    weights <- data.frame(w0 = rep(NA, n), w1 = rep(NA, n))
    weights$w0[data[,delta] == 0] <- w0
    weights$w1[data[,delta] == 0] <- (1 - w0) * s1
    
    beta.est <- c(fit.cox$b, fit.cox$beta)
    names(beta.est) <- c("D:b0", paste0(rep(c("D:b", "T:b"), each = length(Z)), Z))
    
    return(list(beta.est = beta.est, weights= weights))
  }
  # Naive and Complete cases estimation function
  n.c.estimation <- function(f.c, f.a, mfd0, mfd1, k, outcome.model){
    if (outcome.model == "PO") {
      naive.fit <- polr(f.a, data = mfd0, Hess = TRUE)
      cc.fit <- polr(f.c, data = mfd1, Hess = TRUE)
      
      naive.val <- unname(c(naive.fit$zeta, -naive.fit$coefficients))
      cc.val <- unname(c(cc.fit$zeta, -cc.fit$coefficients))
      
      # Compute standard errors using vcov()
      std <- sqrt(diag(vcov(naive.fit)))
      naive.se <- unname(c(tail(std, k-1), head(std, length(std) - k + 1)))
      
      std <- sqrt(diag(vcov(cc.fit)))
      cc.se <- unname(c(tail(std, k-1), head(std, length(std) - k + 1)))
    } 
    else if (outcome.model == "ACAT") {
      naive.fit <- vglm(f.a, acat(reverse = FALSE, parallel = TRUE), data = mfd0)
      cc.fit <- vglm(f.c, acat(reverse = FALSE, parallel = TRUE), data = mfd1)
      
      naive.val <- unname(coef(naive.fit))
      cc.val <- unname(coef(cc.fit))
      
      naive.se <- unname(sqrt(diag(vcov(naive.fit))))
      cc.se <- unname(sqrt(diag(vcov(cc.fit))))
    }
    
    naive.vars <- attr(terms(f.a), "term.labels")
    cc.vars <- attr(terms(f.c), "term.labels")
    
    names(naive.val) <- names(naive.se) <- c(paste0("V.alpha", 1:(k-1)), paste0("V.alpha", naive.vars))
    names(cc.val) <- names(cc.se) <- c(paste0("V.gamma", 1:(k-1)), paste0("V.gamma", cc.vars))
    
    summary.n <- cbind(naive.val, naive.se)
    summary.c <- cbind(cc.val, cc.se)
    colnames(summary.n) <- colnames(summary.c) <- c("Estimate","Std. Error")
    
    return(list(naive = summary.n, cc = summary.c))
  }
  # Initialization values for the full model
  par.initialization <- function(equal.effect, naive.est, cc.est, k, init.case=c("a","c")){
    alpha.vars <- sub("V.alpha","",names(naive.est))
    gamma.vars <- sub("V.gamma","",names(cc.est))
    beta.vars <- setdiff(alpha.vars, all.vars(equal.effect))
    
    if(init.case == "a")  beta.params <- naive.est[paste0("V.alpha", beta.vars)] # beta initial params = alpha initial params
    if(init.case == "c")  beta.params <- cc.est[paste0("V.gamma", beta.vars)] # beta initial params = gamma initial params
    
    par <- c(naive.est, beta.params, cc.est)
    names(par) <- c(names(naive.est), paste0("V.beta", beta.vars), names(cc.est))
    lower.b <- rep(-Inf,length(par))
    lengths <- c(length(alpha.vars), length(beta.vars), length(gamma.vars))
    
    if(outcome.model == "PO" & k > 2){
      lower.b[2:(k-1)] <- lower.b[lengths[1] + (2:(k-1))] <- lower.b[lengths[1] + lengths[2] + (2:(k-1))] <- 1e-8
      par[2:(k-1)] <- par[2:(k-1)]-par[(2:(k-1)) - 1]
      par[lengths[1] + (2:(k-1))] <- par[lengths[1] + (2:(k-1))]-par[lengths[1] + (2:(k-1)) - 1]
      par[lengths[1] + lengths[2] + (2:(k-1))] <- par[lengths[1] + lengths[2] + (2:(k-1))]-par[lengths[1] + lengths[2] + (2:(k-1)) - 1]
    }
    return(list(lengths = lengths, lower.b = lower.b, par.init = par))
  }
  
  # formula.c: V ~ Z1 + Z2 + R + T + T*R
  # equal.effect: ~ Z1 + Z2 + R # need to check if all in 'formula.c'
  
  # V ~ Z1 + Z2 + R + T + T*R ------for gamma group
  # V ~ 1
  # V ~ Z1 + Z2 + R ------for alpha group
  call <- match.call() 
  outcome.model <- match.arg(outcome.model)
  n <- dim(data)[1]
  # --- TO STOP IF: n>7 why cox.ph doesn't converge (it doesn't matter because we change to Weibull)
  if(is.null(formula.a) || is.null(formula.c) 
     || identical(as.character(formula.a), "") || identical(as.character(formula.c), "")){stop("one or both of the models is empty")}
  
  if(!identical(model.extract(model.frame(formula.a, data), "response"),
                model.extract(model.frame(formula.c, data), "response"))){stop("response must be the same for both models")}
  
  if(delta %in% c(all.vars(formula.a), all.vars(formula.c))){stop("delta cannot be in the models")}
  
  if(Tau %in% all.vars(formula.a)){stop("Tau cannot be in model.a, it is not defined for the uncured fraction")}
  
  # Complete case model
  f.c <- update(formula.c, paste("~ . +",Tau,"+",R))
  # Naive case model
  f.a <- update(formula.a, paste("~ . +",R))
  
  if(!all(all.vars(equal.effect) %in% all.vars(f.a))){stop("variables in equal.effect are not in uncured fraction model")}
  
  if(any(vapply(model.frame(f.a, data), function(col) all(col == col[1]), logical(1))) ||
     any(vapply(model.frame(f.c, data), function(col) all(col == col[1]), logical(1)))){stop("one of the variables is constant")}
  
  cvars <- all.vars(f.c)
  mf <- model.frame(f.c, data)
  mfd0 <- mf[data[,delta]==0, !names(mf) %in% Tau]
  mfd1 <- mf[data[,delta]==1,]
  Y <- model.extract(mf, "response")
  if(is.null(Y)){stop("no response variable found")}
  k <- length(levels(Y))
  if(k<3){stop("response must have 3 or more levels")} # I might have to change polr so k=2 is valid too
  
  # call for cox PH cure model weights estimation function
  cox.est <- tryCatch({
    sink("NUL")
    Z <- cvars[!cvars %in% c(colnames(mf)[1],Tau,R)]
    c.cox <- paste0("~", paste(Z, collapse="+"))
    f.cox <- paste0("Surv(",Tau, ",", delta,")", c.cox)
    fit.cox <- smcure(formula = as.formula(f.cox), cureform = as.formula(c.cox), data = data, model = "ph", Var = FALSE)
    sink()
    cox.PH.weights(fit.cox, mfd0, Z, delta, data)},
    error = function(e){
      warning(paste("smcure failure:", e$message))
    })
  
  # call for naive and complete estimation function
  n.c.est <- n.c.estimation(f.c, f.a, mfd0, mfd1, k, outcome.model)
  
  beta.as.alpha.flag <- FALSE
  
  solution <- tryCatch({
    par.init <- par.initialization(equal.effect, n.c.est$naive[,"Estimate"], n.c.est$cc[,"Estimate"], k, init.case = "c")
    l <- par.init$lengths
    
    optim_result <- NULL
    optim_trace <- capture.output({
      optim_result <- optim(par = par.init$par, fn = pseudo.neg.log.like,
                            weights = cox.est$weights, k = k, lengths = l, delta = delta, V = names(mf)[1],
                            outcome.model = outcome.model, data = data, method = "L-BFGS-B",
                            lower = par.init$lower.b, control = list(trace = 1, maxit = 2000))})
    list(optim = optim_result, trace = optim_trace)}, error = function(e) {
    warning(paste("First optim failure:", e$message, "Retrying with beta as alpha"))
    tryCatch({
      beta.as.alpha.flag <- TRUE
      par.init <- par.initialization(equal.effect, n.c.est$naive[,"Estimate"], n.c.est$cc[,"Estimate"], k, init.case = "a")
      l <- par.init$lengths
      
      optim_result <- NULL
      optim_trace <- capture.output({
        optim_result <- optim(par = par.init$par, fn = pseudo.neg.log.like,
                              weights = cox.est$weights, k = k, lengths = l, delta = delta, V = names(mf)[1],
                              outcome.model = outcome.model, data = data, method = "L-BFGS-B",
                              lower = par.init$lower.b, control = list(trace = 1, maxit = 2000))})
      list(optim = optim_result, trace = optim_trace)}, error = function(e2) {
      warning(paste("Second optim failure:", e2$message))
      return(NULL)
    })
  })
  
  if(outcome.model == "PO" & k > 2){
    solution$optim$par[2:(k-1)] <- solution$optim$par[2:(k-1)]+solution$optim$par[(2:(k-1)) - 1]
    solution$optim$par[l[1] + (2:(k-1))] <- solution$optim$par[l[1] + (2:(k-1))]+solution$optim$par[l[1] + (2:(k-1)) - 1]
    solution$optim$par[l[1] + l[2] + (2:(k-1))] <- solution$optim$par[l[1] + l[2] + (2:(k-1))]+solution$optim$par[l[1] + l[2] + (2:(k-1)) - 1]
  }
  
  par.list <- list(d = cox.est$beta.est[1:3],
                   Tau = cox.est$beta.est[4:5],
                   a = solution$optim$par[1:l[1]],
                   b = solution$optim$par[(l[1]+1):(l[1]+l[2])],
                   c = solution$optim$par[(l[1]+l[2]+1):(l[1]+l[2]+l[3])])
  
  return(list(par = par.list, 
              value = solution$optim$value, 
              beta.as.alpha.flag = beta.as.alpha.flag, 
              naive = n.c.est$naive, 
              cc = n.c.est$cc,
              trace = solution$trace))
}

ordcuredfit.param <- function(survform, cureform, formula.a, formula.e, formula.b, formula.c, equal.effect=NULL, Tau, R, delta, data, outcome.model=c("PO","ACAT"), var = var) {
  
  fit.param.Wei.full.1 <- function(data){
    
    neg.loglik.1 <- function(par, Z, w){
      Z.b <- cbind(1,Z) %*% par
      return(-sum(w * Z.b - log1p(exp(Z.b))))
    }
    
    neg.loglik.2 <- function(par, Z, Tau, delta, w){
      beta <- par[1:(length(par)-2)]
      Z.beta <- Z %*% beta
      a <- par[length(par)-1]
      log.a <- log(a)
      k <- par[length(par)]
      alogk <- a*log(k)
      log.T <- log(Tau)
      return(-sum(w * (delta * (Z.beta + log.a + (a-1)*log.T - alogk) - exp(Z.beta + a*log.T - alogk))))
    }
    
    n <- nrow(data)
    d.CC <- data[data[delta]==1,]
    new_survform <- as.formula(paste0("Surv(time = ",Tau,", event = rep(1, nrow(d.CC))) ~ ",
                                      paste(all.vars(survform)[-1], collapse = " + ")))
    
    fit.param <- try(phreg(new_survform, data = d.CC))
    #try(phreg(Surv(time=Tau, event=rep(1, nrow(d.CC))) ~ Z1+Z2, data=d.CC))
    fit.logist <- glm(cureform, family=binomial, data=data)
    #fit.logist <- glm(delta~Z1+Z2, family=binomial, data=data)
    
    # initial value:
    b.m <- fit.logist$coef
    beta.m <- fit.param$coef[1:(length(fit.param$coef)-2)]
    a.m <- exp(fit.param$coef[length(fit.param$coef)])
    k.m <- exp(-fit.param$coef[length(fit.param$coef)-1])
    
    sol <- c(b.m, beta.m, a.m, k.m)
    # parnames.wei <- c("D:b0", "D:bZ1", "D:bZ2",  "T:bZ1", "T:bZ2", "T:sh", "T:sc")
    
    lower.b <- c(rep(-Inf, length(beta.m)), rep(1e-6, 2)) # only for the survival part
    diff <- Inf
    diff.loglik <- curr.loglik <- 1e8
    tol <- 1e-8
    iter <- 0
    repeat{
      w <- rep(1, n)
      # compute pi and s once
      pi <- exp(model.matrix(cureform, data) %*% sol[1:length(b.m)])
      pi <- pi / (1 + pi)
      eta <- exp(model.matrix(update(survform, . ~ . - 1), data) %*% sol[length(b.m)+(1:length(beta.m))])
      s <- exp(-((data[[Tau]] / sol[length(sol)])^sol[length(sol)-1]) * eta)
      # update only relevant entries
      i0 <- data[[delta]] == 0
      w[i0] <- (pi[i0] * s[i0]) / (pi[i0] * s[i0] + (1 - pi[i0]))
      
      sol.1 <- try(optim(par=sol[1:length(b.m)], fn=neg.loglik.1,
                         Z = model.matrix(update(cureform, . ~ . - 1), data = data), w=w, method="BFGS",
                         control=list(trace=0, maxit=2000)))
      
      sol.2 <- try(optim(par=sol[(length(b.m)+1):length(sol)], fn=neg.loglik.2,
                         Z = model.matrix(update(survform, . ~ . - 1), data), Tau = data[[Tau]], delta = data[[delta]], w=w,
                         method="L-BFGS-B", lower=lower.b,
                         control=list(trace=0, maxit=2000)))
      
      # --- Check convergence
      iter <- iter + 1
      upd.loglik <- sol.1$value + sol.2$value
      diff.loglik <- curr.loglik - upd.loglik
      sol.update <- c(sol.1$par, sol.2$par)
      diff <- sqrt(sum((sol-sol.update)^2)) # distance
      sol <- sol.update
      cat("iter =",iter,", loglike =",upd.loglik," ,diff =",diff," ,diff.like =",diff.loglik,"\n")
      curr.loglik <- upd.loglik
      
      if(diff <= tol) break
    }
    
    # compute pi and s once
    pi <- exp(model.matrix(cureform, data) %*% sol[1:length(b.m)])
    pi <- pi / (1 + pi)
    eta <- exp(model.matrix(update(survform, . ~ . - 1), data) %*% sol[length(b.m)+(1:length(beta.m))])
    s <- exp(-((data[[Tau]] / sol[length(sol)])^sol[length(sol)-1]) * eta)
    
    w <- data.frame(w0 = rep(1, n), w1 = rep(1, n))
    w$w0[i0] <- (1-pi)[i0]
    w$w1[i0] <- (pi*s)[i0]

    names(sol) <- c(sub("(Intercept)", "b0", paste0("D:", colnames(model.matrix(cureform, data))), fixed = TRUE), 
                    paste0("T:", colnames(model.matrix(update(survform, . ~ . - 1), data))), "T:sh", "T:sc")

    return(list(est=sol, weights = w))
  }
  # Naive and Complete cases estimation function
  n.c.estimation <- function(f.c, f.a, mfd0, mfd1, k, outcome.model){
    if (outcome.model == "PO" & k>2) {
      naive.fit <- polr(f.a, data = mfd0, Hess = TRUE)
      cc.fit <- polr(f.c, data = mfd1, Hess = TRUE)
      
      naive.val <- unname(c(naive.fit$zeta, -naive.fit$coefficients))
      cc.val <- unname(c(cc.fit$zeta, -cc.fit$coefficients))
      
      # Compute standard errors using vcov()
      std <- rep(0,length(naive.val))
      if(all(diag(vcov(naive.fit))>=0)) {std <- sqrt(diag(vcov(naive.fit)))}
      naive.se <- unname(c(tail(std, k-1), head(std, length(std) - k + 1)))
      
      std <- rep(0,length(cc.val))
      if(all(diag(vcov(cc.fit))>=0)) {std <- sqrt(diag(vcov(cc.fit)))}
      cc.se <- unname(c(tail(std, k-1), head(std, length(std) - k + 1)))
    } 
    else{
      naive.fit <- vglm(f.a, acat(reverse = FALSE, parallel = TRUE), data = mfd0)
      cc.fit <- vglm(f.c, acat(reverse = FALSE, parallel = TRUE), data = mfd1)
      naive.val <- unname(coef(naive.fit))
      cc.val <- unname(coef(cc.fit))
      
      naive.se <- unname(sqrt(diag(vcov(naive.fit))))
      cc.se <- unname(sqrt(diag(vcov(cc.fit))))
    }
    
    naive.vars <- attr(terms(f.a), "term.labels")
    cc.vars <- attr(terms(f.c), "term.labels")
    
    names(naive.val) <- names(naive.se) <- c(paste0("V.alpha.", 1:(k-1)), paste0("V.alpha.", naive.vars))
    names(cc.val) <- names(cc.se) <- c(paste0("V.gamma.", 1:(k-1)), paste0("V.gamma.", cc.vars))
    
    summary.n <- cbind(naive.val, naive.se)
    summary.c <- cbind(cc.val, cc.se)
    colnames(summary.n) <- colnames(summary.c) <- c("Estimate","Std. Error")
    
    return(list(naive = summary.n, cc = summary.c))
  }
  # Initialization values for the full model
  par.initialization <- function(equal.effect, naive.est, cc.est, k){
    alpha.vars <- names(naive.est)
    beta.vars <- sub("V.alpha.","",names(naive.est))
    gamma.vars <- sub("V.gamma.","",names(cc.est))
    equal.vars <- attr(terms(equal.effect), "term.labels")
    
    beta.params <- naive.est[paste0("V.alpha.",beta.vars)]
    names(beta.params) <- paste0("V.beta.",beta.vars)
    
    par <- c(naive.est, beta.params, cc.est)
    lengths <- c(length(alpha.vars), length(beta.vars), length(gamma.vars))
    
    if(length(equal.vars)>0) {
      equal.vars <- equal.vars[which(equal.vars %in% gamma.vars)]
      beta.vars <- beta.vars[which(!beta.vars %in% equal.vars)]
      gamma.vars <- gamma.vars[which(!gamma.vars %in% equal.vars)]
      beta.params <- naive.est[paste0("V.alpha.",beta.vars)]
      names(beta.params) <- paste0("V.beta.",beta.vars)
      equal.params <- cc.est[paste0("V.gamma.", equal.vars)]
      gamma.params <- cc.est[paste0("V.gamma.", gamma.vars)]
      
      par <- c(naive.est, beta.params, gamma.params, equal.params)
      lengths <- c(length(alpha.vars), length(beta.vars), length(gamma.vars), length(equal.vars))
    }
    
    lower.b <- rep(-Inf,length(par))
    # SPESIFIC
    #lower.b[c(4,8,9,13,14,15)] <- -0.5
    if(outcome.model == "PO" & k > 2){
      lower.b[2:(k-1)] <- lower.b[lengths[1] + (2:(k-1))] <- lower.b[lengths[1] + lengths[2] + (2:(k-1))] <- 1e-8
      par[2:(k-1)] <- diff(par[1:(k-1)])
      par[lengths[1] + (2:(k-1))] <- diff(par[lengths[1] + (1:(k-1))])
      par[lengths[1] + lengths[2] + (2:(k-1))] <- diff(par[lengths[1] + lengths[2] + (1:(k-1))])
    }
    return(list(lengths = lengths, lower.b = lower.b, par = par))
  }
  
  # formula.c: V ~ Z1 + Z2 + R + T + T*R
  # equal.effect: ~ Z1 + Z2 + R # need to check if all in 'formula.c'
  
  # V ~ Z1 + Z2 + R + T + T*R ------for gamma group
  # V ~ 1
  # V ~ Z1 + Z2 + R ------for alpha group
  call <- match.call() 
  outcome.model <- match.arg(outcome.model)
  n <- dim(data)[1]
  # --- TO STOP IF: n>7 why cox.ph doesn't converge (it doesn't matter because we change to Weibull)
  if(is.null(formula.a) || is.null(formula.c) 
     || identical(as.character(formula.a), "") || identical(as.character(formula.c), "")){stop("one or both of the models is empty")}
  
  if(!identical(model.extract(model.frame(formula.a, data), "response"),
                model.extract(model.frame(formula.c, data), "response"))){stop("response must be the same for both models")}
  
  if(delta %in% c(all.vars(formula.a), all.vars(formula.c))){stop("delta cannot be in the models")}
  
  if(Tau %in% all.vars(formula.a)){stop("Tau cannot be in model.a, it is not defined for the uncured fraction")}
  
  combine_formulas_clean <- function(..., move_interactions_last = TRUE) {
    formulas <- list(...)
    # Take the LHS from the first formula
    lhs <- deparse(formulas[[1]][[2]])
    
    # Extract all RHS parts
    rhs_all <- unlist(lapply(formulas, function(f) {
      if (length(f) >= 3) deparse(f[[3]]) else "1"
    }))
    
    # Split all RHS strings by '+' and clean up
    rhs_terms <- trimws(unlist(strsplit(paste(rhs_all, collapse = " + "), "\\+")))
    rhs_terms <- rhs_terms[rhs_terms != ""]  # remove empties
    rhs_terms <- unique(rhs_terms)           # remove duplicates
    
    # Handle interactions
    if (move_interactions_last) {
      # Identify terms with ':' or '*' or I( ) (simple heuristic)
      interaction_idx <- grepl("[:*]|I\\(", rhs_terms)
      rhs_terms <- c(rhs_terms[!interaction_idx], rhs_terms[interaction_idx])
    }
    
    # Rebuild formula
    rhs <- paste(rhs_terms, collapse = " + ")
    as.formula(paste(lhs, "~", rhs))
  }
  
  print("----Stop in analysis")
  print(formula.a)
  mf.a <- model.frame(formula.a, data)
  mfd0 <- mf.a[data[,delta]==0, !names(mf.a) %in% Tau]
  print(formula.e)
  print(formula.b)
  print(formula.c)
  f.e.c <- combine_formulas_clean(formula.e,formula.c)
  mf.e.c <- model.frame(f.e.c, data)
  mfd1 <- mf.e.c[data[,delta]==1,]
  k <- length(levels(model.extract(mf.e.c, "response")))
  print(f.e.c)
  cat(paste("Parameters estimation step...\n"))
  param.Weib.est <- fit.param.Wei.full.1(data)
  n.c.est <- n.c.estimation(f.e.c, formula.a, mfd0, mfd1, k, outcome.model)
  print(n.c.est)
  
  extract_params <- function(estimates, formula, include_intercept = FALSE, new_prefix = "gamma") {
    # Extract RHS term labels from formula
    terms <- attr(terms(formula), "term.labels")
    if("Tau:R" %in% terms) terms[terms == "Tau:R"] <- "R:Tau"
    # Escape special characters in term names (like ':')
    terms_escaped <- gsub("([\\.\\+\\-\\*\\^\\(\\)])", "\\\\\\1", terms)
    
    # Build regex pattern for predictors
    pattern <- paste0("^V\\.gamma\\.(", paste(terms_escaped, collapse = "|"), ")$")
    # Subset matching estimates for predictors
    matched <- estimates[grep(pattern, names(estimates))]
    # Reorder to match formula order
    matched <- matched[match(terms, gsub("^V\\.gamma\\.", "", names(matched)))]
    
    # Add intercepts if requested (any V.gamma.<number>)
    if (include_intercept) {
      intercepts <- estimates[grep("^V\\.gamma\\.[0-9]+$", names(estimates))]
      matched <- c(intercepts, matched)
    }
    # Simple rename: replace old prefix with new prefix
    names(matched) <- gsub("^V\\.gamma", paste0("V.", new_prefix), names(matched))
    return(matched)
  }
  
  a.par <- n.c.est$naive[,"Estimate"]
  a.length <- length(a.par)
  e.par <- extract_params(n.c.est$cc[,"Estimate"], formula.e,FALSE,"eta")
  e.length <- length(e.par)
  b.par <- extract_params(n.c.est$cc[,"Estimate"], formula.b, TRUE,"beta")
  b.length <- length(b.par)
  c.par <- extract_params(n.c.est$cc[,"Estimate"], formula.c, TRUE)
  c.length <- length(c.par)
  par <- c(a.par,e.par,b.par,c.par)
  l <- c(a.length,e.length,b.length,c.length)

  lower.b <- rep(-Inf,length(par))

  if(outcome.model == "PO" & k > 2){
    lower.b[2:(k-1)] <- lower.b[l[1] + l[2] + (2:(k-1))] <- lower.b[l[1] + l[2] + l[3] + (2:(k-1))] <- 1e-8
    par[2:(k-1)] <- diff(par[1:(k-1)])
    par[l[1] + l[2] + (2:(k-1))] <- diff(par[l[1] + l[2] + (1:(k-1))])
    par[l[1] + l[2] + l[3] + (2:(k-1))] <- diff(par[l[1] + l[2] + l[3] + (1:(k-1))])
  }
  
  print(par)
  print(l)
  print(lower.b)
  
  optim_result <- optim(par = par, fn = pseudo.neg.log.like,
                        weights = param.Weib.est$weights, k = k, lengths = l, delta = delta, V = names(mf.a)[1],
                        outcome.model = outcome.model, data = data, method = "L-BFGS-B",
                        lower = lower.b, control = list(trace = 3, maxit = 200, factr = 1e7, pgtol = 1e-12, lmm = 20))
  
  par.list <- list(d = param.Weib.est$est[0:length(attr(terms(cureform), "term.labels"))+1],
                   Tau = param.Weib.est$est[(length(attr(terms(cureform), "term.labels"))+2):length(param.Weib.est$est)],
                   a = optim_result$par[1:l[1]],
                   e = optim_result$par[(l[1]+1):(l[1]+l[2])],
                   b = optim_result$par[(l[1]+l[2]+1):(l[1]+l[2]+l[3])],
                   c = optim_result$par[(l[1]+l[2]+l[3]+1):sum(l)])
  if(l[2]==0) par.list$e <- NULL
  
  if(outcome.model == "PO" & k > 2){
    par.list$a[1:(k-1)] <- cumsum(par.list$a[1:(k-1)])
    par.list$b[1:(k-1)] <- cumsum(par.list$b[1:(k-1)])
    par.list$c[1:(k-1)] <- cumsum(par.list$c[1:(k-1)])
  }
  
  if(var){
    cat(paste("Variance estimation step...\n"))
    v.est <- variance.est(est = par.list, outcome.model = outcome.model, data = data, 
                          survform = survform, cureform = cureform, Tau = Tau, delta = delta, V = names(mf.a)[1],
                          lengths = l)
    
    cat(paste("Done!"))
    return(list(par = par.list,
                variance = v.est,
                value = optim_result$value, 
                naive = n.c.est$naive, 
                cc = n.c.est$cc))
  }
  
  cat(paste("Done!"))
  return(list(par = par.list, 
              value = optim_result$value, 
              naive = n.c.est$naive, 
              cc = n.c.est$cc))
  
  stop()
  print("----")
  
  
  
  # Complete case model
  f.c <- formula.c#update(formula.c, paste("~ . +",Tau,"+",R))
  # Naive case model
  f.a <- update(formula.a, paste("~ . +",R))
  
  if(any(vapply(model.frame(f.a, data), function(col) all(col == col[1]), logical(1))) ||
     any(vapply(model.frame(f.c, data), function(col) all(col == col[1]), logical(1)))){stop("one of the variables is constant")}
  
  mf.c <- model.frame(f.c, data)
  mf.a <- model.frame(f.a, data)
  mfd0 <- mf.a[data[,delta]==0, !names(mf.a) %in% Tau]
  mfd1 <- mf.c[data[,delta]==1,]
  Y <- model.extract(mf.c, "response")
  if(is.null(Y)){stop("no response variable found")}
  k <- length(levels(Y))
  if(k<3){stop("response must have 3 or more levels")} # I might have to change polr so k=2 is valid too
  
  cat(paste("Parameters estimation step...\n"))
  
  param.Weib.est <- fit.param.Wei.full.1(data)
  
  # call for naive and complete estimation function
  n.c.est <- n.c.estimation(f.c, f.a, mfd0, mfd1, k, outcome.model)
  
  beta.as.alpha.flag <- FALSE
  
  cat(paste("Optim...\n"))
  solution <- tryCatch({
    par.init <- par.initialization(equal.effect, n.c.est$naive[,"Estimate"], n.c.est$cc[,"Estimate"], k)
    l <- par.init$lengths

    optim_result <- NULL
    optim_result <- optim(par = par.init$par, fn = pseudo.neg.log.like,
                          weights = param.Weib.est$weights, k = k, lengths = l, delta = delta, V = names(mf.c)[1],
                          outcome.model = outcome.model, data = data, method = "L-BFGS-B",
                          lower = par.init$lower.b, control = list(trace = 3, maxit = 200)) #, factr = 1e-12
    list(optim = optim_result)}, error = function(e) {
      warning(paste("Second optim failure:", e$message))
      return(NULL)
    })
  
  if(outcome.model == "PO" & k > 2){
    solution$optim$par[1:(k-1)] <- cumsum(solution$optim$par[1:(k-1)])
    solution$optim$par[l[1] + (1:(k-1))] <- cumsum(solution$optim$par[l[1] + (1:(k-1))])
    solution$optim$par[l[1] + l[2] + (1:(k-1))] <- cumsum(solution$optim$par[l[1] + l[2] + (1:(k-1))])
  }

  par.list <- list(d = param.Weib.est$est[0:length(attr(terms(cureform), "term.labels"))+1],
                   Tau = param.Weib.est$est[(length(attr(terms(cureform), "term.labels"))+2):length(param.Weib.est$est)],
                   a = solution$optim$par[1:l[1]],
                   b = solution$optim$par[(l[1]+1):(l[1]+l[2])],
                   c = solution$optim$par[(l[1]+l[2]+1):(l[1]+l[2]+l[3])])
  if(length(l)==4) par.list$e <- solution$optim$par[(l[1]+l[2]+l[3]+1):sum(l)]
  
  if(var){
    cat(paste("Variance estimation step...\n"))
    v.est <- variance.est(est = par.list, outcome.model = outcome.model, data = as.data.frame(data), 
                          survform = survform, cureform = cureform, Tau = Tau, delta = delta, V = names(mf.c)[1],
                          lengths = l)
    cat(paste("Done!"))
    return(list(par = par.list,
                variance = v.est,
                value = solution$optim$value, 
                beta.as.alpha.flag = beta.as.alpha.flag, 
                naive = n.c.est$naive, 
                cc = n.c.est$cc))
  }
  cat(paste("Done!"))
  return(list(par = par.list, 
              value = solution$optim$value, 
              beta.as.alpha.flag = beta.as.alpha.flag, 
              naive = n.c.est$naive, 
              cc = n.c.est$cc))
}

