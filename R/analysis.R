ordcuredfit.param <- function(survform, cureform, formula.a, formula.e, formula.b, formula.c, Tau, R, delta, data, outcome.model=c("PO","ACAT"), var = var) {
  
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
    
    ## ---- Model validity checks ----
    
    if (missing(formula.a) || missing(formula.c))
      stop("formula.a and formula.c must be supplied")
    
    if (!inherits(formula.a, "formula") || !inherits(formula.c, "formula"))
      stop("formula.a and formula.c must be formula objects")
    
    lhs <- function(f) all.vars(f)[1]
    
    if (lhs(formula.a) != lhs(formula.c))
      stop("Response must be identical in formula.a and formula.c")
    
    if (all.vars(formula.a)[1] != all.vars(formula.c)[1])
      stop("Response must be identical in formula.a and formula.c")
    
    if (delta %in% c(all.vars(formula.a), all.vars(formula.c)))
      stop("delta cannot appear in model formulas")
    
    if (Tau %in% all.vars(formula.a))
      stop("Tau cannot appear in formula.a")
    
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

