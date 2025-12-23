# Main fitting function for ordinal cure models

# ============================================================
# Validation utilities (fitting only)
# ============================================================

validate_ordcured_inputs <- function(
    data, survform, cureform,
    formula.a, formula.e, formula.b, formula.c,
    Tau, R, delta, outcome.model
) {
  
  # ============================================================
  # 1. Outcome model
  # ============================================================
  if (!outcome.model %in% c("PO", "ACAT")) {
    stop("outcome.model must be one of: 'PO', 'ACAT'")
  }
  
  # ============================================================
  # 2. Required variables exist in data
  # ============================================================
  all_forms <- list(survform, cureform, formula.a, formula.e, formula.b, formula.c)
  form_vars <- unique(unlist(lapply(all_forms, all.vars)))
  required_vars <- unique(c(form_vars, Tau, R, delta))
  
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ",
         paste(missing_vars, collapse = ", "))
  }
  
  # ============================================================
  # 3. Response consistency and type
  # ============================================================
  resp <- all.vars(formula.a)[1]
  
  same_resp <- vapply(
    list(formula.e, formula.b, formula.c),
    function(f) all.vars(f)[1] == resp,
    logical(1)
  )
  if (!all(same_resp)) {
    stop("All outcome formulas (a, e, b, c) must share the same response variable")
  }
  
  if (!is.ordered(data[[resp]])) {
    stop("Outcome variable must be an ordered factor")
  }
  
  # ============================================================
  # 4. Forbidden variable placement
  # ============================================================
  outcome_vars <- unique(unlist(lapply(
    list(formula.a, formula.e, formula.b, formula.c),
    all.vars
  )))
  
  if (delta %in% outcome_vars) {
    stop("delta cannot appear in outcome model formulas")
  }
  
  if (Tau %in% all.vars(formula.a)) {
    stop("Tau cannot appear in formula.a")
  }
  
  # ============================================================
  # 5. Basic semantic checks
  # ============================================================
  if (!is.numeric(data[[Tau]]) || any(data[[Tau]] < 0, na.rm = TRUE)) {
    stop("Tau must be a non-negative numeric variable")
  }
  
  if (!is.logical(data[[delta]]) &&
      !all(na.omit(unique(data[[delta]])) %in% c(0, 1))) {
    stop("delta must be logical or binary (0/1)")
  }
  
  invisible(TRUE)
}

# ============================================================
# Complete-Case formula handling
# ============================================================

build_cc_formula <- function(formula.e, formula.c, move_interactions_last = TRUE) {
  
  # ---- checks ----
  if (!inherits(formula.e, "formula") || !inherits(formula.c, "formula"))
    stop("formula.e and formula.c must be formula objects")
  
  lhs.e <- all.vars(formula.e)[1]
  lhs.c <- all.vars(formula.c)[1]
  
  if (lhs.e != lhs.c)
    stop("formula.e and formula.c must have the same response")
  
  # ---- extract RHS terms ----
  terms.e <- attr(terms(formula.e), "term.labels")
  terms.c <- attr(terms(formula.c), "term.labels")
  
  rhs <- unique(c(terms.e, terms.c))
  
  # ---- reorder interactions last ----
  if (move_interactions_last) {
    is_inter <- grepl("[:*]|I\\(", rhs)
    rhs <- c(rhs[!is_inter], rhs[is_inter])
  }
  
  # ---- rebuild formula ----
  if (length(rhs) == 0) {
    as.formula(paste(lhs.e, "~ 1"))
  } else {
    as.formula(paste(lhs.e, "~", paste(rhs, collapse = " + ")))
  }
}

fit_cure_weibull_em <- function(
    data, survform, cureform, Tau, delta,
    tol = 1e-8, maxit = 200, verbose = FALSE
) {
  
  n <- nrow(data)
  i0 <- data[[delta]] == 0

  ## ---- design matrices (cached) ----
  Xc <- model.matrix(cureform, data)
  Xt <- model.matrix(update(survform, . ~ . - 1), data)

  ## ---- initial values ----
  d.CC <- data[data[[delta]] == 1, , drop = FALSE]

  surv_init <- as.formula(
    paste0("Surv(time=", Tau, ", event=rep(1, nrow(d.CC))) ~ ",
           paste(colnames(Xt), collapse = " + "))
  )

  fit_lat <- phreg(surv_init, data = d.CC)
  fit_cure <- glm(cureform, family = binomial, data = data)

  b <- coef(fit_cure)
  beta <- fit_lat$coef[1:(length(fit_lat$coef)-2)]
  a <- exp(fit_lat$coef[length(fit_lat$coef)])
  k <- exp(-fit_lat$coef[length(fit_lat$coef) - 1])

  par <- c(b, beta, a, k)

  ## ---- log-likelihoods ----
  loglik_cure <- function(b, w) {
    eta <- Xc %*% b
    return(-sum(w * eta - log1p(exp(eta))))
  }

  loglik_lat <- function(theta, w) {
    beta <- theta[1:(length(theta)-2)]
    a <- theta[length(theta) - 1]
    k <- theta[length(theta)]

    eta <- Xt %*% beta
    t <- data[[Tau]]

    ll <- data[[delta]] * (eta + log(a) + (a - 1) * log(t) - a * log(k)) -
      exp(eta + a * log(t) - a * log(k))
    return(-sum(w * ll))
  }

  ## ---- EM loop ----
  loglik_prev <- Inf

  for (iter in seq_len(maxit)) {

    ## E-step
    pi <- plogis(Xc %*% b)
    eta <- exp(Xt %*% beta)
    s <- exp(-((data[[Tau]] / k)^a) * eta)

    w <- rep(1, n)
    w[i0] <- (pi[i0] * s[i0]) / (pi[i0] * s[i0] + 1 - pi[i0])

    ## M-step
    fit1 <- optim(b, loglik_cure, w = w, method = "BFGS")
    fit2 <- optim(
      c(beta, a, k), loglik_lat, w = w,
      method = "L-BFGS-B",
      lower = c(rep(-Inf, ncol(Xt)), 1e-6, 1e-6)
    )

    b <- fit1$par
    beta <- fit2$par[1:(length(fit2$par)-2)]
    a <- fit2$par[length(fit2$par) - 1]
    k <- fit2$par[length(fit2$par)]

    loglik <- fit1$value + fit2$value
    diff.L <- abs(loglik_prev - loglik)
    diff   <- sqrt(sum((par - c(b, beta, a, k))^2))
    par    <- c(b, beta, a, k)

    if (verbose)
      message("iter=", iter, " log.L=", -loglik, " diff=", diff, " diff.L=", diff.L)

    if (diff < tol) break
    loglik_prev <- loglik
  }

  ## ---- posterior weights ----
  wdf <- data.frame(
    w0 = ifelse(i0, 1 - pi, 1),
    w1 = ifelse(i0, pi * s, 1)
  )

  ## ---- parameter names ----

  # Cure model (logistic)
  cure_terms <- colnames(Xc)
  cure_terms[cure_terms == "(Intercept)"] <- "0"
  cure_names <- paste0("b.", cure_terms)

  # Latency model (PH part)
  lat_terms <- colnames(Xt)
  lat_names <- paste0("zeta.", lat_terms)

  # Final parameter names
  names(par) <- c(
    cure_names,
    lat_names,
    "shape",
    "scale"
  )

  list(
    par = par,
    weights = wdf,
    loglik = -loglik,
    iter = iter
  )
}


fit_naive_cc_models <- function(
    formula.a, formula.cc, data, delta, outcome.model
) {
  
  ## ---- basic checks ----
  if (!inherits(formula.a, "formula") || !inherits(formula.cc, "formula"))
    stop("formula.a and formula.cc must be formula objects")
  
  if (!outcome.model %in% c("PO", "ACAT"))
    stop("outcome.model must be 'PO' or 'ACAT'")
  
  ## ---- split data ----
  data0 <- data[data[[delta]] == 0, , drop = FALSE]
  data1 <- data[data[[delta]] == 1, , drop = FALSE]
  
  ## ---- outcome info ----
  k <- length(levels(data[[all.vars(formula.cc)[1]]]))
  
  vars.a  <- attr(terms(formula.a),  "term.labels")
  vars.cc <- attr(terms(formula.cc), "term.labels")
  
  ## ---- models fitting ----
  if (outcome.model == "PO" && k > 2) {
    naive.fit <- polr(formula.a, data = data0, Hess = TRUE)
    cc.fit <- polr(formula.cc, data = data1, Hess = TRUE)
    
    est.naive <- unname(c(naive.fit$zeta, -naive.fit$coefficients))
    est.cc <- unname(c(cc.fit$zeta, -cc.fit$coefficients))
    
    # Compute standard errors using vcov()
    std <- rep(0,length(est.naive))
    if(all(diag(vcov(naive.fit))>=0)) {std <- sqrt(diag(vcov(naive.fit)))}
    se.naive <- unname(c(tail(std, k-1), head(std, length(std) - k + 1)))
    
    std <- rep(0,length(est.cc))
    if(all(diag(vcov(cc.fit))>=0)) {std <- sqrt(diag(vcov(cc.fit)))}
    se.cc <- unname(c(tail(std, k-1), head(std, length(std) - k + 1)))
  } 
  else{
    naive.fit <- vglm(formula.a, acat(reverse = FALSE, parallel = TRUE), data = data0)
    cc.fit <- vglm(formula.cc, acat(reverse = FALSE, parallel = TRUE), data = data1)
    est.naive <- unname(coef(naive.fit))
    est.cc <- unname(coef(cc.fit))
    
    se.naive <- unname(sqrt(diag(vcov(naive.fit))))
    se.cc <- unname(sqrt(diag(vcov(cc.fit))))
  }
  
  names(est.naive) <- names(se.naive) <-
    c(paste0("V.alpha.", seq_len(k - 1)),
      paste0("V.alpha.", vars.a))
  
  names(est.cc) <- names(se.cc) <-
    c(paste0("V.gamma.", seq_len(k - 1)),
      paste0("V.gamma.", vars.cc))
  
  summary.n <- cbind(est.naive, se.naive)
  summary.c <- cbind(est.cc, se.cc)
  colnames(summary.n) <- colnames(summary.c) <- c("Estimate","Std. Error")
  
  ## ---- return ----
  list(
    naive = summary.n, 
    cc = summary.c
  )
}

extract_params <- function(estimates, formula, include_intercept = FALSE, prefix) {
  
  terms <- attr(terms(formula), "term.labels")
  terms[terms == "Tau:R"] <- "R:Tau"
  
  pat <- paste0("^V\\.gamma\\.(", paste(terms, collapse = "|"), ")$")
  matched <- estimates[grep(pat, names(estimates))]
  matched <- matched[match(terms, sub("^V\\.gamma\\.", "", names(matched)))]
  
  if (include_intercept) {
    intercepts <- estimates[grep("^V\\.gamma\\.[0-9]+$", names(estimates))]
    matched <- c(intercepts, matched)
  }
  
  names(matched) <- sub("^V\\.gamma", paste0("V.", prefix), names(matched))
  matched
}

build_initial_parameters <- function(nc, formulas, outcome.model, k) {
  
  a.par <- nc$naive[, "Estimate"]
  e.par <- extract_params(nc$cc[, "Estimate"], formulas$e, FALSE, "eta")
  b.par <- extract_params(nc$cc[, "Estimate"], formulas$b, TRUE,  "beta")
  c.par <- extract_params(nc$cc[, "Estimate"], formulas$c, TRUE,  "gamma")
  
  par <- c(a.par, e.par, b.par, c.par)
  lengths <- c(length(a.par), length(e.par), length(b.par), length(c.par))
  
  lower.b <- rep(-Inf, length(par))
  
  if (outcome.model == "PO" && k > 2) {
    idx <- function(i) seq_len(k - 1) + i
    par[idx(0)] <- diff(par[seq_len(k - 1)])
  }
  
  list(par = par, lengths = lengths, lower = lower.b)
}

# ordcure() fits an ordinal regression model with a cured fraction and a parametric survival component, 
# allowing different covariate effects before and after the event time.

ordcure <- function(
    survform, cureform,
    formula.a, formula.e, formula.b, formula.c,
    Tau, R, delta, data,
    outcome.model = c("PO", "ACAT"),
    var = FALSE
) {
  
  outcome.model <- match.arg(outcome.model)
  print("---- validation start ----")
  validate_ordcured_inputs(
    data, survform, cureform,
    formula.a, formula.e, formula.b, formula.c,
    Tau, R, delta, outcome.model
  )
  print("---- validation end ----")
  # ---- combine CC formula ----
  print("---- combine start ----")
  formula.cc <- build_cc_formula(formula.e, formula.c)
  print("---- combine end ----")
  
  # ---- first stage ----
  print("---- first stage start ----")
  stage1 <- fit_cure_weibull_em(data, survform, cureform, Tau, delta, tol = 1e-8, maxit = 200, verbose = FALSE)
  print(stage1$par)
  print("---- first stage end ----")
  
  # ---- naive + CC ----
  print("---- naive and cc start ----")
  nc <- fit_naive_cc_models(formula.a, formula.cc, data, delta, outcome.model)
  print("---- naive and cc end ----")
  
  # ---- init ----
  k <- length(levels(model.extract(model.frame(formula.cc, data), "response")))
  print("---- build initial start ----")
  
  init <- build_initial_parameters(
    nc,
    list(e = formula.e, b = formula.b, c = formula.c),
    outcome.model, k
  )
  print("---- build initial end ----")
  
  # ---- optim ----
  init$par
  print("---- optim start ----")
  opt <- optim(
    par = init$par,
    fn = negloglik_stage2,
    weights = stage1$weights,
    k = k,
    lengths = init$lengths,
    delta = delta,
    response = all.vars(formula.a)[1],
    outcome.model = outcome.model,
    data = data,
    method = "L-BFGS-B",
    lower = init$lower
  )
  print("---- optim end ----")
  l <- init$lengths

  par.list <- list(d = stage1$par[0:length(attr(terms(cureform), "term.labels"))+1],
                   Tau = stage1$par[(length(attr(terms(survform), "term.labels"))+2):length(stage1$par)],
                   a = opt$par[1:l[1]],
                   e = if(l[2]>0) opt$par[(l[1]+1):(l[1]+l[2])] else NULL,
                   b = opt$par[(l[1]+l[2]+1):(l[1]+l[2]+l[3])],
                   c = opt$par[(l[1]+l[2]+l[3]+1):sum(l)])

  
  if(outcome.model == "PO" & k > 2){
    par.list$a[1:(k-1)] <- cumsum(par.list$a[1:(k-1)])
    par.list$b[1:(k-1)] <- cumsum(par.list$b[1:(k-1)])
    par.list$c[1:(k-1)] <- cumsum(par.list$c[1:(k-1)])
  }
  print(par.list)
  out <- list(
    par.list = par.list,
    opt = opt,
    first.stage = stage1,
    naive = nc$naive,
    cc = nc$cc
  )
  
  if (var) {
    print("---- variance estimation start ----")
    out$variance <- variance.est(
      est = par.list,
      outcome.model = outcome.model,
      data = data,
      survform = survform,
      cureform = cureform,
      Tau = Tau,
      delta = delta,
      V = all.vars(formula.a)[1],
      lengths = l
    )
    print("---- variance estimation end ----")
  }
  
  class(out) <- "ordcure"
  out
}
