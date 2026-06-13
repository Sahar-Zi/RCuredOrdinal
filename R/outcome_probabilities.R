# =====================================================================
# Outcome category probabilities for the ordinal cure model.
#
# Two link families are supported:
#   "PO"   - proportional-odds (cumulative logits)
#   "ACAT" - adjacent-category logits
# =====================================================================

# Strip a parameter-type prefix (e.g. "V.alpha.") from coefficient names,
# leaving the bare term labels ("Z1", "Z1:Z2", ...).
#' @noRd
.extract_terms <- function(param_names, prefix) {
  sub(paste0("^", prefix, "\\."), "", param_names)
}

# Build the design matrix for a set of term labels. Each interaction term
# (e.g. "R:Tau") becomes a single column equal to the product of its
# constituent columns. This assumes numeric covariates (continuous or 0/1
# dummies); it does not apply factor contrasts, so categorical predictors
# must be pre-coded as numeric columns in `data`.
#' @noRd
.build_design <- function(data, terms) {
  X <- matrix(1, nrow(data), length(terms))
  colnames(X) <- terms
  
  for (j in seq_along(terms)) {
    vars <- strsplit(terms[j], ":", fixed = TRUE)[[1]]
    X[, j] <- apply(data[, vars, drop = FALSE], 1, prod)
  }
  
  X
}

# Assemble the n x (k-1) linear predictor for one parameter block.
#   - `intercepts` are the k-1 category-specific cut points.
#   - `slopes` are covariate effects, shared across categories
#     (the parallel / proportional assumption), broadcast over the k-1
#     columns.
# Under "PO" the slope and shared-eta contributions are negated so that a
# positive coefficient raises the odds of a higher category.
#' @noRd
.build_lp <- function(intercepts, slopes, data, k, prefix, outcome.model) {
  
  n <- nrow(data)
  lp <- matrix(intercepts, n, k - 1, byrow = TRUE)
  
  if (length(slopes) > 0) {
    terms <- .extract_terms(names(slopes), prefix)
    X <- .build_design(data, terms)
    B <- matrix(slopes, ncol = length(slopes), nrow = k - 1, byrow = TRUE)
    if (outcome.model == "PO")
      B <- -B
    lp <- lp + X %*% t(B)
  }
  
  lp
}

#' Compute outcome probabilities for the ordinal cure model
#'
#' Internal helper. Returns the per-observation category probabilities for
#' the two components of the pseudo-likelihood mixture:
#'   - `D0`: probabilities under the cured ("alpha") parameters, used for
#'     the cured part of censored observations.
#'   - `D1`: probabilities under the "gamma" parameters for observed
#'     (`delta = 1`) rows, and under the "beta" parameters for censored
#'     (`delta = 0`) rows.
#' The shared "eta" covariate effects (when present) are added to both the
#' beta and gamma linear predictors.
#'
#' @param par List with elements `a`, `e`, `b`, `c` (named numeric vectors).
#' @param data A data frame.
#' @param delta Character name of the event-indicator column.
#' @param k Number of ordinal levels.
#' @param outcome.model `"PO"` or `"ACAT"`.
#'
#' @return A list with matrices `D0` and `D1`, each `n x k`, whose rows sum
#'   to one.
#'
#' @noRd
compute_outcome_probs <- function(par,
                                  data,
                                  delta,
                                  k,
                                  outcome.model = c("PO", "ACAT")) {
  
  outcome.model <- match.arg(outcome.model)
  n <- nrow(data)
  delta0 <- data[[delta]] == 0
  
  ## Split a parameter block into its k-1 intercepts and the remaining slopes.
  split_par <- function(p) {
    list(
      intercept = p[seq_len(k - 1)],
      slope     = if (length(p) > (k - 1)) p[-seq_len(k - 1)] else numeric(0)
    )
  }
  
  A <- split_par(par$a)   # cured component ("alpha")
  B <- split_par(par$b)   # censored component ("beta")
  C <- split_par(par$c)   # observed component ("gamma")
  E <- if (is.null(par$e)) numeric(0) else par$e   # shared effects ("eta")
  
  ## ---- linear predictors (n x (k-1)) ----
  eta.a <- .build_lp(A$intercept, A$slope, data, k, "V.alpha", outcome.model)
  eta.b <- .build_lp(B$intercept, B$slope, data, k, "V.beta" , outcome.model)
  eta.c <- .build_lp(C$intercept, C$slope, data, k, "V.gamma", outcome.model)
  
  ## Shared "eta" covariate effects, added to both beta and gamma predictors.
  if (length(E) > 0) {
    terms.e <- .extract_terms(names(E), "V.eta")
    X.e <- .build_design(data, terms.e)
    shared <- matrix(X.e %*% E, ncol = k - 1, nrow = n, byrow = FALSE)
    if (outcome.model == "PO")
      shared <- -shared
    eta.b <- eta.b + shared
    eta.c <- eta.c + shared
  }
  
  ## ---- category probabilities ----
  if (outcome.model == "PO") {
    
    ## Cumulative-logit (proportional-odds) link.
    ## plogis(eta)[, j] = P(Y <= j); differencing the cumulative CDF
    ## (padded with 0 on the left and 1 on the right) gives the per-category
    ## probabilities.
    G0 <- plogis(eta.a)
    G1 <- plogis(eta.c)
    G1[delta0, ] <- plogis(eta.b)[delta0, ]
    
    probs.0 <- cbind(G0, 1) - cbind(0, G0)
    probs.1 <- cbind(G1, 1) - cbind(0, G1)
    
  } else {  # ACAT
    
    ## Adjacent-category link. eta[, j] = log( P(Y = j+1) / P(Y = j) ), so
    ## the unnormalized category probabilities are the cumulative products
    ##   u_1 = 1,   u_j = prod_{l < j} exp(eta_l) = exp(cumsum(eta)[j-1]),
    ## which we then normalize to sum to one. The cumulative sum (rather than
    ## a product of adjacent terms) is essential for k >= 4.
    build_acat <- function(eta) {
      csum <- if (ncol(eta) > 1) t(apply(eta, 1, cumsum)) else eta
      num  <- cbind(1, exp(csum))
      num / rowSums(num)
    }
    
    probs.0 <- build_acat(eta.a)
    probs.1 <- build_acat(eta.c)
    probs.1[delta0, ] <- build_acat(eta.b)[delta0, ]
  }
  
  list(
    D0 = probs.0,
    D1 = probs.1
  )
}