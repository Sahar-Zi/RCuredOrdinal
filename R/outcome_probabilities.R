.extract_terms <- function(param_names, prefix) {
  sub(paste0("^", prefix, "\\."), "", param_names)
}

.build_design <- function(data, terms) {
  X <- matrix(1, nrow(data), length(terms))
  colnames(X) <- terms

  for (j in seq_along(terms)) {
    vars <- strsplit(terms[j], ":", fixed = TRUE)[[1]]
    X[, j] <- Reduce(`*`, data[, vars, drop = FALSE])
  }
  
  X
}

.build_lp <- function(intercepts, slopes, data, k, prefix) {
  
  n <- nrow(data)
  lp <- matrix(intercepts, n, k - 1, byrow = TRUE)
  
  if (length(slopes) > 0) {
    terms <- .extract_terms(names(slopes), prefix)
    X <- .build_design(data, terms)
    B <- matrix(slopes, ncol=length(slopes), nrow = k - 1, byrow = TRUE)
    lp <- lp + X %*% t(B)
  }
  
  lp
}

#' Compute outcome probabilities for ordinal cure model
#'
#' @param par.list List with elements a, e, b, c (named numeric vectors)
#' @param data Data frame
#' @param delta Character name of censoring indicator
#' @param k Number of ordinal levels
#' @param outcome.model "PO" or "ACAT"
#'
#' @return List with matrices D0 and D1 (n x k)
compute_outcome_probs <- function(par.list,
                                  data,
                                  delta,
                                  k,
                                  outcome.model = c("PO", "ACAT")) {
  
  outcome.model <- match.arg(outcome.model)
  n <- nrow(data)
  delta0 <- data[[delta]] == 0
  
  ## ---- split intercepts / slopes ----
  split_par <- function(p) {
    list(
      intercept = p[seq_len(k - 1)],
      slope     = if (length(p) > (k - 1)) p[-seq_len(k - 1)] else numeric(0)
    )
  }
  
  A <- split_par(par.list$a)
  B <- split_par(par.list$b)
  C <- split_par(par.list$c)
  E <- par.list$e %||% numeric(0)
  
  ## ---- linear predictors ----
  eta.a <- .build_lp(A$intercept, A$slope, data, k, "V.alpha")
  eta.b <- .build_lp(B$intercept, B$slope, data, k, "V.beta")
  eta.c <- .build_lp(C$intercept, C$slope, data, k, "V.gamma")
  
  ## shared eta (eta parameters)
  if (length(E) > 0) {
    terms.e <- .extract_terms(names(E), "V.eta")
    X.e <- .build_design(data, terms.e)
    shared <- matrix(X.e %*% E, ncol = k-1, nrow = n, byrow = FALSE)
    eta.b <- eta.b + shared
    eta.c <- eta.c + shared
  }
  
  ## ---- probabilities ----
  if (outcome.model == "PO") {
    
    g <- function(x) 1 / (1 + exp(-x))
    
    G0 <- g(eta.a)
    G1 <- g(eta.c)
    G1[delta0, ] <- g(eta.b)[delta0, ]
    
    probs.0 <- cbind(G0, 1) - cbind(0, G0)
    probs.1 <- cbind(G1, 1) - cbind(0, G1)
    
  } else {  # ACAT
    
    build_acat <- function(x) {
      x <- cbind(1, exp(x))
      out <- matrix(1, n, k)
      out[, 2:k] <- x[, 2:k] * x[, 1:(k - 1)]
      out / rowSums(out)
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
  