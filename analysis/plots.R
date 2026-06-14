############################################################
## Plotting for run_simulation() output (class "ordcure_sim")
##
## plot(sim)                      -> bias + 95% CI by coefficient
## plot(sim, type = "coverage")   -> empirical coverage of 95% CIs
##
## Both compare the Proposed, Naive, and Complete-Case estimators.
## Coefficient grouping and labels are derived from the data, so the
## plots adapt to whatever model the simulation used.
##
## Requires ggplot2 (and scales for the coverage plot's nonlinear axis).
############################################################

## ---- internal helpers -----------------------------------

# Per-coefficient bias (with averaged Wald CI) or empirical coverage.
# `est` and `se` are reps x p matrices with named columns; `true` is a
# named vector aligned to those columns.
.sim_metric <- function(est, se, true, type) {
  
  if (!is.null(names(true)))
    true <- true[colnames(est)]
  
  reps <- nrow(est)
  z <- qnorm(0.975)
  
  if (type == "bias") {
    bias <- colMeans(est, na.rm = TRUE) - true
    lo   <- colMeans(est - z * se, na.rm = TRUE) - true
    hi   <- colMeans(est + z * se, na.rm = TRUE) - true
  } else {
    R  <- matrix(true, reps, length(true), byrow = TRUE)
    cp <- colMeans((R >= est - z * se) & (R <= est + z * se), na.rm = TRUE)
    se_cp <- sqrt(cp * (1 - cp) / reps)
    bias <- cp
    lo   <- pmax(0, cp - z * se_cp)
    hi   <- pmin(1, cp + z * se_cp)
  }
  
  data.frame(coef = colnames(est), y = bias, lo = lo, hi = hi,
             row.names = NULL, stringsAsFactors = FALSE)
}

# Which parameter block a coefficient name belongs to.
.sim_block <- function(name) {
  ifelse(grepl("^(d|Tau)\\.", name), "first",
         ifelse(grepl("V\\.alpha", name),   "alpha",
                ifelse(grepl("V\\.eta",   name),   "eta",
                       ifelse(grepl("V\\.beta",  name),   "beta",
                              ifelse(grepl("V\\.gamma", name),   "gamma", "other")))))
}

# Default plotmath label for a coefficient, e.g. "a.V.alpha.R:Tau" -> gamma[RT].
.sim_label <- function(name) {
  core  <- sub("^[a-z]\\.V\\.", "", name)            # "alpha.R:Tau"
  parts <- strsplit(core, ".", fixed = TRUE)[[1]]
  greek <- parts[1]                                  # alpha/beta/gamma/eta
  suff  <- paste(parts[-1], collapse = ".")
  suff  <- gsub("Tau", "T", suff, fixed = TRUE)
  suff  <- gsub(":", "", suff, fixed = TRUE)
  if (suff == "") greek else sprintf("%s[%s]", greek, suff)
}

# Assemble the long data frame plotted by plot.ordcure_sim().
# Separated out so the alignment / metric logic is testable without ggplot2.
.sim_plot_data <- function(x, type = c("bias", "coverage"),
                           se = c("sandwich", "inv")) {
  
  type <- match.arg(type)
  se   <- match.arg(se)
  
  if (is.null(x$est.se.values))
    stop("plotting needs a run with var = TRUE")
  
  est_se <- if (se == "inv") x$est.se.inv.values else x$est.se.values
  
  ## proposed estimator, restricted to second-stage coefficients
  prop <- .sim_metric(x$est.values, est_se, unlist(x$true.params), type)
  prop <- prop[!grepl("^(d|Tau)\\.", prop$coef), , drop = FALSE]
  prop$method <- "Proposed"
  
  ## x-axis block order: alpha, beta, gamma, then eta on the far right
  coef_order <- prop$coef
  blk <- .sim_block(coef_order)
  coef_order <- coef_order[order(match(blk, c("alpha", "beta", "gamma", "eta",
                                              "first", "other")))]
  prop <- prop[match(coef_order, prop$coef), , drop = FALSE]
  key <- sub("^[a-z]\\.", "", coef_order)            # V.alpha.* / V.eta.* / ...
  
  ## a competitor (naive or CC) contributes only to its own blocks; the
  ## rest of the rows stay NA so the methods line up on a shared x-axis.
  competitor <- function(vals, ses, truth, eligible, label) {
    out <- data.frame(coef = coef_order, y = NA_real_, lo = NA_real_,
                      hi = NA_real_, method = label, stringsAsFactors = FALSE)
    if (!is.null(vals)) {
      kk   <- key[eligible]
      have <- kk %in% colnames(vals)
      kk   <- kk[have]
      if (length(kk) > 0) {
        m <- .sim_metric(vals[, kk, drop = FALSE], ses[, kk, drop = FALSE],
                         truth[kk], type)
        rows <- which(eligible)[have]
        out[rows, c("y", "lo", "hi")] <- m[, c("y", "lo", "hi")]
      }
    }
    out
  }
  
  naive <- competitor(x$naive.est.values, x$naive.se, x$true.params$a,
                      grepl("^V\\.alpha\\.", key), "Naive")
  cc    <- competitor(x$cc.est.values, x$cc.se,
                      c(x$true.params$c, x$true.params$e),
                      grepl("^V\\.(gamma|eta)\\.", key), "Complete-Case")
  
  df <- rbind(prop[, c("coef", "y", "lo", "hi", "method")], naive, cc)
  df$coef   <- factor(df$coef, levels = coef_order)
  df$method <- factor(df$method, levels = c("Proposed", "Naive", "Complete-Case"))
  df
}

## ---- S3 method ------------------------------------------

#' Plot simulation results
#'
#' @param x An object returned by run_simulation().
#' @param type "bias" (default) or "coverage".
#' @param se "sandwich" (default) or "inv" standard errors for the proposed
#'   estimator.
#' @param labels Optional named vector/list of expressions overriding the
#'   default Greek axis labels.
#' @param ... Ignored.
#' @return A ggplot object.
#' @exportS3Method
plot.ordcure_sim <- function(x, type = c("bias", "coverage"),
                             se = c("sandwich", "inv"),
                             labels = NULL, ...) {
  
  type <- match.arg(type)
  se   <- match.arg(se)
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot.ordcure_sim()")
  
  df <- .sim_plot_data(x, type, se)
  coef_order <- levels(df$coef)
  
  ## subtitle describing the run (outcome model, replications, n)
  subtitle <- NULL
  if (!is.null(x$outcome.model)) {
    om <- c(PO = "Proportional-odds", ACAT = "Adjacent-category")[x$outcome.model]
    subtitle <- paste0(om, " outcome model")
    if (!is.null(x$n_success)) subtitle <- paste0(subtitle, ", ", x$n_success, " replications")
    if (!is.null(x$n))         subtitle <- paste0(subtitle, ", n = ", x$n)
  }
  
  ## block shading rectangles (blocks are contiguous along the x-axis)
  blocks <- .sim_block(coef_order)
  runs   <- rle(blocks)
  ends   <- cumsum(runs$lengths)
  rect_df <- data.frame(
    xmin  = c(0, head(ends, -1)) + 0.5,
    xmax  = ends + 0.5,
    block = runs$values
  )
  fills <- c(alpha = "#c6dbef", eta = "#dadaeb", beta = "#fdd0a2",
             gamma = "#d9f0d3", first = "#f0f0f0", other = "#f7f7f7")
  
  ## x-axis labels
  lab_fun <- if (is.null(labels)) {
    function(b) parse(text = vapply(as.character(b), .sim_label, character(1)))
  } else {
    function(b) labels[as.character(b)]
  }
  
  cols   <- c(Proposed = "#0072B2", Naive = "#D55E00", `Complete-Case` = "#CC79A7")
  shapes <- c(Proposed = 16,        Naive = 17,        `Complete-Case` = 15)
  leg    <- c(Proposed = "2-stage Proposed", Naive = "Naive",
              `Complete-Case` = "Complete Case")
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x = coef, y = y,
                                        color = method, shape = method)) +
    ggplot2::geom_rect(
      data = rect_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                   fill = block), alpha = 0.25, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = fills) +
    ggplot2::geom_vline(xintercept = head(ends, -1) + 0.5,
                        linetype = "dashed", color = "red") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi),
                           width = 0.35, linewidth = 0.8, alpha = 0.8,
                           na.rm = TRUE) +
    ggplot2::geom_point(size = 4, na.rm = TRUE) +
    ggplot2::scale_x_discrete(position = "top", labels = lab_fun) +
    ggplot2::scale_color_manual(values = cols, labels = leg) +
    ggplot2::scale_shape_manual(values = shapes, labels = leg) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      legend.position = "top",
      plot.caption = ggplot2::element_text(hjust = 0),
      axis.text.x  = ggplot2::element_text(size = 14),
      legend.background = ggplot2::element_rect(fill = "white", color = "black",
                                                linewidth = 0.5)
    )
  
  if (type == "bias") {
    se_label <- c(sandwich = "sandwich", inv = "inverse-information")[[se]]
    g + ggplot2::geom_hline(yintercept = 0, linetype = "longdash") +
      ggplot2::labs(x = NULL, y = "Bias",
                    title = "Bias and 95% confidence intervals by coefficient",
                    subtitle = subtitle,
                    caption = paste0("Proposed-estimator CIs use ", se_label,
                                     " standard errors."))
  } else {
    if (!requireNamespace("scales", quietly = TRUE))
      stop("scales is required for the coverage plot")
    
    ## stretch the tails (0-0.2 and 0.8-1), compress the middle
    tr <- scales::trans_new(
      "cp_scale",
      transform = function(v) ifelse(v <= 0.2, v * 2,
                                     ifelse(v <= 0.8, 0.4 + (v - 0.2) * 0.5,
                                            0.7 + (v - 0.8) * 2)),
      inverse   = function(v) ifelse(v <= 0.4, v / 2,
                                     ifelse(v <= 0.7, 0.2 + (v - 0.4) / 0.5,
                                            0.8 + (v - 0.7) / 2))
    )
    
    g + ggplot2::geom_hline(yintercept = 0.95, linetype = "longdash") +
      ggplot2::scale_y_continuous(
        transform = tr, limits = c(0, 1),
        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.95, 1)) +
      ggplot2::labs(x = NULL, y = "Coverage probability\n(95% CI)",
                    title = "Empirical coverage of 95% confidence intervals",
                    subtitle = subtitle)
  }
}