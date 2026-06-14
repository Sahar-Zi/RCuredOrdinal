## ---- Setup ------------------------------------------------

source("analysis//_setup.R")

## --------------------------------------------------
## Load local config
## --------------------------------------------------
config_file <- file.path(getwd(), "config.R")

if (!file.exists(config_file)) {
  stop(
    "Missing config.R\n",
    "Create RCuredOrdinal/config.R with:\n",
    "PERTHES_DATA_PATH <- 'full/path/to/data.txt'"
  )
}

source(config_file)

if (!file.exists(PERTHES_DATA_PATH)) {
  stop("Data file not found: ", PERTHES_DATA_PATH)
}

## --------------------------------------------------
## Load data
## --------------------------------------------------

data <- read.table(file = PERTHES_DATA_PATH, sep = "\t", na = "", header = TRUE)
data$age.THR <- pmin(data$age.THR, data$age, na.rm=TRUE)

# modify the data in 32 observations where age.THR=age, shift their age.THR earlier within the same year:
set.seed(17)
data.eps <- data
num0 <- sum(data.eps$THR==1 & data.eps$age.THR == data.eps$age)
data.eps[data.eps$THR==1 & data.eps$age.THR == data.eps$age,"age.THR"] <- 
  data.eps[data.eps$THR==1 & data.eps$age.THR == data.eps$age,"age.THR"] - runif(num0, 0.1, 0.9)

X <- data.eps[, c("age.THR", "age","THR", "sex", "unilateral", "dysplasia",
                  "diag_age6_7", "diag_age8_11", "diag_age11", 
                  "chPtrt.surgery", "chPtrt.activity.restrict", "sf36_RP")]
X$delta <- X$THR

# Make sure it's still ordered
X$sf36_RP <- ordered(X$sf36_RP, levels = sort(unique(X$sf36_RP)))

# Collapse three middle levels
X$sf36_RP <- fct_collapse(X$sf36_RP, "25-75" = c(levels(X$sf36_RP)[2:4]))

## --------------------------------------------------
## Model specification
## --------------------------------------------------

OUTCOME_MODEL <- "ACAT"   # "PO" or "ACAT"
AGE_THR_FIXED <- 35

sf36_levels <- levels(X$sf36_RP)
K <- 3                    # fixed (gradient helpers assume k = 3)

cureform  <- THR ~ sex
survform  <- age.THR ~ sex + unilateral + diag_age6_7 + diag_age8_11 +
  diag_age11 + chPtrt.surgery + chPtrt.activity.restrict + sex:unilateral
formula.a <- sf36_RP ~ age
formula.b <- sf36_RP ~ age
formula.c <- sf36_RP ~ age + age.THR
formula.e <- sf36_RP ~ sex + dysplasia + unilateral + chPtrt.surgery +
  chPtrt.activity.restrict + sex:age

cat("--------------------------------\n", OUTCOME_MODEL, "outcome regression model\n")

fit <- ordcure(
  cureform      = cureform,
  survform      = survform,
  formula.a     = formula.a,
  formula.e     = formula.e,
  formula.b     = formula.b,
  formula.c     = formula.c,
  Tau           = "age.THR",
  R             = "age",
  delta         = "THR",
  data          = X,
  outcome.model = OUTCOME_MODEL,
  var           = TRUE,
  verbose       = TRUE
)

## ---- Summary table ----
print(summary(fit))

## --------------------------------------------------
## Quantities derived from the fit (no hard-coded indices)
## --------------------------------------------------

resp <- fit$vars$response           # outcome name, e.g. "sf36_RP"
Rvar <- fit$vars$R                  # continuous covariate, e.g. "age"
Tauv <- fit$vars$Tau                # censored covariate, e.g. "age.THR"

cf   <- coef(fit)                   # full named vector (d, Tau, a, e, b, c)
nm.a <- grep("^a\\.", names(cf), value = TRUE)   # alpha  (cured)
nm.e <- grep("^e\\.", names(cf), value = TRUE)   # eta    (shared)
nm.b <- grep("^b\\.", names(cf), value = TRUE)   # beta   (uncured, R <= Tau)
nm.c <- grep("^c\\.", names(cf), value = TRUE)   # gamma  (uncured, R >  Tau)

Vhat    <- vcov(fit)                              # name-labelled sandwich
V.alpha <- Vhat[nm.a,          nm.a]              # was [13:15, 13:15]
V.beta  <- Vhat[c(nm.b, nm.e),  c(nm.b, nm.e)]    # was [c(22:24,16:21), ...]
V.gamma <- Vhat[c(nm.c, nm.e),  c(nm.c, nm.e)]    # was [c(25:28,16:21), ...]

len.b <- length(fit$par.list$b)
len.e <- length(fit$par.list$e)
len.c <- length(fit$par.list$c)
nbe   <- len.b + len.e            # length of the (beta, eta) gradient block
ngce  <- len.c + len.e            # length of the (gamma, eta) gradient block

## covariate values (no intercept) for a set of "PREFIX.var" coef names,
## evaluated on one data row, in coefficient order.
covar_row <- function(covar_names, prefix_regex, row) {
  vars <- sub(prefix_regex, "", covar_names)
  mm   <- model.matrix(reformulate(vars), data = row)
  as.numeric(mm[, -1, drop = FALSE])
}

# ----------------------- Create data grid for plots -------------------------

grid.data <- expand.grid(age = seq(from = 18, to = 76, by = 0.2),
                         THR = c(0, 1),
                         sex = c(0, 1),
                         sf36_RP_V = sf36_levels) %>%      # non-constant variables
  arrange(sex, THR, age) %>%
  mutate(
    age.THR                  = AGE_THR_FIXED,
    dysplasia                = 0,
    unilateral               = 0,
    chPtrt.activity.restrict = 0,
    chPtrt.surgery           = 0,
    sf36_RP                  = match(sf36_RP_V, sf36_levels),
    delta                    = as.integer(THR == 1 & age > age.THR)) # constant variables

grid.data$sf36_RP <- ordered(grid.data$sf36_RP)   # keep the outcome ordered
cat.idx <- as.integer(grid.data$sf36_RP)          # category index per row

# ----------------------- 2-stage outcome probabilities (S3 predict) ---------
# delta = "delta" points predict at the "R > Tau" regime column built above,
# rather than the fit's own event column.

P0 <- predict(fit, newdata = grid.data, component = "cured",   delta = "delta")  # alpha
P1 <- predict(fit, newdata = grid.data, component = "uncured", delta = "delta")  # beta/gamma

is0 <- grid.data$THR == 0
grid.data$prob <- NA_real_
grid.data$prob[is0]  <- P0[cbind(which(is0),  cat.idx[is0])]
grid.data$prob[!is0] <- P1[cbind(which(!is0), cat.idx[!is0])]

# ----------------------- Naive & CC models (formulas from the fit) ----------
# naive  = the cured (alpha) outcome model, fit ignoring cure status on THR = 0
# CC     = the uncured outcome model (gamma + eta terms), fit on THR = 1
naive.form <- formula.a
cc.form    <- reformulate(
  c(attr(terms(formula.c), "term.labels"),
    attr(terms(formula.e), "term.labels")),
  response = resp)

fit_outcome <- function(form, dat) {
  if (OUTCOME_MODEL == "PO")   return(polr(form, data = dat))
  if (OUTCOME_MODEL == "ACAT") return(vglm(form, acat(reverse = FALSE, parallel = TRUE), data = dat))
  stop("Unknown OUTCOME_MODEL")
}
outcome_probs <- function(mod, newdata) {
  ty <- if (OUTCOME_MODEL == "PO") "probs" else "response"
  predict(mod, newdata = newdata, type = ty)
}

## Coefficients of a fitted naive/CC outcome model in (intercept1, intercept2,
## covariates) order -- the layout grad_naive_cc() expects on INPUT. polr
## stores intercepts separately in $zeta and lists covariates first; vglm/acat
## already returns intercepts first. (The gradient's OUTPUT order still matches
## each model's vcov(): covariates-last for ACAT, covariates-first for PO.)
coef_int_first <- function(mod) {
  if (OUTCOME_MODEL == "PO") c(mod$zeta, mod$coefficients) else coef(mod)
}

# Naive
grid.data.naive <- grid.data %>% filter(THR == 0)
mod.naive   <- fit_outcome(naive.form, X %>% filter(THR == 0))
naive.probs <- outcome_probs(mod.naive, grid.data.naive)
grid.data.naive$prob <- naive.probs[cbind(seq_len(nrow(grid.data.naive)), grid.data.naive$sf36_RP)]

# Complete case
grid.data.cc <- grid.data %>% filter(THR == 1)
mod.cc   <- fit_outcome(cc.form, X %>% filter(THR == 1))
cc.probs <- outcome_probs(mod.cc, grid.data.cc)
grid.data.cc$prob <- cc.probs[cbind(seq_len(nrow(grid.data.cc)), grid.data.cc$sf36_RP)]

# ==========================================================================
#  Analytic delta-method gradients (kept; generalized to the formulas / link)
# ==========================================================================

## ---- naive / CC predicted-probability gradient (k = 3) --------------------
grad_naive_cc <- function(alpha, x, k) {
  a1 <- alpha[1]
  a2 <- alpha[2]
  ax <- 0
  if (length(alpha) > 2) ax <- alpha[3:length(alpha)]
  
  if (OUTCOME_MODEL == "ACAT") {
    eta1 <- a1 + ax %*% x
    eta2 <- a2 + ax %*% x
    eta3 <- a1 + a2 + 2 * ax %*% x
    D <- 1 + exp(eta1) + exp(eta3)
  } else {
    eta1 <- a1 - ax %*% x
    eta2 <- a2 - ax %*% x
  }
  
  if (k == 1) {
    if (OUTCOME_MODEL == "ACAT") {
      grad <- c(d_1 = -(exp(eta1) + exp(eta3)) / D^2,
                d_2 = -exp(eta3) / D^2,
                d_x = -x * as.vector((exp(eta1) + 2 * exp(eta3)) / D^2))
    } else {
      grad <- c(-x, 1, 0) * as.vector(exp(eta1) / (1 + exp(eta1))^2)
    }
  }
  if (k == 2) {
    if (OUTCOME_MODEL == "ACAT") {
      grad <- c(d_1 = exp(eta1) / D^2,
                d_2 = -(exp(eta1) * exp(eta3)) / D^2,
                d_x = x * as.vector((exp(eta1) - exp(eta1) * exp(eta3)) / D^2))
    } else {
      grad <- c(d_x = x * as.vector((exp(eta1) / (1 + exp(eta1))^2) - (exp(eta2) / (1 + exp(eta2))^2)),
                d_1 = -exp(eta1) / (1 + exp(eta1))^2,
                d_2 =  exp(eta2) / (1 + exp(eta2))^2)
    }
  }
  if (k == 3) {
    if (OUTCOME_MODEL == "ACAT") {
      grad <- c(d_1 = exp(eta3) / D^2,
                d_2 = (exp(eta3) + (exp(eta1) * exp(eta3))) / D^2,
                d_x = x * as.vector((exp(eta1) * exp(eta3) + 2 * exp(eta3)) / D^2))
    } else {
      grad <- c(-x, 0, 1) * as.vector(exp(eta2) / (1 + exp(eta2))^2)
    }
  }
  grad
}

## ---- 2-stage cured (alpha) gradient ---------------------------------------
grad.cured.2stage <- function(par, data, k, resp) {
  a.inter <- par[1:(k - 1)]
  a.cvars <- par[k:length(par)]
  mm.a <- model.matrix(reformulate(sub("V.alpha.", "", names(a.cvars))), data = data)
  
  if (OUTCOME_MODEL == "PO") {
    xi.a <- mm.a %*% t(cbind(a.inter, matrix(-a.cvars, nrow = k - 1, ncol = length(a.cvars), byrow = TRUE)))
    base.g.D0 <- 1 / (1 + exp(xi.a))
    probs.0 <- cbind(1, base.g.D0) - cbind(base.g.D0, 0)
  }
  if (OUTCOME_MODEL == "ACAT") {
    xi.a <- mm.a %*% t(cbind(a.inter, matrix(a.cvars, nrow = k - 1, ncol = length(a.cvars), byrow = TRUE)))
    eta.a <- cbind(1, exp(xi.a))
    out <- matrix(1, nrow(data), k)
    out[, 2:k] <- eta.a[, 2:k] * eta.a[, 1:(k - 1)]
    probs.0 <- out / rowSums(out)
  }
  probs.0[1, as.integer(data[[resp]])]
}

## ---- 2-stage uncured (beta / gamma + eta) gradient ------------------------
grad.uncured.2stage <- function(par, data.m, k, nbe, ngce, Rvar, Tauv, resp) {
  
  split_par <- function(par_block) list(
    intercept = par_block[seq_len(k - 1)],
    covars    = if (length(par_block) > (k - 1)) par_block[-seq_len(k - 1)] else numeric(0))
  
  par.b <- split_par(par[1:nbe])
  par.c <- split_par(par[(nbe + 1):(nbe + ngce)])
  
  build_mm <- function(cvars, data, prefix) {
    if (length(cvars) == 0) return(matrix(1, nrow(data), 1))
    vars <- sub(prefix, "", names(cvars))
    model.matrix(reformulate(vars), data = data)
  }
  
  mm.b <- build_mm(par.b$covars, data.m, "V.(beta|eta).")
  mm.c <- build_mm(par.c$covars, data.m, "V.(gamma|eta).")
  
  compute_xi <- function(mm, intercept, covars) {
    if (OUTCOME_MODEL == "PO")
      coef_mat <- cbind(intercept, matrix(-covars, nrow = k - 1, ncol = length(covars), byrow = TRUE))
    if (OUTCOME_MODEL == "ACAT")
      coef_mat <- cbind(intercept, matrix( covars, nrow = k - 1, ncol = length(covars), byrow = TRUE))
    mm %*% t(coef_mat)
  }
  
  xi.b <- compute_xi(mm.b, par.b$intercept, par.b$covars)
  xi.c <- compute_xi(mm.c, par.c$intercept, par.c$covars)
  
  delta0 <- data.m[[Rvar]] <= data.m[[Tauv]]
  
  if (OUTCOME_MODEL == "PO") {
    g <- plogis
    G <- g(xi.c); G[delta0, ] <- g(xi.b)[delta0, ]
    probs <- cbind(G, 1) - cbind(0, G)
    return(probs[cbind(seq_len(nrow(data.m)), as.integer(data.m[[resp]]))])
  }
  if (OUTCOME_MODEL == "ACAT") {
    g <- function(z) {
      ez <- exp(z); out <- cbind(1, ez)
      out[, 2:k] <- out[, 2:k] * out[, 1:(k - 1)]
      out / rowSums(out)
    }
    G <- g(xi.c); G[delta0, ] <- g(xi.b)[delta0, ]
    return(G[cbind(seq_len(nrow(data.m)), as.integer(data.m[[resp]]))])
  }
  stop("Unknown OUTCOME_MODEL")
}

# ----------------------- CI for naive --------------------------------------
# Delta-method band for the naive fit: gradient evaluated at the naive model's
# OWN coefficients and combined with vcov(mod.naive), so the SE is coherent
# with the point estimate (prob = predict(mod.naive)).
v.a <- vcov(mod.naive)
coef.naive <- coef_int_first(mod.naive)
naive.df.plot <- grid.data.naive[, c("age", "sex", "sf36_RP", "prob")]
naive.df.plot$SE <- naive.df.plot$asymp.LCL <- naive.df.plot$asymp.UCL <- NA_real_

for (i in seq_len(nrow(grid.data.naive))) {
  ki   <- as.integer(grid.data.naive$sf36_RP[i])
  xi   <- covar_row(names(fit$par.list$a)[-(1:(K - 1))], "V.alpha.", grid.data.naive[i, , drop = FALSE])
  gr.a <- grad_naive_cc(coef.naive, x = xi, k = ki)
  se   <- sqrt(as.numeric(t(gr.a) %*% v.a %*% gr.a))
  naive.df.plot$SE[i]        <- se
  naive.df.plot$asymp.LCL[i] <- pmax(0, naive.df.plot$prob[i] - qnorm(0.975) * se)
  naive.df.plot$asymp.UCL[i] <- pmin(1, naive.df.plot$prob[i] + qnorm(0.975) * se)
}
naive.df.plot$df     <- Inf
naive.df.plot$Method <- "Naive"

# ----------------------- CI for CC -----------------------------------------
# Delta-method band for the complete-case fit: gradient evaluated at the CC
# model's OWN coefficients and combined with vcov(mod.cc), matching the point
# estimate (prob = predict(mod.cc)). Evaluating at the two-stage gamma/eta
# values instead would make the SE incoherent with the curve.
v.g <- vcov(mod.cc)
coef.cc <- coef_int_first(mod.cc)
cc.df.plot <- grid.data.cc[, c("age", "sex", "sf36_RP", "prob")]
cc.df.plot$SE <- cc.df.plot$asymp.LCL <- cc.df.plot$asymp.UCL <- NA_real_

cc.covar.names <- c(names(fit$par.list$c), names(fit$par.list$e))[-(1:(K - 1))]
for (i in seq_len(nrow(grid.data.cc))) {
  ki   <- as.integer(grid.data.cc$sf36_RP[i])
  xi   <- covar_row(cc.covar.names, "V.(gamma|eta).", grid.data.cc[i, , drop = FALSE])
  gr.g <- grad_naive_cc(coef.cc, x = xi, k = ki)
  se   <- sqrt(as.numeric(t(gr.g) %*% v.g %*% gr.g))
  cc.df.plot$SE[i]        <- se
  cc.df.plot$asymp.LCL[i] <- pmax(0, cc.df.plot$prob[i] - qnorm(0.975) * se)
  cc.df.plot$asymp.UCL[i] <- pmin(1, cc.df.plot$prob[i] + qnorm(0.975) * se)
}
colnames(cc.df.plot)[4] <- "prob"
cc.df.plot$df     <- Inf
cc.df.plot$Method <- "CC"
cc.df.plot <- cc.df.plot[cc.df.plot$age >= AGE_THR_FIXED, ]

# ----------------------- CI for 2-stage cured ------------------------------
sub.cured <- grid.data[grid.data$THR == 0, ]
cured.df.plot <- do.call(rbind, lapply(seq_len(nrow(sub.cured)), function(i) {
  row     <- sub.cured[i, , drop = FALSE]
  g.alpha <- grad(grad.cured.2stage, x = fit$par.list$a, data = row, k = K, resp = resp)
  se      <- sqrt(as.numeric(t(g.alpha) %*% V.alpha %*% g.alpha))
  data.frame(age = row$age, sex = row$sex, sf36_RP = row$sf36_RP, prob = row$prob,
             SE = se, df = Inf,
             asymp.LCL = pmax(0, row$prob - qnorm(0.975) * se),
             asymp.UCL = pmin(1, row$prob + qnorm(0.975) * se),
             Method = "2stage Cured")
}))

# ----------------------- CI for 2-stage uncured ----------------------------
params <- c(fit$par.list$b, fit$par.list$e, fit$par.list$c, fit$par.list$e)

ci_uncured <- function(sub, which.block) {
  do.call(rbind, lapply(seq_len(nrow(sub)), function(i) {
    row    <- sub[i, , drop = FALSE]
    g.full <- grad(grad.uncured.2stage, x = params, data.m = row, k = K,
                   nbe = nbe, ngce = ngce, Rvar = Rvar, Tauv = Tauv, resp = resp)
    if (which.block == "beta") {
      g <- g.full[1:nbe];                  Vsub <- V.beta
    } else {
      g <- g.full[(nbe + 1):(nbe + ngce)]; Vsub <- V.gamma
    }
    se <- sqrt(as.numeric(t(g) %*% Vsub %*% g))
    data.frame(age = row$age, sex = row$sex, sf36_RP = row$sf36_RP, prob = row$prob,
               SE = se, df = Inf,
               asymp.LCL = pmax(0, row$prob - qnorm(0.975) * se),
               asymp.UCL = pmin(1, row$prob + qnorm(0.975) * se),
               Method = "2stage Uncured")
  }))
}

uncured.df.plot.1 <- ci_uncured(grid.data[grid.data$THR == 1 & grid.data$age <  grid.data$age.THR, ], "beta")
uncured.df.plot.2 <- ci_uncured(grid.data[grid.data$THR == 1 & grid.data$age >= grid.data$age.THR, ], "gamma")

# ----------------------- Plot ----------------------------------------------

plot.df <- rbind(naive.df.plot, cc.df.plot, cured.df.plot, uncured.df.plot.1, uncured.df.plot.2)
plot.df$sf36_RP <- factor(plot.df$sf36_RP, levels = seq_along(sf36_levels), labels = sf36_levels)
plot.df$sex     <- factor(plot.df$sex, levels = c(0, 1), labels = c("Male", "Female"))

method_cols <- c("2stage Uncured" = "red", "2stage Cured" = "blue",
                 "CC" = "green4", "Naive" = "purple")
method_lty  <- c("2stage Uncured" = "solid", "2stage Cured" = "dashed",
                 "CC" = "longdash", "Naive" = "dotdash")

ggplot(plot.df, aes(age, prob, color = Method, linetype = Method)) +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL, fill = Method),
              alpha = 0.15, color = NA) +
  geom_vline(xintercept = AGE_THR_FIXED, linetype = "22", linewidth = 0.6, colour = "grey30") +
  geom_line(linewidth = 1.1) +
  geom_line(aes(y = asymp.LCL), linewidth = 0.37, alpha = 0.4) +
  geom_line(aes(y = asymp.UCL), linewidth = 0.37, alpha = 0.4) +
  scale_x_continuous(breaks = c(20, 30, 35, 40, 50, 60, 70), minor_breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  facet_grid(sex ~ sf36_RP, labeller = \(x) label_both(x, sep = " = ")) +
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  scale_linetype_manual(values = method_lty) +
  labs(x = "Age", y = "Estimated probability",
       colour = "Method:", fill = "Method:", linetype = "Method:",
       title = "Estimated Outcome Probabilities Across Age - Adult Perthes' Study",
       subtitle = paste0(if (OUTCOME_MODEL == "ACAT") "Adjacent categories model" else "Proportional odds model",
                         ", age-at-THR fixed at ", AGE_THR_FIXED,
                         ", remaining covariates at reference value 0")) +
  theme_bw(base_size = 18) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text    = element_text(size = 18, face = "bold"),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(face = "bold"),
    legend.key.width = unit(1.6, "cm"),
    panel.spacing    = unit(1, "cm"),
    axis.line        = element_line(linewidth = 0.4)
  )


# ----------------------- Save the plot -------------------------

#file_name <- sprintf(
#  "%s.est.outcome.probs(age.THR%s).png",
#  OUTCOME_MODEL,
#  AGE_THR_FIXED
#)

#ggsave(
#  filename = file_name,
#  width = 46,
#  height = 28,
#  units = "cm",
#  dpi = 600,
#  bg = "white"
#)
