############################################################
## Main simulation study
## Ordinal Cure Model with Censored Covariates in the Presence of a Cured Fraction
############################################################

## ---- setup ----------------------------------------------

source("analysis/_setup.R")

## ---- true parameters ------------------------------------

## D (incidence / cure) parameters
D.b <- c("b.0" = 0.5, "b.Z1" = 0, "b.Z2" = 1.5)

## T (latency) parameters
T.shape <- 4.5
T.scale <- 1
T.b     <- c("zeta.Z1" = -0.55, "zeta.Z2" = 0)

## V (ordinal outcome) parameters
V.alpha0 <- c("V.alpha.1" = -1, "V.alpha.2" = 1)
V.beta0  <- c("V.beta.1"  = 0,  "V.beta.2"  = 1.5)
V.gamma0 <- c("V.gamma.1" = -2, "V.gamma.2" = 2)

V.alphaZ <- c("V.alpha.Z1" = 0,  "V.alpha.Z2" = 2)
V.betaZ  <- c("V.beta.Z1"  = -2, "V.beta.Z2"  = 0)
V.gammaZ <- c("V.gamma.Z1" = -2, "V.gamma.Z2" = 0)
V.etaZ   <- c("V.eta.Z1"   = -2, "V.eta.Z2"   = 0)

V.alpha.R  <- -1.5
V.beta.R   <-  1
V.eta.R    <-  1
V.gamma.R  <-  1
V.gamma.T  <- -2
V.gamma.TR <-  1.5

true.params <- list(
  d         = D.b,
  Tau       = T.b,
  Tau.shape = T.shape,
  Tau.scale = T.scale,
  a = c(V.alpha0, "V.alpha.R" = V.alpha.R, V.alphaZ),
  e = c("V.eta.R" = V.eta.R, V.etaZ),
  b = c(V.beta0),
  c = c(V.gamma0,
        "V.gamma.Tau"   = V.gamma.T,
        "V.gamma.R:Tau" = V.gamma.TR)
)

## ---- simulation settings --------------------------------

n            <- 2000
replications <- 400
alpha        <- 0.05
outcome.model <- "PO"   # "PO" or "ACAT"

## model formulas
cureform  <- delta ~ Z1 + Z2
survform  <- Tau   ~ Z1 + Z2
formula.a <- V ~ R + Z1 + Z2
formula.e <- V ~ R + Z1 + Z2
formula.b <- V ~ 1
formula.c <- V ~ Tau + R:Tau

## ---- run ------------------------------------------------

sim <- run_simulation(
  seed = seed,
  n = n,
  k = 3,
  par = true.params,
  outcome.model = outcome.model,
  survform = survform,
  cureform = cureform,
  formula.a = formula.a,
  formula.e = formula.e,
  formula.b = formula.b,
  formula.c = formula.c,
  Tau = "Tau",
  R = "R",
  delta = "delta",
  var = TRUE,
  replications = replications,
  alpha = alpha,
  verbose = TRUE
)

## ---- summaries ------------------------------------------

sim$summary.table.est

sim$summary.table.naive

sim$summary.table.cc

plot(sim, type = "bias", se = "sandwich") # bias + 95% CI (sandwich variance estimator), the default
plot(sim, type = "bias", se = "inv")      # bias + 95% CI (inverse-information variance estimator)
plot(sim, type = "coverage")              # empirical coverage

## ---- save (optional) ------------------------------------
# saveRDS(sim, file = "results/simulation_main.rds")