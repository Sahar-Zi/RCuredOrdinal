# Install relevant packages:
library(MASS)
library(VGAM)
library(survival)
library(eha) # for phreg
library(smcure)
library(resample)
library(dplyr)
library(numDeriv)

# Check if the files are in the same directory
source("gen.probs.R")
source("gen.data.R")
source("psudeo.likelihood.R")
source("analysis.R")
source("variance.estimation.R")
source("run.sim.R")

# D parameters: beta_Z0, beta_Z1, beta_Z2
D.b <- c(0.5,0,1.5)

# T parameters: Scale, Shape, beta_Z1, beta_Z2
T.scale <- 1
T.shape <- 4.5
T.b <- c(-0.55,0)

# True V parameters
V.alpha0 <- c("V.alpha.1" = -1, "V.alpha.2" = 1)
V.beta0 <-  c("V.beta.1" = 0 , "V.beta.2" = 1.5)
V.gamma0 <- c("V.gamma.1" = -2, "V.gamma.2" = 2)
V.alphaZ <- c("V.alpha.Z1" = 0 , "V.alpha.Z2" = 2)
V.betaZ <- c("V.beta.Z1" = -2 , "V.beta.Z2" = 0)
V.gammaZ <- c("V.gamma.Z1" = -2 , "V.gamma.Z2" = 0)
V.etaZ <- c("V.eta.Z1" = -2 , "V.eta.Z2" = 0)
V.alpha.R <- -1.5
V.beta.R <- 1
V.gamma.R <- 1
V.gamma.T <- -2
V.gamma.TR <- 1.5

true.params <- list(d = c(D.b),
                    Tau = c(T.b),
                    shape = T.shape,
                    scale = T.scale,
                    a = c(V.alpha0, "V.alpha.R" = V.alpha.R, V.alphaZ),
                    e = c("V.eta.R" = V.beta.R, V.etaZ),
                    b = c(V.beta0),
                    c = c(V.gamma0, "V.gamma.Tau" = V.gamma.T, "V.gamma.R:Tau" = V.gamma.TR))

# Replications, Bootstrap simulations and sample size
seed <- 2*2213208
reps <- 400
#B <- 100
n <- 2000

# Outcome model
v.model <- "ACAT"        # Choose between "PO" / "ACAT"

cureform <- as.formula(delta~Z1+Z2)
survform <- as.formula(Tau~Z1+Z2)
formula.a <- as.formula(V~R+Z1+Z2)
formula.e <- as.formula(V~R+Z1+Z2)
formula.b <- as.formula(V~1)
formula.c <- as.formula(V~Tau+R:Tau)

sim <- run.sim.param(seed = seed, n = n, k = 3, par = true.params, outcome.model = v.model,
                     survform = survform, cureform = cureform, formula.a = formula.a, formula.e = formula.e,
                     formula.b = formula.b, formula.c = formula.c, equal.effect = equal.effect,
                     Tau="Tau", R="R", delta="delta", var = T,
                     replications = 400, alpha = 0.05, save.iter = T)
