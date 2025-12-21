############################################################
## Global setup
## Ordinal Cure Model with Censored Covariates
############################################################

## ---- Clean session --------------------------------------

rm(list = ls())
gc()

## ---- Reproducibility ------------------------------------

seed <- 2 * 2213208
set.seed(seed)

## ---- Packages -------------------------------------------

library(MASS)
library(VGAM)
library(survival)
library(eha)
library(smcure)
library(dplyr)
library(numDeriv)
library(matrixStats)

## ---- Source core functions -------------------------------

source("R//gen.probs.R")
source("R//gen.data.R")
source("R//pseudo.likelihood.R")
source("R//variance.estimation.R")
source("R//analysis.R")
source("R//simulation_utils.R")
source("R//run.sim.R")
