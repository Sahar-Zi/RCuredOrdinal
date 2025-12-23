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

source("R//outcome_probabilities.R")
source("R//gen_data.R")
source("R//ordinal_pseudo_likelihood.R")
source("R//variance.estimation.R")
source("R//ordcure.R")
source("R//simulation_utils.R")
source("R//run_simulation.R")
