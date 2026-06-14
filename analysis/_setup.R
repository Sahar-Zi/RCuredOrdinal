############################################################
## analysis/_setup.R
## Setup for the ordcure simulation & application scripts.
## Sourced by main_simulation.R (and the application script).
##
## This is ANALYSIS scaffolding, not package code: it lives in
## analysis/ and is excluded from the build via .Rbuildignore.
############################################################

## ---- clean session --------------------------------------

rm(list = ls())
gc()

## ---- reproducibility ------------------------------------

seed <- 2 * 2213208
set.seed(seed)

## ---- load the model code --------------------------------
## DEVELOPMENT (current): source the package functions directly.
## Source order is irrelevant (R resolves functions at call time),
## but core functions are listed first for readability.
##
## Once the package is built, replace this whole block with one line:
##   library(ordcure)            # installed package
##   # or, during development:  devtools::load_all(".")

core_files <- c(
  "outcome_probabilities.R",  # shared PO/ACAT category probabilities
  "pseudo_likelihood.R",      # stage-2 objective
  "cure_model.R",             # stage-1 Weibull mixture-cure EM
  "initialization.R",         # naive + complete-case starting values
  "variance.R",               # sandwich variance
  "ordcure.R",                # exported driver
  "gen_data.R",               # simulation data generator
  "methods.R"
)
invisible(lapply(file.path("R", core_files), source))

## simulation drivers (analysis side)
sim_files <- c(
  "simulation_utils.R",
  "run_simulation.R",
  "plots.R"
)
invisible(lapply(file.path("analysis", sim_files), source))

## ---- packages -------------------------------------------
## NOTE: MASS, VGAM, eha and numDeriv are referenced with :: in the
## model code, so they only need to be INSTALLED, not attached here.
## survival is attached because the stage-1 latency initialization
## builds a formula that calls Surv() unqualified.

library(survival)

## Analysis-only packages, used directly (unqualified) by the
## simulation / application scripts rather than by the model code.
## Prune any you don't actually use.
library(matrixStats)   # colSds() in simulation_utils.R
library(dplyr)
library(forcats)
library(emmeans)
library(ggplot2)
library(smcure)