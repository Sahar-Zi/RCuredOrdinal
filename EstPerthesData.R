# ----
library(tidyverse)
library(ggplot2)
library(MASS)
library(VGAM)
library(survival)
library(eha) # for phreg
library(smcure)
library(resample)
library(dplyr)
library(numDeriv)

source("gen.probs.R")
source("psudeo.likelihood.R")
source("analysis.R")
source("variance.estimation.R")

get.probs <- function(par.list, data, R, Tau, k, outcome.model=c("PO","ACAT")){
  
  # --------------------------------------------------
  # what if par is a list of vectors:
  #
  # a: {V.alpha1  | V.alpha2  | V.alphaZ1 | V.alphaZ2 | V.alphaR}
  # b: {V.beta1   | V.beta2   | V.betaZ1}
  # c: {V.gamma1  | V.gamma2  | V.gammaZ1 | V.gammaZ2 | V.gammaR | V.gammaT}
  # --------------------------------------------------
  a.inter <- par.list$a[1:(k-1)]
  a.cvars <- par.list$a[k:length(par.list$a)]
  b.inter <- par.list$b[1:(k-1)]
  b.cvars <- par.list$b[k:length(par.list$b)]
  if ("e" %in% names(par.list)) b.cvars <- c(b.cvars, par.list$e)
  c.inter <- par.list$c[1:(k-1)]
  c.cvars <- par.list$c[k:length(par.list$c)]
  if ("e" %in% names(par.list)) c.cvars <- c(c.cvars, par.list$e)
  
  mm.a <- model.matrix(as.formula(paste("~",paste(sub("V.alpha.","",names(a.cvars)),collapse="+"))), data=data)
  mm.b <- model.matrix(as.formula(paste("~",paste(sub("V.(beta|gamma).","",names(b.cvars)),collapse="+"))), data=data)
  mm.c <- model.matrix(as.formula(paste("~",paste(sub("V.gamma.","",names(c.cvars)),collapse="+"))), data=data)
  print(mm.c)
  xi.a <- exp(mm.a %*% t(cbind(a.inter, matrix(a.cvars, nrow = k-1, ncol = length(a.cvars), byrow=TRUE))))
  xi.b <- exp(mm.b %*% t(cbind(b.inter, matrix(b.cvars, nrow = k-1, ncol = length(b.cvars), byrow=TRUE))))
  xi.c <- exp(mm.c %*% t(cbind(c.inter, matrix(c.cvars, nrow = k-1, ncol = length(c.cvars), byrow=TRUE))))
  
  condition1 <- (data[[R]]<=data[[Tau]])
  
  if(outcome.model=="PO"){
    # for parameters alpha: 
    # => D = 0, delta = 0
    base.g.D0 <- 1/(1+xi.a)
    probs.0 <- cbind(1,base.g.D0)-cbind(base.g.D0,0)
    
    # for parameters beta/gamma:
    # => D = 1, delta = 1
    base.g.D1 <- 1/(1+xi.c)
    
    # => D = 1, delta = 0
    base.g.D1[condition1,] <- 1/(1+xi.b)[condition1,]
    probs.1 <- cbind(1,base.g.D1)-cbind(base.g.D1,0)
  }
  if(outcome.model=="ACAT"){
    xi.a <- cbind(1,xi.a)
    xi.b <- cbind(1,xi.b)
    xi.c <- cbind(1,xi.c)
    g.D0 <- c.xi.b <- g.D1 <- matrix(1, nrow = nrow(data), ncol = k)
    for(j in 2:k){
      g.D0[,j] <- xi.a[,j] * xi.a[,j-1]
      c.xi.b[,j] <- xi.b[,j] * xi.b[,j-1]
      g.D1[,j] <- xi.c[,j] * xi.c[,j-1]
    }
    # for parameters alpha: 
    probs.0 <- g.D0/rowSums(g.D0)
    
    # for parameters beta/gamma:
    g.D1[condition1,] <- c.xi.b[condition1,]
    probs.1 <- g.D1/rowSums(g.D1)
  }
  return(list("D0" = unname(probs.0), "D1" = unname(probs.1)))
}

data <- read.table(file = "data-2stage-April26.txt", sep = "\t", na = "", header = TRUE)
data$age.THR <- pmin(data$age.THR, data$age, na.rm=TRUE)

# modify the data in 32 observations where age.THR=age, shift their age.THR earlier within the same year:
set.seed(17)
data.eps <- data
num0 <- sum(data.eps$THR==1 & data.eps$age.THR == data.eps$age)
data.eps[data.eps$THR==1 & data.eps$age.THR == data.eps$age,"age.THR"] <- 
  data.eps[data.eps$THR==1 & data.eps$age.THR == data.eps$age,"age.THR"] - runif(num0, 0.1, 0.9)

#plot(data.eps$age, data.eps$HOOS_DL)
#plot(data.eps$age, data.eps$sf36_PF)

# there is a problem with the smcure function.
# we need, first, to choose the variables we use in a separate data set and pass it to smcure!!!
ddd <- data.eps[, c("age.THR", "age","THR", "obesity", "sex", "unilateral", "dysplasia",
                    "diag_age6_7", "diag_age8_11", "country.ind", "diag_age11", 
                    "chPtrt.surgery", "chPtrt.bracing", "chPtrt.casting",
                    "chPtrt.weight.bearing", "chPtrt.activity.restrict",
                    "chPtrt.phys.therapy", "chPtrt.walking.device",
                    "chPtrt.not.treated", "sf36_RP")]
# current obesity as a proxy for childhood obesity.....
# I will not use obesity!!!!!!
# "diag_ageUNK", removed

ddd$sf36_RP <- fct_collapse(
  ddd$sf36_RP,
  "25-75" = c(levels(ddd$sf36_RP)[2], levels(ddd$sf36_RP)[3], levels(ddd$sf36_RP)[4])
)

# Make sure it's still ordered
ddd$sf36_RP <- ordered(ddd$sf36_RP, levels = unique(ddd$sf36_RP))

# ddd.female <- ddd[ddd$sex==1,c("age.THR", "age", "THR", 
#                              "unilateral","dysplasia", "diag_age6_7", "diag_age8_11", "diag_age11",
#                              "chPtrt.surgery", "chPtrt.bracing", "chPtrt.casting",
#                              "chPtrt.weight.bearing", "chPtrt.activity.restrict",
#                              "chPtrt.phys.therapy", "chPtrt.walking.device",
#                              "chPtrt.not.treated", "sf36_RP")]
# ddd.female$sf36_RP


cureform <- as.formula(THR~sex)
survform <- as.formula(age.THR~sex+unilateral+diag_age6_7+diag_age8_11+diag_age11+
                         chPtrt.surgery+chPtrt.activity.restrict)
formula.c <- as.formula(sf36_RP~age+age.THR)
formula.a <- as.formula(sf36_RP~age)
equal.effect <- as.formula(~age+diag_age6_7)

mod <- ordcuredfit.param(survform, cureform, formula.c, formula.a, equal.effect,
                  Tau="age.THR", R="age", delta="THR", ddd, 
                  outcome.model="PO", var = F)

mod$value
s.t <- data.frame(est = unlist(mod$par),
              SD.sandwich = sqrt(diag(mod$variance$stacked.v.est)),
              Z = abs(unlist(mod$par)/sqrt(diag(mod$variance$stacked.v.est))))
#              naive.cc.est = NA, naive.cc.SD = NA)

#s.t[10:14,4:5] <- mod$naive
#s.t[20:26,4:5] <- mod$cc

s.t

summary(ddd.female)

ddd.female$sf36_RP <- ordered(ddd.female$sf36_RP)
# Proposed method outcome probabilities

sf36_levels <- c(0, "25-75", 100)
pred.df1 <- expand.grid(
  age = seq(from = 18, to = 76, by = 0.2),
  THR = c(0, 1),
  sf36_RP_V = c(0, "25-75", 100)
) %>%
  arrange(THR, age) %>%
  mutate(
    age.THR = 35,  # or consider updating this based on `THR`
    diag_age6_7 = 1,
    dysplasia = 1,
    chPtrt.activity.restrict = 0,
    sf36_RP = match(sf36_RP_V, sf36_levels)
  )

pred.df1$sf36_RP <- ordered(pred.df1$sf36_RP)

v.probs.est <- get.probs(mod$par, pred.df1, R="age", Tau="age.THR", k=3, outcome.model="PO")

pred.df1$phi <- NA
pred.df1$Method <- "2stage Proposed"

pred.df1$phi[pred.df1$THR==1] <- v.probs.est$D1[cbind(which(pred.df1$THR==1), pred.df1$sf36_RP[pred.df1$THR==1])]

pred.df1$phi[pred.df1$THR==0] <- v.probs.est$D0[cbind(which(pred.df1$THR==0), pred.df1$sf36_RP[pred.df1$THR==0])]
pred.df1[pred.df1$THR==0,]
# predicted.df.est$naive.pred <- NA
# predicted.df.est$cc.pred <- NA

## Naive
pred_naive <- pred.df1 %>% filter(THR == 0) 
mod.naive <- polr(formula.a, ddd.female %>% filter(THR == 0))
naive_probs <- predict(mod.naive, newdata = pred_naive, type = "probs")
pred_naive$phi <- naive_probs[cbind(1:nrow(pred_naive), pred_naive$sf36_RP)]
pred_naive$Method <- "Naive"

## CC
pred_cc <- pred.df1 %>% filter(THR == 1) 
mod.cc <- polr(formula.c, ddd.female %>% filter(THR == 1))
cc_probs <- predict(mod.cc, newdata = pred_cc, type = "probs")
pred_cc$phi <- cc_probs[cbind(1:nrow(pred_cc), pred_cc$sf36_RP)]
pred_cc$Method <- "Complete Case"

combined_df <- rbind(pred.df1, pred_naive, pred_cc)


## ---- Cumulative probabilities plot
# Compute cumulative probabilities
df_plot <- combined_df %>%
  filter(!(Method == "Complete Case" & age < 35)) %>%
  group_by(Method, THR, age) %>%
  arrange(sf36_RP, .by_group = TRUE) %>%
  mutate(cum_prob = cumsum(phi)) %>%
  ungroup()

# Facet labels
facet_labels <- c("0" = "Cured", "1" = "Uncured")

# Plot
ggplot(df_plot, aes(x = age, y = cum_prob, color = factor(sf36_RP), linetype = Method)) +
  geom_line(size = 1) +
  geom_vline(data = df_plot %>% filter(THR == 1),
             aes(xintercept = 35),
             linetype = "dashed", color = "black") +
  scale_x_continuous(
    breaks = c(20, 30, 35, 40, 50, 60, 70),
    minor_breaks = c(25,35,45,55,65)
  ) +
  #scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ THR, ncol = 1, labeller = labeller(THR = facet_labels)) +
  labs(
    x = "Age",
    y = expression("Cumulative Probability of  " * phi(V)),
    color = "sf36_RP Score",
    linetype = "Method",
    title = "Cumulative Probabilities by Age and THR Status",
    subtitle = expression("Predicted from 2stage Proposed, Naive, and Complete Case models")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12),
    legend.position = "top")

## ---- probability plot

# Filter for Complete Case age >= 35
probabilities_df <- combined_df %>%
  mutate(Method = case_when(
    Method == "2stage Proposed" & THR == 0 ~ "2stage Proposed D=0",
    Method == "2stage Proposed" & THR == 1 ~ "2stage Proposed D=1",
    TRUE ~ Method  # Leave other method names unchanged
  ))

df_plot <- probabilities_df %>%
  filter(!(Method == "Complete Case" & age < 35))

# Plot individual probabilities
ggplot(df_plot, aes(x = age, y = phi, color = Method)) +
  geom_line(size = 1) +
  geom_vline(aes(xintercept = 35), linetype = "dashed", color = "black") +
  scale_x_continuous(
    breaks = c(20, 30, 35, 40, 50, 60, 70),
    minor_breaks = c(25,35,45,55,65)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ sf36_RP, ncol = 5, labeller = label_both) +
  labs(
    x = "Age",
    y = expression("Probability of " * phi(V)),
    color = "Method",
    linetype = "THR",
    title = "Predicted Category Probabilities by Age",
    subtitle = "Faceted by sf36_RP level; shows both THR = 0 and THR = 1"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12),
    legend.position = "top"
  )




