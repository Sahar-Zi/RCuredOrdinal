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
library(emmeans)

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

ddd$sf36_RP <- ordered(ddd$sf36_RP)

plot.df <- ddd %>%
  group_by(THR, sf36_RP) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(percent = 100 * count / sum(count))

ggplot(plot.df, aes(x = factor(sf36_RP), y = count, fill = factor(THR))) +
  ylim(c(0,400)) +
  geom_col(position = position_dodge(width = 0.35), width = 1.4) +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.35),
            vjust = 1.5, size = 3.5) +
  geom_text(aes(label = paste0(sprintf("%.1f", percent), "%")),
            position = position_dodge(width = 0.35),
            vjust = 3, size = 3.2, color = "black") +
  geom_vline(xintercept = c(1.5, 4.5), color = "red", linetype = "dashed", size = 1.25) +
  scale_fill_manual(
    values = c("0" = "#A9A9FF", "1" = "#FFD6A5"),
    name = expression(delta~": Total Hip Replacement"),
    labels = c("0: No", "1: Yes")
  ) +
  annotate("text", x = 3, y = max(plot.df$count) * 0.75,
           label = "These outcome categories were grouped",
           color = "red", size = 4, fontface = "bold")+
  labs(x = "SF-36 RP", y = "Count", title = expression("Number of Participants by"~delta)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5))

# current obesity as a proxy for childhood obesity.....
# I will not use obesity!!!!!!
# "diag_ageUNK", removed
ddd$sf36_RP <- fct_collapse(
  ddd$sf36_RP,
  "(25,50,75)" = c(levels(ddd$sf36_RP)[2], levels(ddd$sf36_RP)[3], levels(ddd$sf36_RP)[4])
)

# Make sure it's still ordered
ddd$sf36_RP <- ordered(ddd$sf36_RP, levels = unique(ddd$sf36_RP))
# ddd.female <- ddd[ddd$sex==1,c("age.THR", "age", "THR", 
#                                "unilateral","dysplasia", "diag_age6_7", "diag_age8_11", "diag_age11",
#                                "chPtrt.surgery", "chPtrt.bracing", "chPtrt.casting",
#                                "chPtrt.weight.bearing", "chPtrt.activity.restrict",
#                                "chPtrt.phys.therapy", "chPtrt.walking.device",
#                                "chPtrt.not.treated", "sf36_RP")]
# ddd.female$sf36_RP
cureform <- as.formula(THR~sex)
survform <- as.formula(age.THR~sex+unilateral+diag_age6_7+diag_age8_11+diag_age11+
                         chPtrt.surgery+chPtrt.activity.restrict)
formula.c <- as.formula(sf36_RP~age+age.THR)
formula.a <- as.formula(sf36_RP~age)
equal.effect <- as.formula(~1)

mod <- ordcuredfit.param(survform, cureform, formula.c, formula.a, equal.effect,
                         Tau="age.THR", R="age", delta="THR", ddd, 
                         outcome.model="PO", var = T)

s.t <- data.frame(est = unlist(mod$par),
                  SD.sandwich = sqrt(diag(mod$variance$stacked.v.est)),
                  Z = abs(unlist(mod$par)/sqrt(diag(mod$variance$stacked.v.est))))
#              SD.inv = sqrt(-diag(mod$variance$G.tilde.inv)),
#              naive.cc.est = NA, naive.cc.SD = NA)

#s.t[10:14,4:5] <- mod$naive
#s.t[20:26,4:5] <- mod$cc

s.t

summary(ddd)

ddd$sf36_RP <- ordered(ddd$sf36_RP)
# Proposed method outcome probabilities

sf36_levels <- c("0", "(25,50,75)", "100")
pred.df1 <- expand.grid(
  age = seq(from = 18, to = 76, by = 0.2),
  THR = c(0, 1),
  sf36_RP_V = c("0", "(25,50,75)", "100")
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
mod.naive <- polr(formula.a, ddd %>% filter(THR == 0))
naive_probs <- predict(mod.naive, newdata = pred_naive, type = "probs")
pred_naive$phi <- naive_probs[cbind(1:nrow(pred_naive), pred_naive$sf36_RP)]
pred_naive$Method <- "Naive"

## CC
pred_cc <- pred.df1 %>% filter(THR == 1) 
mod.cc <- polr(formula.c, ddd %>% filter(THR == 1))
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
  geom_line(linewidth = 1) +
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
  geom_line(linewidth = 1) +
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


###### -----------------------------------

grad.cured <- function(par, data, k){
  a.cvars <- NULL
  a.inter <- par[1:(k-1)]
  a.cvars <- par[k:length(par)]
  mm.a <- model.matrix(as.formula(paste("~",paste(sub("V.alpha.","",names(a.cvars)),collapse="+"))), data=data)
  xi.a <- exp(mm.a %*% t(cbind(a.inter, matrix(a.cvars, nrow = k-1, ncol = length(a.cvars), byrow=TRUE))))
  base.g.D0 <- 1/(1+xi.a)
  probs.0 <- cbind(1,base.g.D0)-cbind(base.g.D0,0)
  outcome.probs <- probs.0[1,data$sf36_RP]
  return(outcome.probs)
}

grad.uncured <- function(par, data, k){
  b.cvars <- NULL
  c.cvars <- NULL
  par.b <- par[1:3]
  par.c <- par[4:7]
  b.inter <- par.b[1:(k-1)]
  
  if(length(par.b)>=k)  b.cvars <- par.b[k:length(par.b)]
  
  c.inter <- par.c[1:(k-1)]
  c.cvars <- par.c[k:length(par.c)]
  
  mm.b <- model.matrix(as.formula(paste("~",paste(sub("V.(beta|gamma).","",names(b.cvars)),collapse="+"))), data=data)
  mm.c <- model.matrix(as.formula(paste("~",paste(sub("V.gamma.","",names(c.cvars)),collapse="+"))), data=data)
  
  xi.b <- exp(mm.b %*% t(cbind(b.inter, matrix(b.cvars, nrow = k-1, ncol = length(b.cvars), byrow=TRUE))))
  xi.c <- exp(mm.c %*% t(cbind(c.inter, matrix(c.cvars, nrow = k-1, ncol = length(c.cvars), byrow=TRUE))))
  
  condition1 <- (data["age"]<=data["age.THR"])
  # => D = 1, delta = 1
  base.g.D1 <- 1/(1+xi.c)
  
  # => D = 1, delta = 0
  base.g.D1[condition1,] <- 1/(1+xi.b)[condition1,]
  probs.1 <- cbind(1,base.g.D1)-cbind(base.g.D1,0)
  outcome.probs <- probs.1[1,data$sf36_RP]
  return(outcome.probs)
}

params <- c(mod$par$b, mod$par$e, mod$par$c, mod$par$e)

# keep age in the results
naive.df.plot <- as.data.frame(emmeans(mod.naive, ~ sf36_RP | age,
                                       at = list(age = 18:76),
                                       mode = "prob"))


cc.df.plot <- as.data.frame(emmeans(mod.cc, ~ sf36_RP | age,
                                    at = list(age = 35:76),
                                    mode = "prob"))

naive.df.plot$Method <- "Naive"
cc.df.plot$Method <- "CC"

subset.cured <- pred.df1[pred.df1$THR == 0, ]
cured.df.plot <- do.call(rbind, apply(subset.cured, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.alpha <- grad(grad.cured, x = mod$par$a, data = data_row, k=3)
  v.alpha <- mod$variance$stacked.v.est[12:14,12:14]
  se <- sqrt(g.alpha %*% v.alpha %*% g.alpha)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age","sf36_RP_V","phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Cured")
}))
colnames(cured.df.plot) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")
cured.df.plot <- cured.df.plot %>% mutate(sf36_RP = ifelse(is.na(sf36_RP), "(25,50,75)", sf36_RP))

subset.uncured.1 <- pred.df1[pred.df1$THR == 1 & pred.df1$age <= pred.df1$age.THR, ]
uncured.df.plot.1 <- do.call(rbind, apply(subset.uncured.1, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.beta <- grad(grad.uncured, x = params, data = data_row, k=3)[1:3]
  v.beta <- mod$variance$stacked.v.est[c(15:17),c(15:17)]
  se <- sqrt(g.beta %*% v.beta %*% g.beta)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age","sf36_RP_V","phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Uncured")
}))
colnames(uncured.df.plot.1) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")
uncured.df.plot.1 <- uncured.df.plot.1 %>% mutate(sf36_RP = ifelse(is.na(sf36_RP), "(25,50,75)", sf36_RP))

subset.uncured.2 <- pred.df1[pred.df1$THR == 1 & pred.df1$age > pred.df1$age.THR, ]
uncured.df.plot.2 <- do.call(rbind, apply(subset.uncured.2, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.gamma <- grad(grad.uncured, x = params, data = data_row, k=3)[4:7]
  v.gamma <- mod$variance$stacked.v.est[c(18:21),c(18:21)]
  se <- sqrt(g.gamma %*% v.gamma %*% g.gamma)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age","sf36_RP_V","phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Uncured")
}))
colnames(uncured.df.plot.2) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")
uncured.df.plot.2 <- uncured.df.plot.2 %>% mutate(sf36_RP = ifelse(is.na(sf36_RP), "(25,50,75)", sf36_RP))

plot.df <- rbind(naive.df.plot, cc.df.plot, cured.df.plot, uncured.df.plot.1, uncured.df.plot.2)


ggplot(plot.df, aes(x = age, y = prob, color = Method, fill = Method, group = Method, linetype = Method)) +
  # CI shading
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL, fill = Method), alpha = 0.1, color = NA) +
  # main line
  geom_line(aes(y = prob, color = Method, linetype = Method), size = 1) +
  # CI bounds
  geom_line(aes(y = asymp.LCL, color = Method, linetype = Method), size = 0.5) +
  geom_line(aes(y = asymp.UCL, color = Method, linetype = Method), size = 0.5) +
  geom_vline(aes(xintercept = 35), linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = c(20, 30, 35, 40, 50, 60, 70),
                     minor_breaks = c(25,35,45,55,65)) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ sf36_RP, ncol = 3, labeller = label_both) +
  labs(y = "Estimated probability", x = "Age", color = "Method", fill = "Method", linetype = "Method", title = "Estimated Outcome Probabilities by Age\n When age of THR is 35") +
  scale_color_manual(values = c(
    "2stage Uncured" = "red",
    "2stage Cured"  = "blue",
    "CC"            = "green4",
    "Naive"         = "purple"
  )) +
  scale_fill_manual(values = c(
    "2stage Uncured" = "red",
    "2stage Cured"  = "blue",
    "CC"            = "green4",
    "Naive"         = "purple"
  )) +
  scale_linetype_manual(values = c(
    "2stage Uncured" = "solid",
    "2stage Cured"  = "dashed",   # small dashes
    "CC"            = "longdash", # medium dashes
    "Naive"         = "dotdash"   # longer dashes
  )) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12),
    legend.position = "top"
  )




