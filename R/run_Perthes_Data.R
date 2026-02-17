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
## Fit model
## --------------------------------------------------

## ---- Global constants ----
OUTCOME_MODEL <- "ACAT"   # "PO" or "ACAT"
AGE_THR_FIXED <- 35

<<<<<<< HEAD
sf36_levels <- levels(X$sf36_RP)
=======
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
K <- n_distinct(X$sf36_RP)
cureform <- THR ~ sex
survform <- age.THR ~ sex + unilateral + diag_age6_7 + diag_age8_11 + 
  diag_age11 + chPtrt.surgery + chPtrt.activity.restrict + sex:unilateral
formula.a <- sf36_RP ~ age
formula.b <- sf36_RP ~ age
formula.c <- sf36_RP ~ age + age.THR
formula.e <- sf36_RP ~ sex + dysplasia + unilateral + chPtrt.surgery + chPtrt.activity.restrict + sex:age

cat("--------------------------------\n",OUTCOME_MODEL,"outcome regression model")

fit <- ordcure(
  cureform   = cureform,
  survform   = survform,
  formula.a  = formula.a,
  formula.e  = formula.e,
  formula.b  = formula.b,
  formula.c  = formula.c,
  Tau        = "age.THR",
  R          = "age",
  delta      = "THR",
  data       = X,
  outcome.model = OUTCOME_MODEL,
  var = TRUE
)

cat("negative log-loikelihood=",fit$opt$value,"\n")
print(data.frame("par.est" = unlist(fit$par.list), 
            "sandwich.se" = sqrt(diag(fit$variance$stacked.v.est)),
            "Z" = abs(unlist(fit$par.list)/sqrt(diag(fit$variance$stacked.v.est))),
            "p-value" = round(dnorm(abs(unlist(fit$par.list)/sqrt(diag(fit$variance$stacked.v.est)))),3),
            "95% CI" = paste("[",round(unlist(fit$par.list) - qnorm(0.975) * sqrt(diag(fit$variance$stacked.v.est)),3),", ", 
                         round(unlist(fit$par.list) + qnorm(0.975) * sqrt(diag(fit$variance$stacked.v.est)),3),"]",sep = ""),
            check.names = F))
             
<<<<<<< HEAD
# ----------------------- Create data grid for plots -------------------------
=======
print("--------------------------------")
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef

grid.data <- expand.grid(age = seq(from = 18, to = 76, by = 0.2),
                         THR = c(0, 1),
                         sex = c(0, 1),
                         sf36_RP_V = sf36_levels) %>% # Non-constant variables
  arrange(sex, THR, age) %>%
  mutate(
    age.THR = AGE_THR_FIXED, 
<<<<<<< HEAD
=======
    sex = 1,
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
    dysplasia = 0,
    unilateral = 0,
    chPtrt.activity.restrict = 0,
    chPtrt.surgery = 0,
    sf36_RP = match(sf36_RP_V, sf36_levels),
    delta = as.integer(THR == 1 & age > age.THR)) # Constant variables

grid.data$sf36_RP <- ordered(grid.data$sf36_RP) # Keep the outcome ordered

# ----------------------- Estimate 2stage outcome probabilities -------------------------

v.probs.est <- compute_outcome_probs( 
  par           = fit$par,
  data          = grid.data,
  delta         = "delta",
  k             = K,
<<<<<<< HEAD
  outcome.model = OUTCOME_MODEL)
=======
  outcome.model = OUTCOME_MODEL
)
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef

grid.data$phi <- NA
grid.data$phi[grid.data$THR==1] <- v.probs.est$D1[cbind(which(grid.data$THR==1), grid.data$sf36_RP[grid.data$THR==1])]
grid.data$phi[grid.data$THR==0] <- v.probs.est$D0[cbind(which(grid.data$THR==0), grid.data$sf36_RP[grid.data$THR==0])]

# ----------------------- Estimate naive outcome probabilities -------------------------

<<<<<<< HEAD
grid.data.naive <- grid.data %>% filter(THR == 0) 
if(OUTCOME_MODEL == "PO"){
  mod.naive <- polr(sf36_RP ~ age, data = X %>% filter(THR == 0))
  naive.probs <- predict(mod.naive, newdata = grid.data.naive, type = "probs")
=======
pred.df1$phi[pred.df1$THR==0] <- v.probs.est$D0[cbind(which(pred.df1$THR==0), pred.df1$sf36_RP[pred.df1$THR==0])]
pred.df1[pred.df1$THR==0,]
# predicted.df.est$naive.pred <- NA
# predicted.df.est$cc.pred <- NA

## Naive
pred_naive <- pred.df1 %>% filter(THR == 0) 
# mod.naive <- polr(sf36_RP ~ age, data = X %>% filter(THR == 0))
mod.naive <- vglm(sf36_RP ~ age, acat(reverse = FALSE, parallel = TRUE), data = X %>% filter(THR == 0))
naive_probs <- predict(mod.naive, newdata = pred_naive, type = "response")
pred_naive$phi <- naive_probs[cbind(1:nrow(pred_naive), pred_naive$sf36_RP)]
pred_naive$Method <- "Naive"

## CC
pred_cc <- pred.df1 %>% filter(THR == 1) 
#mod.cc <- polr(sf36_RP ~ age + age.THR + sex + dysplasia + unilateral + chPtrt.surgery + chPtrt.activity.restrict + sex:age, data = X %>% filter(THR == 1))
mod.cc <- vglm(sf36_RP ~ age + age.THR + sex + dysplasia + unilateral + chPtrt.surgery + chPtrt.activity.restrict + sex:age, acat(reverse = FALSE, parallel = TRUE), data = X %>% filter(THR == 1))
cc_probs <- predict(mod.cc, newdata = pred_cc, type = "response")
pred_cc$phi <- cc_probs[cbind(1:nrow(pred_cc), pred_cc$sf36_RP)]
pred_cc$Method <- "Complete Case"

combined_df <- rbind(pred.df1, pred_naive, pred_cc)


# ## ---- Cumulative probabilities plot
# # Compute cumulative probabilities
# df_plot <- combined_df %>%
#   filter(!(Method == "Complete Case" & age < 35)) %>%
#   group_by(Method, THR, age) %>%
#   arrange(sf36_RP, .by_group = TRUE) %>%
#   mutate(cum_prob = cumsum(phi)) %>%
#   ungroup()
# 
# # Facet labels
# facet_labels <- c("0" = "Cured", "1" = "Uncured")
# 
# # Plot
# ggplot(df_plot, aes(x = age, y = cum_prob, color = factor(sf36_RP), linetype = Method)) +
#   geom_line(linewidth = 1) +
#   geom_vline(data = df_plot %>% filter(THR == 1),
#              aes(xintercept = 35),
#              linetype = "dashed", color = "black") +
#   scale_x_continuous(
#     breaks = c(20, 30, 35, 40, 50, 60, 70),
#     minor_breaks = c(25,35,45,55,65)
#   ) +
#   #scale_y_continuous(limits = c(0, 1)) +
#   facet_wrap(~ THR, ncol = 1, labeller = labeller(THR = facet_labels)) +
#   labs(
#     x = "Age",
#     y = expression("Cumulative Probability of  " * phi(V)),
#     color = "sf36_RP Score",
#     linetype = "Method",
#     title = "Cumulative Probabilities by Age and THR Status",
#     subtitle = expression("Predicted from 2stage Proposed, Naive, and Complete Case models")
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     strip.text = element_text(size = 12),
#     legend.position = "top")
# 
# ## ---- probability plot
# 
# # Filter for Complete Case age >= 35
# probabilities_df <- combined_df %>%
#   mutate(Method = case_when(
#     Method == "2stage Proposed" & THR == 0 ~ "2stage Proposed D=0",
#     Method == "2stage Proposed" & THR == 1 ~ "2stage Proposed D=1",
#     TRUE ~ Method  # Leave other method names unchanged
#   ))
# 
# df_plot <- probabilities_df %>%
#   filter(!(Method == "Complete Case" & age < 35))
# 
# # Plot individual probabilities
# ggplot(df_plot, aes(x = age, y = phi, color = Method)) +
#   geom_line(linewidth = 1) +
#   geom_vline(aes(xintercept = 35), linetype = "dashed", color = "black") +
#   scale_x_continuous(
#     breaks = c(20, 30, 35, 40, 50, 60, 70),
#     minor_breaks = c(25,35,45,55,65)
#   ) +
#   scale_y_continuous(limits = c(0, 1)) +
#   facet_wrap(~ sf36_RP, ncol = 5, labeller = label_both) +
#   labs(
#     x = "Age",
#     y = expression("Probability of " * phi(V)),
#     color = "Method",
#     linetype = "THR",
#     title = "Predicted Category Probabilities by Age",
#     subtitle = "Faceted by sf36_RP level; shows both THR = 0 and THR = 1"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     strip.text = element_text(size = 12),
#     legend.position = "top"
#   )
# 

###### -----------------------------------

grad.cured.po <- function(par, data, k){
  a.cvars <- NULL
  a.inter <- par[1:(k-1)]
  a.cvars <- par[k:length(par)]
  mm.a <- model.matrix(as.formula(paste("~",paste(sub("V.alpha.","",names(a.cvars)),collapse="+"))), data=data)
  xi.a <- exp(mm.a %*% t(cbind(a.inter, matrix(a.cvars, nrow = k-1, ncol = length(a.cvars), byrow=TRUE))))
  base.g.D0 <- 1/(1+xi.a)
  probs.0 <- cbind(1,base.g.D0)-cbind(base.g.D0,0)
  outcome.probs <- probs.0[1,data$sf36_RP]
  return(outcome.probs)
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
}
if(OUTCOME_MODEL == "ACAT"){
  mod.naive <- vglm(sf36_RP ~ age, acat(reverse = FALSE, parallel = TRUE), data = X %>% filter(THR == 0))
  naive.probs <- predict(mod.naive, newdata = grid.data.naive, type = "response")
}
grid.data.naive$phi <- naive.probs[cbind(1:nrow(grid.data.naive), grid.data.naive$sf36_RP)]
grid.data.naive$Method <- "Naive"

# ----------------------- Estimate CC outcome probabilities -------------------------

grid.data.cc <- grid.data %>% filter(THR == 1)
if(OUTCOME_MODEL == "PO"){
  mod.cc <- polr(sf36_RP ~ age + age.THR + sex + dysplasia + unilateral + chPtrt.surgery + chPtrt.activity.restrict + sex:age, data = X %>% filter(THR == 1))
  cc.probs <- predict(mod.cc, newdata = grid.data.cc, type = "probs")
}
if(OUTCOME_MODEL == "ACAT"){
  mod.cc <- vglm(sf36_RP ~ age + age.THR + sex + dysplasia + unilateral + chPtrt.surgery + chPtrt.activity.restrict + sex:age, acat(reverse = FALSE, parallel = TRUE), data = X %>% filter(THR == 1))
  cc.probs <- predict(mod.cc, newdata = grid.data.cc, type = "response")
}
grid.data.cc$phi <- cc.probs[cbind(1:nrow(grid.data.cc), grid.data.cc$sf36_RP)]
grid.data.cc$Method <- "Complete Case"

# ----------------------- CI for naive -------------------------

grad_naive_cc <- function(alpha, x, k) {
  a1 <- alpha[1]
  a2 <- alpha[2]
  ax <- 0
  if(length(alpha)>k-1) ax <- alpha[3:length(alpha)]

  eta1 <- a1 + ax %*% x
  eta2 <- a2 + ax %*% x
  eta3 <- a1 + a2 + 2 * ax %*% x
  
  D <- 1 + exp(eta1) + exp(eta3)
  
  if(k==1){
<<<<<<< HEAD
    if(OUTCOME_MODEL == "ACAT"){
      grad <- c(
        d_1 = -(exp(eta1) + exp(eta3)) / D^2,
        d_2 = -exp(eta3) / D^2,
        d_x = -x * as.vector((exp(eta1) + 2 * exp(eta3)) / D^2))
    } else {
      grad <- c(1, 0, x) * as.vector(1/(1+exp(-eta1)))
    }
  
  }
  if(k==2){
    if(OUTCOME_MODEL == "ACAT"){
      grad <- c(
        d_1 = exp(eta1) / D^2,
        d_2 = -(exp(eta1) * exp(eta3)) / D^2,
        d_x = x * as.vector((exp(eta1) - exp(eta1) * exp(eta3)) / D^2))
    } else {
      grad <- c(
        d_1 = -exp(-eta1)/(1+exp(-eta1))^2,
        d_2 = -exp(-eta2)/(1+exp(-eta2))^2,
        d_x = -x * as.vector((exp(-eta2)/(1+exp(-eta2))^2) - (exp(-eta2)/(1+exp(-eta2))^2))
      )
    }
  }
  if(k==3){
    if(OUTCOME_MODEL == "ACAT"){
      grad <- c(
        d_1 = exp(eta3) / D^2,
        d_2 = (exp(eta3) + (exp(eta1) * exp(eta3))) / D^2,
        d_x = x * as.vector((exp(eta1) * exp(eta3) + 2 * exp(eta3)) / D^2))
    } else {
      grad <- c(1, 0, x) * as.vector(exp(-eta2)/(1+exp(-eta2))^2)
    }
=======
    grad <- c(
      d_alpha_1 = -(e1 + e2) / D^2,
      d_alpha_2 = -e2 / D^2,
      d_alpha_age = - (a * e1 + 2 * a * e2) / D^2)
  }
  if(k==2){
    grad <- c(
      d_alpha_1 = e1 / D^2,
      d_alpha_2 = -(e1 * e2) / D^2,
      d_alpha_age = a * (e1 - e1 * e2) / D^2)
  }
  if(k==3){
    grad <- c(
      d_alpha_1 = e2 / D^2,
      d_alpha_2 = (e2 + (e1 * e2)) / D^2,
      d_alpha_age = (a * e1 * e2 + 2 * a * e2) / D^2)
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
  }
  
  return(grad)
}

<<<<<<< HEAD
naive.df.plot <- grid.data.naive[,c("age","sex","sf36_RP","phi")]
naive.df.plot$df <- naive.df.plot$SE <- Inf
v.a <- vcov(mod.naive)

for (i in 1:nrow(grid.data.naive)) {
=======
grad_g <- function(gamma, age, age.THR, sex, sex.age, k) {
  a1 <- gamma[1]
  a2 <- gamma[2]
  b1 <- gamma[3]
  b2 <- gamma[4]
  b3 <- gamma[5]
  b4 <- gamma[6]
  
  eta1 <- a1 + b1 * age + b2 * age.THR + b3 * sex + b4 * sex.age
  eta2 <- a1 + a2 + 2 * (b1 * age + b2 * age.THR + b3 * sex + b4 * sex.age)
  
  e1 <- exp(eta1)
  e2 <- exp(eta2)
  
  D <- 1 + e1 + e2
  
  if(k==1){
    grad <- c(
      d_1 = -(e1 + e2) / D^2,
      d_2 = -e2 / D^2,
      d_age = -age * (e1 + 2 * e2) / D^2,
      d_age.THR = -age.THR * (e1 + 2 * e2) / D^2,
      d_sex = -sex * (e1 + 2 * e2) / D^2,
      d_sex.age = -sex.age * (e1 + 2 * e2) / D^2)
  }
  if(k==2){
    grad <- c(
      d_1 = e1 / D^2,
      d_2 = -(e1 * e2) / D^2,
      d_age = age * (e1 - e1 * e2) / D^2,
      d_age.THR = age.THR * (e1 - e1 * e2) / D^2,
      d_sex = sex * (e1 - e1 * e2) / D^2,
      d_sex.age = sex.age * (e1 - e1 * e2) / D^2)
  }
  if(k==3){
    grad <- c(
      d_1 = e2 / D^2,
      d_2 = (e2 + (e1 * e2)) / D^2,
      d_age = age * (e1 * e2 + 2 * e2) / D^2,
      d_age.THR = age.THR * (e1 * e2 + 2 * e2) / D^2,
      d_sex = sex * (e1 * e2 + 2 * e2) / D^2,
      d_sex.age = sex.age * (e1 * e2 + 2 * e2) / D^2)
  }
  
  return(grad)
}

v.a <- vcov(mod.naive)
for (i in 1:nrow(pred_naive)) {
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
  if (i %% 3 == 1) {
    gr.a <- grad_naive_cc(fit$par.list$a, x = grid.data.naive[i,c("age")], k = 1)
  } else if (i %% 3 == 2) {
    gr.a <- grad_naive_cc(fit$par.list$a, x = grid.data.naive[i,c("age")], k = 2)
  } else {
    gr.a <- grad_naive_cc(fit$par.list$a, x = grid.data.naive[i,c("age")], k = 3)
  }
  se <- sqrt(t(gr.a) %*% v.a %*% gr.a)
  naive.df.plot$SE[i] <- se
  naive.df.plot$asymp.LCL[i] <- pmax(0, naive.df.plot$phi[i] - qnorm(0.975) * se)
  naive.df.plot$asymp.UCL[i] <- pmin(1, naive.df.plot$phi[i] + qnorm(0.975) * se)
}
colnames(naive.df.plot) <- c("age", "sex", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
naive.df.plot$Method <- "Naive"

<<<<<<< HEAD
# ----------------------- CI for CC -------------------------

cc.df.plot <- grid.data.cc[,c("age","sex","sf36_RP","phi")]
cc.df.plot$df <- cc.df.plot$SE <- Inf
v.g <- vcov(mod.cc)

for (i in 1:nrow(grid.data.cc)) {
  if (i %% 3 == 1) {
    gr.g <- grad_naive_cc(c(fit$par$c, fit$par$e), x = c(grid.data.cc[i,"age"],                               # Age
                                                  AGE_THR_FIXED,                                       # Age.THR
                                                  grid.data.cc[i,"sex"], 0, 0, 0, 0,                   # Sex, dysp, uni, surg, act
                                                  grid.data.cc[i,"sex"]*grid.data.cc[i,"age"]), k = 1) # Age:Sex
  } else if (i %% 3 == 2) {
    gr.g <- grad_naive_cc(c(fit$par$c, fit$par$e), x = c(grid.data.cc[i,"age"],                               # Age
                                                  AGE_THR_FIXED,                                       # Age.THR
                                                  grid.data.cc[i,"sex"], 0, 0, 0, 0,                   # Sex, dysp, uni, surg, act
                                                  grid.data.cc[i,"sex"]*grid.data.cc[i,"age"]), k = 2) # Age:Sex
  } else {
    gr.g <- grad_naive_cc(c(fit$par$c, fit$par$e), x = c(grid.data.cc[i,"age"],                               # Age
                                                  AGE_THR_FIXED,                                       # Age.THR
                                                  grid.data.cc[i,"sex"], 0, 0, 0, 0,                   # Sex, dysp, uni, surg, act
                                                  grid.data.cc[i,"sex"]*grid.data.cc[i,"age"]), k = 3) # Age:Sex  
=======
v.g <- vcov(mod.cc)[c(1:5,10),c(1:5,10)]
for (i in 1:nrow(pred_cc)) {
  if (i %% 3 == 1) {
    gr.g <- grad_g(params[c(10:14,19)], age = pred_cc$age[i], age.THR = AGE_THR_FIXED, sex = 1, sex.age = pred_cc$age[i], k = 1)
  } else if (i %% 3 == 2) {
    gr.g <- grad_g(params[c(10:14,19)], age = pred_cc$age[i], age.THR = AGE_THR_FIXED, sex = 1, sex.age = pred_cc$age[i], k = 2)
  } else {
    gr.g <- grad_g(params[c(10:14,19)], age = pred_cc$age[i], age.THR = AGE_THR_FIXED,  sex = 1, sex.age = pred_cc$age[i], k = 3)
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
  }
  se <- sqrt(t(gr.g) %*% v.g %*% gr.g)
  cc.df.plot$SE[i] <- se
  cc.df.plot$asymp.LCL[i] <- pmax(0, cc.df.plot$phi[i] - qnorm(0.975) * se)
  cc.df.plot$asymp.UCL[i] <- pmin(1, cc.df.plot$phi[i] + qnorm(0.975) * se)
}
<<<<<<< HEAD
colnames(cc.df.plot) <- c("age", "sex", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
=======
colnames(cc.df.plot) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
cc.df.plot <- cc.df.plot %>% filter(age>=35)
# newdata <- expand.grid(
#   age = 35:76,
#   age.THR = 35,
#   sex = 1,
#   dysplasia = 0,
#   unilateral = 0,
#   chPtrt.surgery = 0,
#   chPtrt.activity.restrict = 0,
#   "age:sex" = 0
# )
# newdata$`age:sex` <- newdata$age
# 
# cc.df.plot <- acat_prob_ci_df(mod.cc, newdata)
#
# naive.df.plot <- as.data.frame(emmeans(mod.naive, ~ sf36_RP | age,
#                                        at = list(age = 18:76),
#                                        mode = "prob"))
# cc.df.plot <- as.data.frame(emmeans(mod.cc, ~ sf36_RP | age + age.THR + sex + dysplasia + unilateral + chPtrt.surgery + chPtrt.activity.restrict,
#                                     at = list(age = 35:76,
#                                               age.THR = 35,
#                                               sex = 1,
#                                               dysplasia = 0,
#                                               unilateral = 0,
#                                               chPtrt.surgery = 0,
#                                               chPtrt.activity.restrict = 0),
#                                     mode = "prob"))

naive.df.plot$Method <- "Naive"
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
cc.df.plot$Method <- "CC"
cc.df.plot <- cc.df.plot[cc.df.plot$age >= AGE_THR_FIXED,]

# ----------------------- CI for 2stage cured -------------------------

grad.cured.2stage <- function(par, data, k){
  a.cvars <- NULL
  a.inter <- par[1:(k-1)]
  a.cvars <- par[k:length(par)]
  mm.a <- model.matrix(as.formula(paste("~",paste(sub("V.alpha.","",names(a.cvars)),collapse="+"))), data=data)
  xi.a <- mm.a %*% t(cbind(a.inter, matrix(a.cvars, nrow = k-1, ncol = length(a.cvars), byrow=TRUE)))
  if(OUTCOME_MODEL == "PO"){
    base.g.D0 <- 1/(1+exp(xi.a))
    probs.0 <- cbind(1,base.g.D0)-cbind(base.g.D0,0)
  }
  if(OUTCOME_MODEL == "ACAT"){
    eta.a <- cbind(1, exp(xi.a))
    out <- matrix(1, nrow(data), k)
    out[, 2:k] <- eta.a[, 2:k] * eta.a[, 1:(k - 1)]
    probs.0 <- out / rowSums(out)
  }
  outcome.probs <- probs.0[1,data$sf36_RP]
  return(outcome.probs)
}

subset.cured <- grid.data[grid.data$THR == 0, -4]
variance.cured <- fit$variance$stacked.v.est[13:15,13:15] # The sub-matrix of the alpha parameters
cured.df.plot <- do.call(rbind, apply(subset.cured, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.alpha <- grad(grad.cured.2stage, x = fit$par$a, data = data_row, k=3)
  se <- sqrt(g.alpha %*% variance.cured %*% g.alpha)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age","sex","phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Cured")
})) %>% mutate(grid.data[grid.data$THR == 0,"sf36_RP"], .before = 3)
colnames(cured.df.plot) <- c("age", "sex", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")

# ----------------------- CI for 2stage uncured -------------------------

grad.uncured.2stage <- function(par, data.m, k) {
  
  ## ---- split parameters ---------------------------------------------------
  split_par <- function(par_block) {
    list(
      intercept = par_block[seq_len(k - 1)],
      covars    = if (length(par_block) > (k - 1)) par_block[-seq_len(k - 1)] else numeric(0)
    )
  }
  
  par.b <- split_par(par[1:9])
  par.c <- split_par(par[10:19])
  
  ## ---- helper: build model matrix ----------------------------------------
  build_mm <- function(cvars, data, prefix) {
    if (length(cvars) == 0) return(matrix(1, nrow(data), 1))
    
    vars <- sub(prefix, "", names(cvars))
    form <- reformulate(vars)
    model.matrix(form, data = data)
  }
  
  mm.b <- build_mm(par.b$covars, data.m, "V.(beta|eta).")
  mm.c <- build_mm(par.c$covars, data.m, "V.(gamma|eta).")
  
  ## ---- helper: compute linear predictors ---------------------------------
  compute_xi <- function(mm, intercept, covars) {
    coef_mat <- cbind(
      intercept,
      matrix(covars, nrow = k - 1, ncol = length(covars), byrow = TRUE)
    )
    mm %*% t(coef_mat)
  }
  
  xi.b <- compute_xi(mm.b, par.b$intercept, par.b$covars)
  xi.c <- compute_xi(mm.c, par.c$intercept, par.c$covars)
  
  ## ---- switch rule --------------------------------------------------------
  delta0 <- data.m$age <= data.m$age.THR
  
  ## ---- link functions -----------------------------------------------------
  if (OUTCOME_MODEL == "PO") {
    
    g <- plogis  # faster than manual logistic
    
    G <- g(xi.c)
    G[delta0, ] <- g(xi.b)[delta0, ]
    
    probs <- cbind(G, 1) - cbind(0, G)
    return(probs[cbind(seq_len(nrow(data.m)), data.m$sf36_RP)])
  }
  
  if (OUTCOME_MODEL == "ACAT") {
    
    g <- function(x) {
      ex <- exp(x)
      out <- cbind(1, ex)
      out[, 2:k] <- out[, 2:k] * out[, 1:(k - 1)]
      out / rowSums(out)
    }
    
    G <- g(xi.c)
    G[delta0, ] <- g(xi.b)[delta0, ]
    
    return(G[cbind(seq_len(nrow(data.m)), data.m$sf36_RP)])
  }
  
  stop("Unknown OUTCOME_MODEL")
}

params <- c(fit$par$b, fit$par$e, fit$par$c, fit$par$e)

subset.uncured.1 <- grid.data[grid.data$THR == 1 & grid.data$age < grid.data$age.THR, -4]
variance.uncured.beta <- fit$variance$stacked.v.est[c(22:24,16:21),c(22:24,16:21)] # The sub-matrix of the beta parameters
uncured.df.plot.1 <- do.call(rbind, apply(subset.uncured.1, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.beta <- grad(grad.uncured.2stage, x = params, data = data_row, k=3)[1:9]
  se <- sqrt(g.beta %*% variance.uncured.beta %*% g.beta)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age", "sex", "phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Uncured")
})) %>% mutate(grid.data[grid.data$THR == 1 & grid.data$age < grid.data$age.THR, "sf36_RP"], .before = 3)
colnames(uncured.df.plot.1) <- c("age", "sex", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")

<<<<<<< HEAD

subset.uncured.2 <- grid.data[grid.data$THR == 1 & grid.data$age >= grid.data$age.THR, -4]
variance.uncured.gamma <- fit$variance$stacked.v.est[c(25:28,16:21),c(25:28,16:21)] # The sub-matrix of the gamma parameters
=======
subset.uncured.2 <- pred.df1[pred.df1$THR == 1 & pred.df1$age > pred.df1$age.THR, -3]
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
uncured.df.plot.2 <- do.call(rbind, apply(subset.uncured.2, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.gamma <- grad(grad.uncured.2stage, x = params, data = data_row, k=3)[10:19]
  se <- sqrt(g.gamma %*% variance.uncured.gamma %*% g.gamma)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age", "sex", "phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Uncured")
})) %>% mutate(grid.data[grid.data$THR == 1 & grid.data$age >= grid.data$age.THR, "sf36_RP"], .before = 3)
colnames(uncured.df.plot.2) <- c("age", "sex", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")

# ----------------------- Plot -------------------------

plot.df <- rbind(naive.df.plot, cc.df.plot, cured.df.plot, uncured.df.plot.1, uncured.df.plot.2)
plot.df$sf36_RP <- factor(
  plot.df$sf36_RP,
  levels = c(1, 2, 3),
  labels = c(0, "25-75", 100)
)
plot.df$sex <- factor(
  plot.df$sex,
  levels = c(0, 1),
  labels = c("Male", "Female")
)

# define aesthetics
method_cols <- c(
  "2stage Uncured" = "red",
  "2stage Cured"   = "blue",
  "CC"             = "green4",
  "Naive"          = "purple"
)

method_lty <- c(
  "2stage Uncured" = "solid",
  "2stage Cured"   = "dashed",
  "CC"             = "longdash",
  "Naive"          = "dotdash"
)

ggplot(
  plot.df,
  aes(age, prob, color = Method, linetype = Method)) +
  # CI ribbon
  geom_ribbon(
    aes(ymin = asymp.LCL, ymax = asymp.UCL, fill = Method),
    alpha = 0.15,
    color = NA
    ) +
  # reference line
  geom_vline(xintercept = AGE_THR_FIXED,
             linetype = "22",
             linewidth = 0.6,
             colour = "grey30"
             ) +
  # main curve
  geom_line(linewidth = 1.1) +
  # CI boundaries
  geom_line(aes(y = asymp.LCL), linewidth = 0.37, alpha = 0.4) +
  geom_line(aes(y = asymp.UCL), linewidth = 0.37, alpha = 0.4) +
  # axes
  scale_x_continuous(
    breaks = c(20, 30, 35, 40, 50, 60, 70),
    minor_breaks = NULL
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  # facets
<<<<<<< HEAD
  facet_grid(
    sex ~ sf36_RP,
    labeller = \(x) label_both(x, sep = " = ")
  ) +
=======
  facet_wrap(~ sf36_RP, ncol = 3, labeller = \(x) label_both(x, sep = " = ")) +
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
  # manual scales (defined once above)
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  scale_linetype_manual(values = method_lty) +
  # labels
  labs(
    x = "Age",
    y = "Estimated probability",
    colour = "Method:",
    fill = "Method:",
    linetype = "Method:",
    title = paste(
      "Estimated", OUTCOME_MODEL,
      "Outcome Probabilities by Age\nWhen age of THR is",AGE_THR_FIXED
    )
  ) +
  # theme
<<<<<<< HEAD
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 18, face = "bold"),
=======
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 20, face = "bold"),
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.6, "cm"),
<<<<<<< HEAD
    panel.spacing = unit(1, "cm"),
    axis.line = element_line(linewidth = 0.4)
  )

## ----------------------- Save the plot -------------------------
# 
# file_name <- sprintf(
#   "%s.Pr.ageTHR%s.with.interaction.png",
#   OUTCOME_MODEL,
#   AGE_THR_FIXED
# )
# 
# ggsave(
#   filename = file_name,
#   width = 11.69,
#   height = 8.27,
#   units = "in",
#   dpi = 300,
#   bg = "white"
# )
=======
    axis.line = element_line(linewidth = 0.4)
  )

file_name <- sprintf(
  "%s.Pr.ageTHR%s.with.interaction.png",
  OUTCOME_MODEL,
  AGE_THR_FIXED
)

ggsave(
  filename = file_name,
  width = 11.69,
  height = 8.27,
  units = "in",
  dpi = 300,
  bg = "white"
)
>>>>>>> c08abff8f3d9e699dcf011b89b1a097c97d9eaef
