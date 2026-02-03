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

X$sf36_RP <- fct_collapse(
  X$sf36_RP,
  "25-75" = c(levels(X$sf36_RP)[2], levels(X$sf36_RP)[3], levels(X$sf36_RP)[4])
)

## --------------------------------------------------
## Run model
## --------------------------------------------------

fit <- ordcure(
  cureform   = THR ~ sex,
  survform   = age.THR ~ sex + unilateral + diag_age6_7 + diag_age8_11 + 
    diag_age11 + chPtrt.surgery + chPtrt.activity.restrict + sex * unilateral,
  formula.a  = sf36_RP ~ age,
  formula.e  = sf36_RP ~ sex + dysplasia + unilateral + chPtrt.surgery + chPtrt.activity.restrict + sex:age,
  formula.b  = sf36_RP ~ age,
  formula.c  = sf36_RP ~ age + age.THR,
  Tau        = "age.THR",
  R          = "age",
  delta      = "THR",
  data       = X,
  outcome.model = "ACAT",
  var = TRUE
)

print(fit$opt$value)
print(cbind("par.est" = unlist(fit$par.list), 
            "sandwich.se" = sqrt(diag(fit$variance$stacked.v.est)),
            "Z" = abs(unlist(fit$par.list)/sqrt(diag(fit$variance$stacked.v.est)))))


sf36_levels <- levels(X$sf36_RP)
pred.df1 <- expand.grid(
  age = seq(from = 18, to = 76, by = 0.2),
  THR = c(0, 1),
  sf36_RP_V = sf36_levels
) %>%
  arrange(THR, age) %>%
  mutate(
    age.THR = 35,  # or consider updating this based on `THR`
    sex = 1,
    dysplasia = 0,
    unilateral = 0,
    chPtrt.activity.restrict = 0,
    chPtrt.surgery = 0,
    sf36_RP = match(sf36_RP_V, sf36_levels),
    delta = as.integer(THR == 1 & age > age.THR)
  )

pred.df1$sf36_RP <- ordered(pred.df1$sf36_RP)

v.probs.est <- compute_outcome_probs(
  par           = fit$par,
  data          = pred.df1,
  delta         = "delta",
  k             = 3,
  outcome.model = "ACAT"
)

pred.df1$phi <- NA
#pred.df1$Method <- "2stage Proposed"

pred.df1$phi[pred.df1$THR==1] <- v.probs.est$D1[cbind(which(pred.df1$THR==1), pred.df1$sf36_RP[pred.df1$THR==1])]

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
# mod.cc <- polr(sf36_RP ~ age + age.THR + sex + dysplasia + unilateral + chPtrt.surgery + chPtrt.activity.restrict + sex:age, data = X %>% filter(THR == 1))
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
}

grad.uncured.po <- function(par, data.m, k){
  b.cvars <- NULL
  c.cvars <- NULL
  par.b <- par[1:9]
  par.c <- par[10:19]
  b.inter <- par.b[1:(k-1)]
  
  if(length(par.b)>=k)  b.cvars <- par.b[k:length(par.b)]
  
  c.inter <- par.c[1:(k-1)]
  c.cvars <- par.c[k:length(par.c)]
  mm.b <- model.matrix(as.formula(paste("~",paste(sub("V.(beta|eta).","",names(b.cvars)),collapse="+"))), data=data.m)
  mm.c <- model.matrix(as.formula(paste("~",paste(sub("V.(gamma|eta).","",names(c.cvars)),collapse="+"))), data=data.m)

  xi.b <- mm.b %*% t(cbind(b.inter, matrix(b.cvars, nrow = k-1, ncol = length(b.cvars), byrow=TRUE)))
  xi.c <- mm.c %*% t(cbind(c.inter, matrix(c.cvars, nrow = k-1, ncol = length(c.cvars), byrow=TRUE)))
  
  delta0 <- (data.m["age"]<=data.m["age.THR"])
  
  g <- function(x) 1 / (1 + exp(-x))
  
  G1 <- g(xi.c)
  G1[delta0, ] <- g(xi.b)[delta0, ]
  
  probs.1 <- cbind(G1, 1) - cbind(0, G1)
  
  # # => D = 1, delta = 1
  # base.g.D1 <- 1/(1+xi.c)
  # 
  # # => D = 1, delta = 0
  # base.g.D1[condition1,] <- 1/(1+xi.b)[condition1,]
  # probs.1 <- cbind(1,base.g.D1)-cbind(base.g.D1,0)
  # outcome.probs <- probs.1[1,data$sf36_RP]
  return(probs.1[,data.m$sf36_RP])
}

grad.cured.acat <- function(par, data, k){
  a.cvars <- NULL
  a.inter <- par[1:(k-1)]
  a.cvars <- par[k:length(par)]
  mm.a <- model.matrix(as.formula(paste("~",paste(sub("V.alpha.","",names(a.cvars)),collapse="+"))), data=data)
  xi.a <- mm.a %*% t(cbind(a.inter, matrix(a.cvars, nrow = k-1, ncol = length(a.cvars), byrow=TRUE)))
  eta.a <- cbind(1, exp(xi.a))
  out <- matrix(1, nrow(data), k)
  out[, 2:k] <- eta.a[, 2:k] * eta.a[, 1:(k - 1)]
  probs.0 <- out / rowSums(out)
  outcome.probs <- probs.0[,data$sf36_RP]
  return(outcome.probs)
}

grad.uncured.acat <- function(par, data.m, k){
  b.cvars <- NULL
  c.cvars <- NULL
  par.b <- par[1:9]
  par.c <- par[10:19]
  b.inter <- par.b[1:(k-1)]
  
  if(length(par.b)>=k)  b.cvars <- par.b[k:length(par.b)]
  
  c.inter <- par.c[1:(k-1)]
  c.cvars <- par.c[k:length(par.c)]
  mm.b <- model.matrix(as.formula(paste("~",paste(sub("V.(beta|eta).","",names(b.cvars)),collapse="+"))), data=data.m)
  mm.c <- model.matrix(as.formula(paste("~",paste(sub("V.(gamma|eta).","",names(c.cvars)),collapse="+"))), data=data.m)
  
  xi.b <- mm.b %*% t(cbind(b.inter, matrix(b.cvars, nrow = k-1, ncol = length(b.cvars), byrow=TRUE)))
  xi.c <- mm.c %*% t(cbind(c.inter, matrix(c.cvars, nrow = k-1, ncol = length(c.cvars), byrow=TRUE)))
  
  delta0 <- (data.m["age"]<=data.m["age.THR"])
  
  g <- function(x) {
    x <- cbind(1, exp(x))
    out <- matrix(1, 1, k)
    out[, 2:k] <- x[, 2:k] * x[, 1:(k - 1)]
    out / rowSums(out)
  }
  
  G1 <- g(xi.c)
  G1[delta0, ] <- g(xi.b)[delta0, ]
  
  outcome.probs <- G1[,data.m$sf36_RP]
  return(outcome.probs)
}


params <- c(fit$par$b, fit$par$e, fit$par$c, fit$par$e)

naive.df.plot <- pred_naive[,c("age","sf36_RP","phi")]
naive.df.plot$df <- naive.df.plot$SE <- Inf

cc.df.plot <- pred_cc[,c("age","sf36_RP","phi")]
cc.df.plot$df <- cc.df.plot$SE <- Inf

grad_a <- function(alpha, a, k) {
  a1 <- alpha[1]
  a2 <- alpha[2]
  ax <- alpha[3]
  
  eta1 <- a1 + ax * a
  eta2 <- a1 + a2 + 2 * ax * a
  
  e1 <- exp(eta1)
  e2 <- exp(eta2)
  
  D <- 1 + e1 + e2
  
  if(k==1){
    grad <- c(
      d_alpha_1 = -(e1 + e2) / D^2,
      d_alpha_2 = -e2 / D^2,
      d_alpha_age = - a * (e1 + 2 * e2) / D^2)
  }
  if(k==2){
    grad <- c(
      d_alpha_1 = e1 / D^2,
      d_alpha_2 = -(e1*e2) / D^2,
      d_alpha_age = a * ((e1 / D) - (e1 * (e1 + 2 * e2) / D^2)))
  }
  if(k==3){
    grad <- c(
      d_alpha_1 = e2 / D^2,
      d_alpha_2 = (e2 / D^2) * (1 - (e2 / D^2)),
      d_alpha_age = a * ((2 * e2 / D) - (e2 * (e1 + 2 * e2) / D^2)))
  }
  
  return(grad)
}

grad_g <- function(alpha, x1, x2, k) {
  a1 <- alpha[1]
  a2 <- alpha[2]
  b1 <- alpha[3]
  b2 <- alpha[4]
  
  eta1 <- a1 + b1 * x1 + b2 * x2
  eta2 <- a1 + a2 + 2 * (b1 * x1 + b2 * x2)
  
  e1 <- exp(eta1)
  e2 <- exp(eta2)
  
  D <- 1 + e1 + e2
  
  if (k == 1) {
    grad <- c(
      d_alpha_1 = -(e1 + e2) / D^2,
      d_alpha_2 = -e2 / D^2,
      d_beta_1  = -x1 * (e1 + 2 * e2) / D^2,
      d_beta_2  = -x2 * (e1 + 2 * e2) / D^2
    )
  }
  
  if (k == 2) {
    grad <- c(
      d_alpha_1 = e1 / D^2,
      d_alpha_2 = -(e1 * e2) / D^2,
      d_beta_1  = x1 * ((e1 / D) - (e1 * (e1 + 2 * e2) / D^2)),
      d_beta_2  = x2 * ((e1 / D) - (e1 * (e1 + 2 * e2) / D^2))
    )
  }
  
  if (k == 3) {
    grad <- c(
      d_alpha_1 = e2 / D^2,
      d_alpha_2 = (e2 / D^2) * (1 - (e2 / D^2)),
      d_beta_1  = x1 * ((2 * e2 / D) - (e2 * (e1 + 2 * e2) / D^2)),
      d_beta_2  = x2 * ((2 * e2 / D) - (e2 * (e1 + 2 * e2) / D^2))
    )
  }
  
  return(grad)
}

v.a <- fit$variance$stacked.v.est[13:15,13:15]
for (i in 1:nrow(pred_naive)) {
  if (i %% 3 == 1) {
    gr.a <- grad_a(fit$par.list$a, a = pred_naive$age[i], k = 1)
  } else if (i %% 3 == 2) {
    gr.a <- grad_a(fit$par.list$a, a = pred_naive$age[i], k = 2)
  } else {
    gr.a <- grad_a(fit$par.list$a, a = pred_naive$age[i], k = 3)
  }
  se <- sqrt(t(gr.a) %*% v.a %*% gr.a)
  naive.df.plot$SE[i] <- se
  naive.df.plot$asymp.LCL[i] <- pmax(0, naive.df.plot$phi[i] - qnorm(0.975) * se)
  naive.df.plot$asymp.UCL[i] <- pmin(1, naive.df.plot$phi[i] + qnorm(0.975) * se)
}
colnames(naive.df.plot) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")

v.g <- fit$variance$stacked.v.est[25:28,25:28]
for (i in 1:nrow(pred_cc)) {
  if (i %% 3 == 1) {
    gr.g <- grad_g(fit$par.list$c, x1 = pred_cc$age[i], x2 = 35, k = 1)
  } else if (i %% 3 == 2) {
    gr.g <- grad_g(fit$par.list$c, x1 = pred_cc$age[i], x2 = 35, k = 2)
  } else {
    gr.g <- grad_g(fit$par.list$c, x1 = pred_cc$age[i], x2 = 35, k = 3)
  }
  se <- sqrt(t(gr.g) %*% v.g %*% gr.g)
  cc.df.plot$SE[i] <- se
  cc.df.plot$asymp.LCL[i] <- pmax(0, cc.df.plot$phi[i] - qnorm(0.975) * se)
  cc.df.plot$asymp.UCL[i] <- pmin(1, cc.df.plot$phi[i] + qnorm(0.975) * se)
}
colnames(cc.df.plot) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")

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
cc.df.plot$Method <- "CC"

subset.cured <- pred.df1[pred.df1$THR == 0, -3]
cured.df.plot <- do.call(rbind, apply(subset.cured, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.alpha <- grad(grad.cured.acat, x = fit$par$a, data = data_row, k=3)
  v.alpha <- fit$variance$stacked.v.est[13:15,13:15]
  se <- sqrt(g.alpha %*% v.alpha %*% g.alpha)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age","phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Cured")
})) %>% mutate(pred.df1[pred.df1$THR == 0,"sf36_RP"], .before = 2)
colnames(cured.df.plot) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")
#cured.df.plot <- cured.df.plot %>% mutate(sf36_RP = ifelse(is.na(sf36_RP), "(25,50,75)", sf36_RP))

subset.uncured.1 <- pred.df1[pred.df1$THR == 1 & pred.df1$age <= pred.df1$age.THR, -3]
uncured.df.plot.1 <- do.call(rbind, apply(subset.uncured.1, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.beta <- grad(grad.uncured.acat, x = params, data = data_row, k=3)[1:9]
  v.beta <- fit$variance$stacked.v.est[c(22:24,16:21),c(22:24,16:21)]
  se <- sqrt(g.beta %*% v.beta %*% g.beta)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age","phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Uncured")
})) %>% mutate(pred.df1[pred.df1$THR == 1 & pred.df1$age <= pred.df1$age.THR, "sf36_RP"], .before = 2)
colnames(uncured.df.plot.1) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")
#uncured.df.plot.1 <- uncured.df.plot.1 %>% mutate(sf36_RP = ifelse(is.na(sf36_RP), "(25,50,75)", sf36_RP))

grad.uncured.acat(par = params, data = subset.uncured.2[14,], k=3)

subset.uncured.2 <- pred.df1[pred.df1$THR == 1 & pred.df1$age > pred.df1$age.THR, -3]
uncured.df.plot.2 <- do.call(rbind, apply(subset.uncured.2, 1, function(row) {
  data_row <- as.data.frame(lapply(as.data.frame(t(row)), as.numeric))
  g.gamma <- grad(grad.uncured.acat, x = params, data = data_row, k=3)[10:19]
  v.gamma <- fit$variance$stacked.v.est[c(25:28,16:21),c(25:28,16:21)]
  se <- sqrt(g.gamma %*% v.gamma %*% g.gamma)
  LCL <- pmax(0, data_row["phi"] - qnorm(0.975) * se)
  UCL <- pmin(1, data_row["phi"] + qnorm(0.975) * se)
  data.frame(data_row[c("age","phi")], SE=se, df=Inf, LCL=LCL, UCL=UCL, Method = "2stage Uncured")
})) %>% mutate(pred.df1[pred.df1$THR == 1 & pred.df1$age > pred.df1$age.THR, "sf36_RP"], .before = 2)
colnames(uncured.df.plot.2) <- c("age", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL", "Method")
#uncured.df.plot.2 <- uncured.df.plot.2 %>% mutate(sf36_RP = ifelse(is.na(sf36_RP), "(25,50,75)", sf36_RP))

plot.df <- rbind(naive.df.plot, cc.df.plot, cured.df.plot, uncured.df.plot.1, uncured.df.plot.2)
plot.df$sf36_RP <- factor(
  plot.df$sf36_RP,
  levels = c(1, 2, 3),
  labels = c(0, "25-75", 100)
)

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
  labs(y = "Estimated probability", x = "Age", color = "Method", fill = "Method", linetype = "Method", 
       title = "Estimated PO Outcome Probabilities by Age\n When age of THR is 35",
       subtitle = "in the selected model with sex:age interaction",
       caption ="The outcome probabilities are P(sf36_RP|age.THR=35,sex=female,dysplasia=0,unilateral=0,activity.restrict=0,surgery=0)") +
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

