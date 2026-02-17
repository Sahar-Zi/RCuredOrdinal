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

sf36_levels <- levels(X$sf36_RP)
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
             
# ----------------------- Create data grid for plots -------------------------

grid.data <- expand.grid(age = seq(from = 18, to = 76, by = 0.2),
                         THR = c(0, 1),
                         sex = c(0, 1),
                         sf36_RP_V = sf36_levels) %>% # Non-constant variables
  arrange(sex, THR, age) %>%
  mutate(
    age.THR = AGE_THR_FIXED, 
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
  outcome.model = OUTCOME_MODEL)

grid.data$phi <- NA
grid.data$phi[grid.data$THR==1] <- v.probs.est$D1[cbind(which(grid.data$THR==1), grid.data$sf36_RP[grid.data$THR==1])]
grid.data$phi[grid.data$THR==0] <- v.probs.est$D0[cbind(which(grid.data$THR==0), grid.data$sf36_RP[grid.data$THR==0])]

# ----------------------- Estimate naive outcome probabilities -------------------------

grid.data.naive <- grid.data %>% filter(THR == 0) 
if(OUTCOME_MODEL == "PO"){
  mod.naive <- polr(sf36_RP ~ age, data = X %>% filter(THR == 0))
  naive.probs <- predict(mod.naive, newdata = grid.data.naive, type = "probs")
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
  }
  
  return(grad)
}

naive.df.plot <- grid.data.naive[,c("age","sex","sf36_RP","phi")]
naive.df.plot$df <- naive.df.plot$SE <- Inf
v.a <- vcov(mod.naive)

for (i in 1:nrow(grid.data.naive)) {
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
  }
  se <- sqrt(t(gr.g) %*% v.g %*% gr.g)
  cc.df.plot$SE[i] <- se
  cc.df.plot$asymp.LCL[i] <- pmax(0, cc.df.plot$phi[i] - qnorm(0.975) * se)
  cc.df.plot$asymp.UCL[i] <- pmin(1, cc.df.plot$phi[i] + qnorm(0.975) * se)
}
colnames(cc.df.plot) <- c("age", "sex", "sf36_RP", "prob", "SE", "df", "asymp.LCL", "asymp.UCL")
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


subset.uncured.2 <- grid.data[grid.data$THR == 1 & grid.data$age >= grid.data$age.THR, -4]
variance.uncured.gamma <- fit$variance$stacked.v.est[c(25:28,16:21),c(25:28,16:21)] # The sub-matrix of the gamma parameters
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
  facet_grid(
    sex ~ sf36_RP,
    labeller = \(x) label_both(x, sep = " = ")
  ) +
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
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.6, "cm"),
    panel.spacing = unit(1, "cm"),
    axis.line = element_line(linewidth = 0.4)
  )

# # ----------------------- Save the plot -------------------------
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
