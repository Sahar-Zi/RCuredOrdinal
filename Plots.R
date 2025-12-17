# Load required libraries
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(gridExtra)
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

# Need to load the RData file

# ---- Simulation plots
Bias.CI <- function(est, se, true, reps){
  u <- est + qnorm(0.975) * se
  l <- est - qnorm(0.975) * se
  r <- matrix(true, nrow = reps, ncol = length(true), byrow = T)
  return(data.frame(bias = colMeans(est, na.rm=T)-true, l = colMeans(l,na.rm=T)-true, u = colMeans(u,na.rm=T)-true))
}

# 1. proposed_bias assumed given:
proposed.bias <- cbind(Bias.CI(est.values, est.se.values, unlist(true.params), reps = 131),Method = "Proposed")[-1:-7,]

# # 2. Prepare full covariate list (from proposed)
all_covariates <- rownames(proposed.bias)

# # 3. Create naive_bias with NA for missing covariates
naive.bias <- proposed.bias
naive.bias[,] <- NA
naive.bias$Method <- "Naive"
naive.subset <- Bias.CI(naive.est.values, naive.est.se.values, unlist(true.params$a), reps = 131)
naive.bias[intersect(rownames(naive.bias), paste0("a.", rownames(naive.subset))),-4] <- naive.subset

# # 4. Create cc_bias similarly
CC.bias <- proposed.bias
CC.bias[,] <- NA
CC.bias$Method <- "Complete-Case"
CC.subset <- Bias.CI(cc.est.values, cc.est.se.values, unlist(true.params$c), reps = 131)
CC.bias[intersect(rownames(CC.bias), paste0("c.", rownames(CC.subset))),-4] <- CC.subset

# # 5. Combine all
proposed.bias$coefficient <- rownames(proposed.bias)
rownames(proposed.bias) <- NULL
naive.bias$coefficient <- rownames(naive.bias)
rownames(naive.bias) <- NULL
CC.bias$coefficient <- rownames(CC.bias)
rownames(CC.bias) <- NULL
bias_all <- bind_rows(proposed.bias, naive.bias, CC.bias)
bias_all$coefficient <- factor(bias_all$coefficient, levels = unique(bias_all$coefficient))
#bias_all$Method <- factor(bias_all$Method, levels = unique(bias_all$Method))

# # 6. Plot
ggplot(bias_all %>% arrange(Method), aes(x = coefficient, y = bias, color = Method, shape = Method)) +
  annotate("rect", xmin = 0.5, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "#c6dbef", alpha = 0.06) +
  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = -Inf, ymax = Inf, fill = "#fdd0a2", alpha = 0.06) +
  annotate("rect", xmin = 7.5, xmax = 14.5, ymin = -Inf, ymax = Inf, fill = "#d9f0d3", alpha = 0.06) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "black") +
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.35, size = 0.8, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 4, alpha = 1, na.rm = TRUE) +
  geom_vline(xintercept = 5.5, linetype = "longdash", color = "salmon") +
  geom_vline(xintercept = 7.5, linetype = "longdash", color = "salmon") +
  ylim(-2, 2) +
  scale_x_discrete(position = "top", labels = c(
    "a.V.alpha1"     = expression(alpha[1]),
    "a.V.alpha2"     = expression(alpha[2]),
    "a.V.alphaZ1"    = expression(alpha["Z"[1]]),
    "a.V.alphaZ2"    = expression(alpha["Z"[2]]),
    "a.V.alphaR"     = expression(alpha[R]),
    "b.V.beta1"      = expression(beta[1]),
    "b.V.beta2"      = expression(beta[2]),
    "c.V.gamma1"     = expression(gamma[1]),
    "c.V.gamma2"     = expression(gamma[2]),
    "c.V.gammaZ1"    = expression(gamma["Z"[1]]),
    "c.V.gammaZ2"    = expression(gamma["Z"[2]]),
    "c.V.gammaR"     = expression(gamma[R]),
    "c.V.gammaTau"   = expression(gamma["T"]),
    "c.V.gammaR:Tau" = expression(gamma[R * "T"])
  )) +
  labs(x = NULL, y = "Bias", 
       title = "Bias and 95% Confidence Intervals for Estimated Coefficients",
       subtitle = expression("Results based on 400 replications with sample size n = 2000 of " * 
                               bold("Proportional Odds") * " outcome model"),
       caption = "Standard errors and confidence intervals were computed using the sandwich variance estimator specifically developed for the two-stage proposed method.") +
  theme_minimal(base_size = 15) +
  theme(plot.caption = element_text(hjust = 0),
        axis.text.x = element_text(size = 14),
        legend.background = element_rect(
          fill = "white",       # or "gray95" or any color you like
          color = "black",      # color of the border
          linewidth = 0.5,      # thickness of the border
          linetype = "solid"    # style of the border (solid, dashed, etc.)
        ),
        legend.position = "top"
  ) +
  scale_color_manual(values = c("Proposed" = "#0072B2", "Naive" = "#D55E00", "Complete-Case" = "#CC79A7"),
                     labels = c("Proposed" = "2stage Proposed", "Naive" = "Naive", "Complete-Case" = "Complete Case")) +
  scale_shape_manual(values = c("Proposed" = 16, "Naive" = 17, "Complete-Case" = 15),
                     labels = c("Proposed" = "2stage Proposed", "Naive" = "Naive", "Complete-Case" = "Complete Case"))




# ---- Simulation plots 2  
#plot_CP_with_CIs <- function(sim, title = "") {
CI.CP <- function(est, se, true, reps) {
u <- est + qnorm(0.975) * se
l <- est - qnorm(0.975) * se
r <- matrix(rep(true, each = reps), nrow = reps)

cp <- colMeans(r >= l & r <= u)
se_cp <- sqrt(cp * (1 - cp) / reps)

cp_l <- pmax(0, cp - qnorm(0.975) * se_cp)
cp_u <- pmin(1, cp + qnorm(0.975) * se_cp)

return(data.frame(cp = cp, l = cp_l, u = cp_u))
}

# 1. proposed_bias assumed given:
proposed.CP <- cbind(CI.CP(est.values, est.se.values, unlist(true.params), reps = 400),Method = "Proposed")[-1:-7,]

# # 2. Prepare full covariate list (from proposed)
all_covariates <- rownames(proposed.CP)

# # 3. Create naive_bias with NA for missing covariates
naive.CP <- proposed.CP
naive.CP[,] <- NA
naive.CP$Method <- "Naive"
naive.subset <- CI.CP(naive.est.values, naive.est.se.values, unlist(true.params$a), reps = 400)
naive.CP[intersect(rownames(naive.CP), paste0("a.", rownames(naive.subset))),-4] <- naive.subset

# # 4. Create cc_bias similarly
CC.CP <- proposed.CP
CC.CP[,] <- NA
CC.CP$Method <- "Complete-Case"
CC.subset <- CI.CP(cc.est.values, cc.est.se.values, unlist(true.params$c), reps = 400)
CC.CP[intersect(rownames(CC.CP), paste0("c.", rownames(CC.subset))),-4] <- CC.subset

# # 5. Combine all
proposed.CP$coefficient <- rownames(proposed.CP)
rownames(proposed.CP) <- NULL
naive.CP$coefficient <- rownames(naive.CP)
rownames(naive.CP) <- NULL
CC.CP$coefficient <- rownames(CC.CP)
rownames(CC.CP) <- NULL
cp_all <- bind_rows(proposed.CP, naive.CP, CC.CP)
cp_all$coefficient <- factor(cp_all$coefficient, levels = unique(cp_all$coefficient))

custom_transform <- trans_new(
  name = "triple_scale_2",
  transform = function(x) {
    ifelse(
      x <= 0.2, x * 2,                  # stretch 0-0.2
      ifelse(
        x <= 0.8, 0.4 + (x - 0.2) * 0.5,  # compress 0.2-0.8
        0.4 + 0.3 + (x - 0.8) * 2         # stretch 0.8-1
      )
    )
  },
  inverse = function(x) {
    ifelse(
      x <= 0.4, x / 2,                       # inverse for 0-0.2
      ifelse(
        x <= 0.7, 0.2 + (x - 0.4)/0.5,      # inverse for 0.2-0.8
        0.8 + (x - 0.7)/2                    # inverse for 0.8-1
      )
    )
  }
)

# # 6. Plot
ggplot(cp_all %>% arrange(Method), aes(x = coefficient, y = cp, color = Method, shape = Method)) +
annotate("rect", xmin = 0.5, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "#c6dbef", alpha = 0.06) +
annotate("rect", xmin = 5.5, xmax = 7.5, ymin = -Inf, ymax = Inf, fill = "#fdd0a2", alpha = 0.06) +
annotate("rect", xmin = 7.5, xmax = 14.5, ymin = -Inf, ymax = Inf, fill = "#d9f0d3", alpha = 0.06) +
geom_hline(yintercept = 0.95, linetype = "longdash", color = "black") +
geom_point(size = 4, alpha = 1, na.rm = TRUE) +
geom_errorbar(aes(ymin = l, ymax = u), width = 0.35, size = 0.8, alpha = 0.8, na.rm = TRUE) +
geom_vline(xintercept = 5.5, linetype = "longdash", color = "salmon") +
geom_vline(xintercept = 7.5, linetype = "longdash", color = "salmon") +
scale_x_discrete(position = "top", labels = c(
  "a.V.alpha1"     = expression(alpha[1]),
  "a.V.alpha2"     = expression(alpha[2]),
  "a.V.alphaZ1"    = expression(alpha["Z"[1]]),
  "a.V.alphaZ2"    = expression(alpha["Z"[2]]),
  "a.V.alphaR"     = expression(alpha[R]),
  "b.V.beta1"      = expression(beta[1]),
  "b.V.beta2"      = expression(beta[2]),
  "c.V.gamma1"     = expression(gamma[1]),
  "c.V.gamma2"     = expression(gamma[2]),
  "c.V.gammaZ1"    = expression(gamma["Z"[1]]),
  "c.V.gammaZ2"    = expression(gamma["Z"[2]]),
  "c.V.gammaR"     = expression(gamma[R]),
  "c.V.gammaTau"   = expression(gamma["T"]),
  "c.V.gammaR:Tau" = expression(gamma[R * "T"])
)) +
scale_y_continuous(
  transform = custom_transform,
  limits = c(0, 1),
  breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1),
  minor_breaks = c(0.15,0.85)) +
labs(x = NULL, y = "Coverage Probability \n95% CI", 
     title = "Coverage Probability of 95% Confidence Intervals for Estimated Coefficients",
     subtitle = expression("Results based on 400 replications with sample size n = 2000 of " * 
                             bold("Propotional Odds") * " outcome models"),
     caption = "Coverage probabilities (CP) are based on 400 simulation replications.") +
theme_minimal(base_size = 15) +
theme(plot.caption = element_text(hjust = 0),
      axis.text.x = element_text(size = 14),
      legend.background = element_rect(
        fill = "white",       # or "gray95" or any color you like
        color = "black",      # color of the border
        linewidth = 0.5,      # thickness of the border
        linetype = "solid"    # style of the border (solid, dashed, etc.)
      ),
      legend.position = "top"
) +
scale_color_manual(values = c("Proposed" = "#0072B2", "Naive" = "#D55E00", "Complete-Case" = "#CC79A7"),
                   labels = c("Proposed" = "2stage Proposed", "Naive" = "Naive", "Complete-Case" = "Complete Case")) +
scale_shape_manual(values = c("Proposed" = 16, "Naive" = 17, "Complete-Case" = 15),
                   labels = c("Proposed" = "2stage Proposed", "Naive" = "Naive", "Complete-Case" = "Complete Case"))
