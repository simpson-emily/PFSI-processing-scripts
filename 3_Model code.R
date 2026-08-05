
###########################
# Supplementary to: Simpson et al. 2026: Fire - Moisture Interactions Shape Recovery Trajectories Across Temperate Vegetation Communities, Global Change Biology
# Code takes dataset [produced from JavaScripts 1 & 2] and generates model for PFSI in relation to vegetation community, fire severity and soil moisture
# your_path_here.csv >> insert provided dataset ["Final_Dataset"]
###########################

library(cmdstanr)
library(dplyr)
library(brms)
library(corrplot)
library(ggplot2)

#To install cmdstanr
# install.packages(
#  "cmdstanr",
#  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
# )
# library(cmdstanr)
# cmdstanr::install_cmdstan(cores = 2)
# cmdstanr::check_cmdstan_toolchain() #check it was installed properly

df_final <- read.csv("your_path_here.csv")

# **** Determine if soil moisture should be lagged ****

selected_cols <- df_final[, c("SM", "SMLag1", "SMLag1Seas", "PFSI")]
cor_matrix <- cor(selected_cols, use = "complete.obs", method = "pearson")
corrplot(cor_matrix, method = "number", type = "upper", tl.col = "black", tl.srt = 45)

df_final$SM <- scale(df_final$SM, center = TRUE, scale = TRUE)
df_final$SM_lag1 <- scale(df_final$SMLag1, center = TRUE, scale = TRUE)
df_final$SM_lag1seas <- scale(df_final$SMLag1Seas, center = TRUE, scale = TRUE)

# Model with no lag
sm_model <- lm(
  PFSI ~ vegTypeN + SM + fireSevN,
  data = df_final
)
summary(sm_model)

# Model with yearly lag
year_model <- lm(
  PFSI ~ vegTypeN + SMLag1 + fireSevN,
  data = df_final
)
summary(year_model)

# Model with seasonal lag
seas_model <- lm(
  PFSI ~ vegTypeN + SMLag1Seas + fireSevN,
  data = df_final
)
summary(seas_model)

# To get values for comparison table

model_list <- list(sm_model = sm_model, year_model = year_model, seas_model = seas_model)

model_comparison <- data.frame(
  Model = names(model_list),
  k = sapply(model_list, function(m) length(coef(m))),           # number of parameters
  LogLik = sapply(model_list, logLik),                           # log-likelihood
  AIC = sapply(model_list, AIC)                                  # AIC values
)

# Delta AIC
min_aic <- min(model_comparison$AIC)
model_comparison$Delta_AIC <- model_comparison$AIC - min_aic

# AIC weights
model_comparison$AIC_weight <- exp(-0.5 * model_comparison$Delta_AIC) / sum(exp(-0.5 * model_comparison$Delta_AIC))

model_comparison <- model_comparison %>%
  mutate(
    LogLik = round(LogLik, 2),
    AIC = round(AIC, 2),
    Delta_AIC = round(Delta_AIC, 2)
  )

print(model_comparison)

# **** Best model was soil moisture with yearly lag ****
# **** Prepare data for modelling ****

# Set dataset and reference levels 

df_model <- df_final %>% 
  mutate(
    Year_c      = year - 2019.75,                  
    fireSevN   = factor(fireSevN), 
    vegTypeN   = factor(vegTypeN),
    SMLag1 = SMLag1             
  )

df_model$fireSevN <- relevel(df_model$fireSevN, ref = "Unburnt")
df_model$vegTypeN <- relevel(df_model$vegTypeN, ref = "Dry Sclerophyll Forest")

# Test model priors

bf_recov <- bf(
  PFSI ~
    s(Year_c, by = fireSevN, k = 6) +
    fireSevN * vegTypeN * SMLag1 +
    (1 | sample_id),
  family = student(),
  decomp = "QR",
  center = TRUE
)

priors <- c(
  set_prior("student_t(3,  0, 10)", class = "b"),     
  set_prior("normal(0, 80)", class = "Intercept"),      
  set_prior("normal(115, 30)", class = "sigma", lb = 0), 
  set_prior("gamma(2, 0.7)", class = "nu"),              
  set_prior("student_t(3,  0,  6)", class = "sd")       
)

prior_model <- brm(
  bf_recov,
  data    = df_model,
  prior   = priors,
  sample_prior = "only",
  iter    = 4000,
  warmup  = 1000,
  chains = 2,
  seed = 123,
  control = list(adapt_delta = 0.95)
)

pp_check(prior_model, type = "dens_overlay") +
  xlim(c(-2000,2000))

# **** Run model ****

# Model is large; can be slow to run. Running chains parallel can speed up convergence.

options(
  mc.cores     = parallel::detectCores(),
  brms.backend = "cmdstanr"
)

final_model <- brm(
  bf_recov,
  data    = df_model,
  prior   = priors,
  iter    = 6000,
  warmup  = 2000,
  chains = 4,
  cores = 4,
  threads = threading(2),
  seed = 123,
  control = list(adapt_delta = 0.995, max_treedepth = 25) 
)

summary(final_model)

#pp_check(final_model, type = "dens_overlay") +
#xlim(c(-2000,2000))
