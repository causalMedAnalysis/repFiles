# Preliminaries
chapter <- "ch3"
title <- "table_3-7"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/table_3-7.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/JOBSII/Jobs-NoMiss-Binary.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/linmed.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/impcde.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwmed.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwcde.R

# Outputs:     .../code/ch3/_LOGS/table_3-7_log.txt

# Description: Replicates Chapter 3, Table 3-7: Total, Direct, and Indirect 
#              Effects of Job Training on Employment as Estimated from JOBSII.
#-------------------------------------------------------------------------------


#-------------#
#  LIBRARIES  #
#-------------#
library(mediation)
library(tidyverse)
library(haven)




#-----------------------------#
#  LOAD CAUSAL MED FUNCTIONS  #
#-----------------------------#
# utilities
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/utils_bare.R")
# product-of-coefficients estimator, based on linear models
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/linmed.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/linmed.R")
# regression imputation CDE estimator
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/impcde.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/impcde.R")
# IPW estimator
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwmed.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/ipwmed.R")
# IPW CDE estimator
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwcde.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/ipwcde.R")




#------------------#
#  SPECIFICATIONS  #
#------------------#
# outcome
Y <- "work1"

# exposure
D <- "treat"

# mediator
M <- "job_seek"

# baseline confounder(s)
C <- c(
  "econ_hard",
  "sex",
  "age",
  "nonwhite",
  "educ",
  "income"
)

# mediator value for CDE
m <- 4

# number of simulations for simulation estimator
n_sims <- 1000

# number of bootstrap replications
n_reps <- 2000




#----------------#
#  PREPARE DATA  #
#----------------#
jobs_raw <- read_stata(
  #file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/JOBSII/Jobs-NoMiss-Binary.dta"
  file = "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Data/JOBSII/Jobs-NoMiss-Binary.dta"
)

jobs <- jobs_raw |>
  zap_labels()




#--------------------------#
#  LINEAR MODEL ESTIMATOR  #
#--------------------------#
# Linear model with D x M interaction

out_lin <- linmed(
  data = jobs,
  D = D,
  M = M,
  Y = Y,
  C = C,
  m = m,
  interaction_DM = TRUE,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
  # ^ Note that parallelizing the bootstrap is optional, but requires that you 
  # have installed the following R packages: doParallel, doRNG, foreach.
  # (You do not need to load those packages beforehand, with the library 
  # function.)
  # If you choose not to parallelize the bootstrap (by setting the boot_parallel 
  # argument to FALSE), the results may differ slightly, due to simulation 
  # variance (even if you specify the same seed).
)




#-------------------------------------------------#
#  SIMULATION & REGRESSION IMPUTATION ESTIMATORS  #
#-------------------------------------------------#
# M model: Additive linear model
# Y model: Logit model with D x M interaction

# Mediator model formula
predictors_M <- paste(c(D,C), collapse = " + ")
formula_M_string <- paste(M, "~", predictors_M)
formula_M_string

# Outcome model formula
## main effects
predictors_Y <- paste(c(D,M,C), collapse = " + ")
## D x M interaction
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(D, M, sep = ":", collapse = " + ")
)
## full formula
formula_Y_string <- paste(Y, "~", predictors_Y)
formula_Y_string

# Fit mediator and outcome models
mod_M <- lm(
  as.formula(formula_M_string),
  data = jobs
)
mod_Y <- glm(
  #as.formula(formula_Y_string),
  work1 ~ treat*job_seek + econ_hard + sex + age + nonwhite + educ + income,
  # ^ scoping issues with the update function (used in the impcde bootstrap) 
  # require us to directly specify the formula, rather than reference the 
  # formula_Y_string object
  family = binomial(link = "logit"),
  data = jobs
)

# Estimate ATE(1,0), NDE(1,0), and NIE(1,0) by simulation estimator
set.seed(3308004)
out_sim <- mediate(
  model.m = mod_M,
  model.y = mod_Y,
  sims = n_sims,
  boot = TRUE,
  treat = D,
  mediator = M,
  parallel = "multicore"
)

# Estimate CDE(1,0,4) by regression imputation estimator
out_imp_cde <- impcde(
  data = jobs,
  model_y = mod_Y,
  D = D,
  M = M,
  m = m,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
  # ^ See note above about parallelizing the bootstrap.
)




#-----------------#
#  IPW ESTIMATOR  #
#-----------------#
# Additive logit models

# D model 1 formula: f(D|C)
predictors1_D <- paste(C, collapse = " + ")
formula1_D_string <- paste(D, "~", predictors1_D)
formula1_D_string

# D model 2 formula: s(D|C,M)
predictors2_D <- paste(c(M,C), collapse = " + ")
formula2_D_string <- paste(D, "~", predictors2_D)
formula2_D_string

# M model formula: g(M|C,D)
formula_M_string # defined above

# Estimate ATE(1,0), NDE(1,0), NIE(1,0)
out_ipw <- ipwmed(
  data = jobs,
  D = D,
  M = M,
  Y = Y,
  formula1_string = formula1_D_string,
  formula2_string = formula2_D_string,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
  # ^ See note above about parallelizing the bootstrap.
)

# Estimate CDE(1,0,0)
# out_ipw_cde <- ipwcde(
#   data = jobs,
#   D = D,
#   M = M,
#   Y = Y,
#   formula_D_string = formula1_D_string,
#   formula_M_string = formula_M_string,
#   boot = TRUE,
#   boot_reps = n_reps,
#   boot_seed = 3308004,
#   boot_parallel = TRUE
#   # ^ See note above about parallelizing the bootstrap.
# )




#-------------------#
#  COLLATE RESULTS  #
#-------------------#
master <- data.frame(
  param = c("ATE(1,0)", "NDE(1,0)", "NIE(1,0)", "CDE(1,0,4)"),
  
  # linear models
  lin_est = c(
    out_lin$ATE,
    out_lin$NDE,
    out_lin$NIE,
    out_lin$CDE
  ),
  lin_ci_low = c(
    out_lin$ci_ATE[1],
    out_lin$ci_NDE[1],
    out_lin$ci_NIE[1],
    out_lin$ci_CDE[1]
  ),
  lin_ci_high = c(
    out_lin$ci_ATE[2],
    out_lin$ci_NDE[2],
    out_lin$ci_NIE[2],
    out_lin$ci_CDE[2]
  ),
  
  # simulation and regression imputation
  sim_est = c(
    out_sim$tau.coef,
    out_sim$z0,
    out_sim$d1,
    out_imp_cde$CDE
  ),
  sim_ci_low = c(
    out_sim$tau.ci[1],
    out_sim$z0.ci[1],
    out_sim$d1.ci[1],
    out_imp_cde$ci_CDE[1]
  ),
  sim_ci_high = c(
    out_sim$tau.ci[2],
    out_sim$z0.ci[2],
    out_sim$d1.ci[2],
    out_imp_cde$ci_CDE[2]
  ),
  
  # IPW
  ipw_est = c(
    out_ipw$ATE,
    out_ipw$NDE,
    out_ipw$NIE,
    NA_real_
  ),
  ipw_ci_low = c(
    out_ipw$ci_ATE[1],
    out_ipw$ci_NDE[1],
    out_ipw$ci_NIE[1],
    NA_real_
  ),
  ipw_ci_high = c(
    out_ipw$ci_ATE[2],
    out_ipw$ci_NDE[2],
    out_ipw$ci_NIE[2],
    NA_real_
  )
)

master |>
  mutate(
    across(
      .cols = !param,
      .fns = \(x) round(x, 3)
    )
  )


# Close log
sink()

