# Preliminaries
chapter <- "ch3"
title <- "table_3-3"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/table_3-3.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/impcde.R

# Outputs:     .../code/ch3/_LOGS/table_3-3_log.txt

# Description: Replicates Chapter 3, Table 3-3: Total, Direct, and Indirect 
#              Effects of College Attendance on CES-D Scores as Estimated from 
#              the NLSY using the Simulation and Imputation Approach.
#-------------------------------------------------------------------------------


#-------------#
#  LIBRARIES  #
#-------------#
library(tidyverse)
library(haven)




#-----------------------------#
#  LOAD CAUSAL MED FUNCTIONS  #
#-----------------------------#
# simulation estimator
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R")
# regression imputation CDE estimator
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/impcde.R")




#------------------#
#  SPECIFICATIONS  #
#------------------#
# outcome
Y <- "std_cesd_age40"

# exposure
D <- "att22"

# mediator
M <- "ever_unemp_age3539"

# baseline confounder(s)
C <- c(
  "female",
  "black",
  "hispan",
  "paredu",
  "parprof",
  "parinc_prank",
  "famsize",
  "afqt3"
)

# key variables
key_vars <- c(
  "cesd_age40", # unstandardized version of Y
  D,
  M,
  C
)

# mediator value for CDE
m <- 0

# number of simulations
n_sims <- 2000




#----------------#
#  PREPARE DATA  #
#----------------#
nlsy_raw <- read_stata(
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
)

nlsy <- nlsy_raw[complete.cases(nlsy_raw[,key_vars]),] |>
  mutate(
    std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
  )




#-------------#
#  VERSION 1  #
#-------------#
# M model: Additive logit model
# Y model: Additive linear model

# Mediator model formula
predictors1_M <- paste(c(D,C), collapse = " + ")
formula1_M_string <- paste(M, "~", predictors1_M)
formula1_M_string

# Outcome model formula
predictors1_Y <- paste(c(D,M,C), collapse = " + ")
formula1_Y_string <- paste(Y, "~", predictors1_Y)
formula1_Y_string

# Define model specifications
out1_specs <- list(
  list(
    func = "glm",
    formula = as.formula(formula1_M_string),
    args = list(family = "binomial")
  ),
  list(
    func = "lm",
    formula = as.formula(formula1_Y_string)
  )
)

# Estimate ATE(1,0), NDE(1,0), and NIE(1,0) by simulation estimator
out1 <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = NULL,
  model_spec = out1_specs,
  seed = 3308004
)

# Estimate CDE(1,0,0) by regression imputation estimator
mod1_Y <- lm(
  as.formula(formula1_Y_string),
  data = nlsy
)
out1_cde <- impcde(
  data = nlsy,
  model_y = mod1_Y,
  D = D,
  M = M,
  m = m
)




#-------------#
#  VERSION 2  #
#-------------#
# M model: Logit model with C x D interactions
# Y model: Linear model with D x M, C x D, C x M interactions

# Mediator model formula
## main effects
predictors2_M <- paste(c(D,C), collapse = " + ")
## D x C interactions
predictors2_M <- paste(
  predictors2_M,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## full formula
formula2_M_string <- paste(M, "~", predictors2_M)
formula2_M_string

# Outcome model formula
## main effects
predictors2_Y <- paste(c(D,M,C), collapse = " + ")
## D x M interaction
predictors2_Y <- paste(
  predictors2_Y,
  "+",
  paste(D, M, sep = ":", collapse = " + ")
)
## D x C interactions
predictors2_Y <- paste(
  predictors2_Y,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## M x C interactions
predictors2_Y <- paste(
  predictors2_Y,
  "+",
  paste(M, C, sep = ":", collapse = " + ")
)
## full formula
formula2_Y_string <- paste(Y, "~", predictors2_Y)
formula2_Y_string

# Define model specifications
out2_specs <- list(
  list(
    func = "glm",
    formula = as.formula(formula2_M_string),
    args = list(family = "binomial")
  ),
  list(
    func = "lm",
    formula = as.formula(formula2_Y_string)
  )
)

# Estimate ATE(1,0), NDE(1,0), and NIE(1,0)
out2 <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = NULL,
  model_spec = out2_specs,
  seed = 3308004
)

# Estimate CDE(1,0,0) by regression imputation estimator
mod2_Y <- lm(
  as.formula(formula2_Y_string),
  data = nlsy
)
out2_cde <- impcde(
  data = nlsy,
  model_y = mod2_Y,
  D = D,
  M = M,
  m = m
)




#---------------------#
#  COLLATE ESTIMATES  #
#---------------------#
master <- data.frame(
  param = c("ATE(1,0)", "NDE(1,0)", "NIE(1,0)", "CDE(1,0,0)"),
  est_v1 = c(
    out1[[1]],
    out1[[2]],
    out1[[3]],
    out1_cde
  ),
  est_v2 = c(
    out2[[1]],
    out2[[2]],
    out2[[3]],
    out2_cde
  )
)

master |>
  mutate(
    across(
      .cols = starts_with("est_"),
      .fns = \(x) round(x, 3)
    )
  )


# Close log
sink()

