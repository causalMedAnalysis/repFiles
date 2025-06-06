# Preliminaries
chapter <- "ch4"
title <- "table_4-3"
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch4/table_4-3.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R

# Outputs:     .../code/ch4/_LOGS/table_4-3_log.txt

# Description: Replicates Chapter 4, Table 4.3: Interventional Effects of 
#              College Attendance on CES-D Scores as Estimated from the NLSY 
#              Using the Simulation Approach.
#-------------------------------------------------------------------------------

#------------------------#
#  INSTALL DEPENDENCIES  #
#------------------------#
# The following packages are used to parallelize the bootstrap.
dependencies <- c("doParallel", "doRNG", "foreach")

#install.packages(dependencies)
# ^ Uncomment this line above to install these packages.

# And note that, once you have installed these packages, there is no need for 
# you to load these packages with the library function to run the code in this 
# script.

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

#------------------#
#  SPECIFICATIONS  #
#------------------#
# outcome
Y <- "std_cesd_age40"

# exposure
D <- "att22"

# mediator
M <- "log_faminc_adj_age3539"

# exposure-induced confounder
L <- "ever_unemp_age3539"

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
  L,
  C
)

# mediator value for CDE
m <- log(5e4)

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
# L model: Additive logit model
# M model: Additive linear model
# Y model: Linear model with D x M interaction

# L and M model formulae
predictors1_LM <- paste(c(D,C), collapse = " + ")
(formula1_L_string <- paste(L, "~", predictors1_LM))
(formula1_M_string <- paste(M, "~", predictors1_LM))

# Y model formula
## main effects
predictors1_Y <- paste(c(D,M,L,C), collapse = " + ")
## D x M interaction
predictors1_Y <- paste(
  predictors1_Y,
  "+",
  paste(D, M, sep = ":", collapse = " + ")
)
## full formula
(formula1_Y_string <- paste(Y, "~", predictors1_Y))

# Define model specifications
out1_specs <- list(
  ## L model
  list(
    func = "glm",
    formula = as.formula(formula1_L_string),
    args = list(family = "binomial")
  ),
  ## M model
  list(
    func = "lm",
    formula = as.formula(formula1_M_string)
  ),
  ## Y model
  list(
    func = "lm",
    formula = as.formula(formula1_Y_string)
  )
)

# Estimate interventional effects
out1 <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = M,
  model_spec = out1_specs,
  seed = 3308004
)

# Estimate CDE(1,0,ln(50K))
out1_cde <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = paste0(M,"=m"),
  model_spec = out1_specs,
  seed = 3308004
)

#-------------#
#  VERSION 2  #
#-------------#
# L model: Logit model with D x C interactions
# M model: Linear model with D x C interactions
# Y model: Linear model with D x M, D x C, M x C, M x L interactions

# L and M model formulae
## main effects
predictors2_LM <- paste(c(D,C), collapse = " + ")
## D x C interactions
predictors2_LM <- paste(
  predictors2_LM,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## full formula
(formula2_L_string <- paste(L, "~", predictors2_LM))
(formula2_M_string <- paste(M, "~", predictors2_LM))

# Y model formula
## main effects
predictors2_Y <- paste(c(D,M,L,C), collapse = " + ")
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
## M x L interaction
predictors2_Y <- paste(
  predictors2_Y,
  "+",
  paste(M, L, sep = ":", collapse = " + ")
)
## full formula
(formula2_Y_string <- paste(Y, "~", predictors2_Y))

# Define model specifications
out2_specs <- list(
  ## L model
  list(
    func = "glm",
    formula = as.formula(formula2_L_string),
    args = list(family = "binomial")
  ),
  ## M model
  list(
    func = "lm",
    formula = as.formula(formula2_M_string)
  ),
  ## Y model
  list(
    func = "lm",
    formula = as.formula(formula2_Y_string)
  )
)

# Estimate interventional effects
out2 <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = M,
  model_spec = out2_specs,
  seed = 3308004
)

# Estimate CDE(1,0,ln(50K))
out2_cde <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = paste0(M,"=m"),
  model_spec = out2_specs,
  seed = 3308004
)

#---------------------#
#  COLLATE ESTIMATES  #
#---------------------#
master <- data.frame(
  param = c("OE(1,0)", "IDE(1,0)", "IIE(1,0)", "CDE(1,0,ln(50K))"),
  est_v1 = c(
    out1[[3]],
    out1[[1]],
    out1[[2]],
    out1_cde[[1]]
  ),
  est_v2 = c(
    out2[[3]],
    out2[[1]],
    out2[[2]],
    out2_cde[[1]]
  )
)

# Open log
sink(log_path, split = TRUE)

master |>
  mutate(
    across(
      .cols = starts_with("est_"),
      .fns = \(x) round(x, 3)
    )
  )

# Close log
sink()

