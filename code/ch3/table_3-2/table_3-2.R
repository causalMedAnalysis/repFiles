# Preliminaries
chapter <- "ch3"
title <- "table_3-2"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/table_3-2.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/linmed.R

# Outputs:     .../code/ch3/_LOGS/table_3-2_log.txt

# Description: Replicates Chapter 3, Table 3-2: Total, Direct, and Indirect 
#              Effects of College Attendance on CES-D Scores as Estimated from 
#              Linear Models Fit to the NLSY.
#-------------------------------------------------------------------------------


#-------------#
#  LIBRARIES  #
#-------------#
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




#----------------#
#  PREPARE DATA  #
#----------------#
nlsy_raw <- read_stata(
  #file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
  file = "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Data/NLSY79/nlsy79BK_ed2.dta"
)

nlsy <- nlsy_raw[complete.cases(nlsy_raw[,key_vars]),] |>
  mutate(
    std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
  )




#-------------------------#
#  ADDITIVE LINEAR MODEL  #
#-------------------------#
out1 <- linmed(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  C = C,
  m = m
)




#--------------------------------#
#  INTERACTIVE MODEL: Version A  #
#--------------------------------#
# Linear model with D x M interaction
out2 <- linmed(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  C = C,
  m = m,
  interaction_DM = TRUE
)




#--------------------------------#
#  INTERACTIVE MODEL, Version B  #
#--------------------------------#
# Linear model with D x M, C x D, C x M interactions
out3 <- linmed(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  C = C,
  m = m,
  interaction_DM = TRUE,
  interaction_DC = TRUE,
  interaction_MC = TRUE
)




#---------------------#
#  COLLATE ESTIMATES  #
#---------------------#
master <- data.frame(
  param = c("ATE(1,0)", "NDE(1,0)", "NIE(1,0)", "CDE(1,0,0)"),
  est_add = c(
    out1$ATE,
    out1$NDE,
    out1$NIE,
    out1$CDE
  ),
  est_intX = c(
    out2$ATE,
    out2$NDE,
    out2$NIE,
    out2$CDE
  ),
  est_intXX = c(
    out3$ATE,
    out3$NDE,
    out3$NIE,
    out3$CDE
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

