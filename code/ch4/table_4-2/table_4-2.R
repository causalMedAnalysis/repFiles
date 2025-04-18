# Preliminaries
chapter <- "ch4"
title <- "table_4-2"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch4/table_4-2.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R

# Outputs:     .../code/ch4/_LOGS/table_4-2_log.txt

# Description: Replicates Chapter 4, Table 4.2: Interventional Effects of 
#              College Attendance on CES-D Scores as Estimated from the NLSY 
#              Using Regression-with-Residuals.
#-------------------------------------------------------------------------------


#------------------------#
#  INSTALL DEPENDENCIES  #
#------------------------#
# This script uses the rwrmed R package, which is available to install from 
# GitHub.
# To install the package directly, you must first have installed the devtools 
# package (which is available on CRAN).

#install.packages("devtools")
# ^ Uncomment this line above to install the devtools package, if you have not 
# already done so.

#devtools::install_github("xiangzhou09/rwrmed")
# ^ Uncomment this line above to install the rwrmed package from GitHub.




#-------------#
#  LIBRARIES  #
#-------------#
library(tidyverse)
library(haven)




#-----------------------------#
#  LOAD CAUSAL MED FUNCTIONS  #
#-----------------------------#
# utilities
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")
# regression-with-residuals estimator
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R")
# ^ Note that rwrlite() is a wrapper for two functions from the rwrmed package. 
# It requires that you have installed rwrmed. (But you do not need to load the 
# rwrmed package beforehand, with the library function.)




#------------------#
#  SPECIFICATIONS  #
#------------------#
# outcome
Y <- "std_cesd_age40"

# exposure
D <- "att22"

# mediator
M <- "faminc_adj_age3539"
ln_M <- "log_faminc_adj_age3539"

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
m <- 5e4
ln_m <- log(m)




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
# RWR with D x M interaction

# L and M model formulae
predictors1_LM <- paste(c(D,C), collapse = " + ")
(formula1_L_string <- paste(L, "~", predictors1_LM))
(formula1_M_string <- paste(M, "~", predictors1_LM))
formula1_L <- as.formula(formula1_L_string)
formula1_M <- as.formula(formula1_M_string)

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
formula1_Y <- as.formula(formula1_Y_string)

# Estimate effects
out1 <- rwrlite(
  data = nlsy,
  D = D,
  C = C,
  m = m,
  Y_formula = formula1_Y,
  M_formula = formula1_M,
  L_formula_list = list(formula1_L)
)




#-------------#
#  VERSION 2  #
#-------------#
# RWR with ln(M) and D x ln(M) interaction

# L and M model formulae
(formula2_M_string <- paste(ln_M, "~", predictors1_LM))
formula2_L <- formula1_L
formula2_M <- as.formula(formula2_M_string)

# Y model formula
## main effects
predictors2_Y <- paste(c(D,ln_M,L,C), collapse = " + ")
## D x M interaction
predictors2_Y <- paste(
  predictors2_Y,
  "+",
  paste(D, ln_M, sep = ":", collapse = " + ")
)
## full formula
(formula2_Y_string <- paste(Y, "~", predictors2_Y))
formula2_Y <- as.formula(formula2_Y_string)

# Estimate effects
out2 <- rwrlite(
  data = nlsy,
  D = D,
  C = C,
  m = ln_m,
  Y_formula = formula2_Y,
  M_formula = formula2_M,
  L_formula_list = list(formula2_L)
)




#-------------#
#  VERSION 3  #
#-------------#
# RWR with ln(M) and two-way interactions
# L model:     With D x C interactions
# ln(M) model: With D x C interactions
# Y model:     With D x ln(M), D x C, ln(M) x C, ln(M) x L interactions

# L and M model formulae
## main effects
predictors3_LM <- paste(c(D,C), collapse = " + ")
## D x C interactions
predictors3_LM <- paste(
  predictors3_LM,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## full formula
(formula3_L_string <- paste(L, "~", predictors3_LM))
(formula3_M_string <- paste(ln_M, "~", predictors3_LM))
formula3_L <- as.formula(formula3_L_string)
formula3_M <- as.formula(formula3_M_string)

# Y model formula
## main effects
predictors3_Y <- paste(c(D,ln_M,L,C), collapse = " + ")
## D x ln(M) interaction
predictors3_Y <- paste(
  predictors3_Y,
  "+",
  paste(D, ln_M, sep = ":", collapse = " + ")
)
## D x C interactions
predictors3_Y <- paste(
  predictors3_Y,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## ln(M) x C interactions
predictors3_Y <- paste(
  predictors3_Y,
  "+",
  paste(ln_M, C, sep = ":", collapse = " + ")
)
## ln(M) x L interaction
predictors3_Y <- paste(
  predictors3_Y,
  "+",
  paste(ln_M, L, sep = ":", collapse = " + ")
)
## full formula
(formula3_Y_string <- paste(Y, "~", predictors3_Y))
formula3_Y <- as.formula(formula3_Y_string)

# Estimate effects
out3 <- rwrlite(
  data = nlsy,
  D = D,
  C = C,
  m = ln_m,
  Y_formula = formula3_Y,
  M_formula = formula3_M,
  L_formula_list = list(formula3_L)
)




#---------------------#
#  COLLATE ESTIMATES  #
#---------------------#
master <- data.frame(
  param = c("OE(1,0)", "IDE(1,0)", "IIE(1,0)", "CDE(1,0,50K)"),
  est_v1 = c(
    out1$OE,
    out1$IDE,
    out1$IIE,
    out1$CDE
  ),
  est_v2 = c(
    out2$OE,
    out2$IDE,
    out2$IIE,
    out2$CDE
  ),
  est_v3 = c(
    out3$OE,
    out3$IDE,
    out3$IIE,
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

