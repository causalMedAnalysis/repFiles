# Preliminaries
chapter <- "ch3"
title <- "table_3-4"

#dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"
dir_root <- "C:/Users/Geoffrey Wodtke/Desktop/repFiles-Dev"

dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")

# Ensure all necessary directories exist under your root folder
# if not, the function will create folders for you

create_dir_if_missing <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message("Created directory: ", dir)
  } else {
    message("Directory already exists: ", dir)
  }
}

create_dir_if_missing(dir_root)
create_dir_if_missing(dir_log)

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/table_3-4.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch3/_LOGS/table_3-4_log.txt

# Description: Replicates Chapter 3, Table 3-4: Total, Direct, and Indirect 
#              Effects of College Attendance on CES-D Scores as Estimated from 
#              the NLSY using Inverse Probability Weighting.
#-------------------------------------------------------------------------------

#---------------------------------------------------------#
#  INSTALL DEPENDENCIES AND LOAD CAUSAL MED FUNCTIONS     #
#---------------------------------------------------------#
packages <-
  c(
    "tidyverse",
    "haven",
    "devtools"
  )

install_and_load <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing missing package: ", pkg)
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

install_and_load(packages)

install_github("causalMedAnalysis/cmedR")

library(cmedR)

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
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
)

nlsy <- nlsy_raw[complete.cases(nlsy_raw[,key_vars]),] |>
  mutate(
    std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
  )

#--------------------#
#  ESTIMATE EFFECTS  #
#--------------------#
# Additive logit models

# D model 1 formula: f(D|C)
predictors1_D <- paste(C, collapse = " + ")
formula1_D <- as.formula(paste(D, "~", predictors1_D))

# D model 2 formula: s(D|C,M)
predictors2_D <- paste(c(M,C), collapse = " + ")
formula2_D <- as.formula(paste(D, "~", predictors2_D))

# M model formula: g(M|C,D)
predictors_M <- paste(c(D,C), collapse = " + ")
formula_M <- as.formula(paste(M, "~", predictors_M))

# Estimate ATE(1,0), NDE(1,0), NIE(1,0)
out1 <- ipwmed(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  D_C_model = formula1_D,
  D_CM_model = formula2_D
)

# Estimate CDE(1,0,0)
out1_cde <- ipwcde(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  m = m,
  D_C_model = formula1_D,
  M_CD_model = formula_M
)

#---------------------#
#  COLLATE ESTIMATES  #
#---------------------#
master <- data.frame(
  param = c("ATE(1,0)", "NDE(1,0)", "NIE(1,0)", "CDE(1,0,0)"),
  est = c(
    out1$ATE,
    out1$NDE,
    out1$NIE,
    out1_cde$CDE
  )
)

# Open log
sink(log_path, split = TRUE)

master |>
  mutate(
    est = round(est, 3)
  )

# Close log
sink()
