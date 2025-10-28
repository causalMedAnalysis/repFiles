# Preliminaries
chapter <- "ch3"
title <- "table_3-5"
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

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

# Script:      .../code/ch3/table_3-5.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch3/_LOGS/table_3-5_log.txt

# Description: Replicates Chapter 3, Table 3-5: Inferential Statistics for 
#              Total, Direct, and Indirect Effects of College Attendance on 
#              CES-D Scores Computed from the NLSY using Inverse Probability 
#              Weighting and the Nonparametric Bootstrap.
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#
packages <-
  c(
    "tidyverse",
    "haven",
    "doParallel",
    "doRNG", 
    "foreach",
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

# baseline confounders
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

# number of bootstrap replications
n_reps <- 2000

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

#----------------------------------------#
#  ESTIMATE EFFECTS & PERFORM BOOTSTRAP  #
#----------------------------------------#
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
predictors_M <- paste(c(D,C), collapse = " + ")
formula_M_string <- paste(M, "~", predictors_M)
formula_M_string

# Estimate ATE(1,0), NDE(1,0), NIE(1,0)
out1 <- ipwmed(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  formula1_string = formula1_D_string,
  formula2_string = formula2_D_string,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

# Estimate CDE(1,0,0)
out1_cde <- ipwcde(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  m = m,
  formula_D_string = formula1_D_string,
  formula_M_string = formula_M_string,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------#
#  COLLATE RESULTS  #
#-------------------#
master <- data.frame(
  param = c("ATE(1,0)", "NDE(1,0)", "NIE(1,0)", "CDE(1,0,0)"),
  ci_low = c(
    out1$ci_ATE[1],
    out1$ci_NDE[1],
    out1$ci_NIE[1],
    out1_cde$ci_CDE[1]
  ),
  ci_high = c(
    out1$ci_ATE[2],
    out1$ci_NDE[2],
    out1$ci_NIE[2],
    out1_cde$ci_CDE[2]
  ),
  pvalue = c(
    out1$pvalue_ATE,
    out1$pvalue_NDE,
    out1$pvalue_NIE,
    out1_cde$pvalue_CDE
  )
)

# Open log
sink(log_path, split = TRUE)

master |>
  mutate(
    across(
      .cols = !param,
      .fns = \(x) round(x, 3)
    )
  )

# Close log
sink()
