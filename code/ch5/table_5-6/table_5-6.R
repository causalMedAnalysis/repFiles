#----------Preliminaries----------#
rm(list = ls())
chapter <- "ch5"
title <- "table_5-6"

# Specify the root directory:
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text" 

# Define subdirectories for logs and figures:
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")

# Ensure all necessary directories exist under your root folder
# if not, the following function will create folders for you

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

# Script:      .../code/ch5/table_5-6.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/pathimp.R

# Outputs:     .../code/ch5/_LOGS/table_5-6_log.txt

# Description: Replicates Chapter 5, Table 5.6: Total and Path-Specific Effects 
#              of College Attendance on CES-D Scores as Estimated using Regression 
#              Imputation with NLSY.
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL DEPENDENCIES and LOAD RERUIRED PACKAGES
#------------------------------------------------#

# The following packages are required for replicate results:
packages <-
  c(
    "tidyverse", 
    "paths" # Main Package for this exercise
  )

# Function below will automatically download the packages you need
# Otherwise simply load the required packages
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

#-----------------------------#
#  LOAD CAUSAL MED FUNCTIONS  #
#-----------------------------#

source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/pathimp.R")

#------------------#
#  SPECIFICATIONS  #
#------------------#

# outcome
Y <- "std_cesd_age40"

# exposure
D <- "att22"

# mediators
M <- list(
  "ever_unemp_age3539",
  "log_faminc_adj_age3539"
)

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
  unlist(M),
  C
)

# number of bootstrap replications
n_reps <- 2000

# set seed:
boot_seed <- 02138

#-----------------------------#
#        PREPARE DATA         #
#-----------------------------#

nlsy_raw <- read_stata(
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
)

df <- 
  nlsy_raw[complete.cases(nlsy_raw[,key_vars]),] |>
  mutate(std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40))

#------------------------------------------------------------------------------#
#                            REPLICATE TABLE 5.6                               #
#------------------------------------------------------------------------------#

# Specify form of the models for the outcome and exposure

# E(Y|D,C):
glm_m0 <- glm(
  std_cesd_age40 ~ female + black + hispan + paredu + parprof + 
    parinc_prank + famsize + afqt3 + att22,
  data = df
)

# E(Y|D,C,M1):
glm_m1 <- glm(
  std_cesd_age40 ~ female + black + hispan + paredu + parprof + 
    parinc_prank + famsize + afqt3 + att22 + 
    ever_unemp_age3539, 
  data = df
)

# E(Y|D,C,M1,M2):
glm_m2 <- glm(
  std_cesd_age40 ~ female + black + hispan + paredu + parprof + 
    parinc_prank + famsize + afqt3 + att22 + 
    ever_unemp_age3539 + log_faminc_adj_age3539, 
  data = df
)

# E(D|C):
glm_ps <- glm(
  att22 ~ female + black + hispan + paredu + parprof + 
    parinc_prank + famsize + afqt3, 
  family = binomial("logit"),
  data = df
)

# Store outcome models in a list:
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

# Compute effect estimates
glm_paths <-
  pathimp(
    D = D,
    Y = Y,
    M = M,
    Y_models = glm_ymodels,
    D_model = glm_ps,
    data = df,
    boot_reps = 2000,
    boot_seed = boot_seed,
    boot_parallel = "multicore",
    round_decimal = 3,
    out_ipw = TRUE
  )

# Open log
sink(log_path, split = TRUE)

print(glm_paths$summary_df)

# Close log
sink()

