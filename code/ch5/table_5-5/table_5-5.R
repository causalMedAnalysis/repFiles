# Preliminaries
chapter <- "ch5"
title <- "table_5-5"
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")

# Ensure all necessary directories exist under your root folder
# If not, the function below will create folders for you

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

# Open log
sink(log_path, split = TRUE)

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch5/table_5-5.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch5/_LOGS/table_5-5_log.txt

# Description: Replicates Chapter 5, Table 5.5: Total and Path-Specific Effects 
#              of College Attendance on CES-D Scores as Estimated from Linear 
#              Models and Inverse Probability Weighting with the NLSY.
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

# mediators
M <- c(
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
  M,
  C
)

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

#-------------------------------------#
#  LINEAR MODEL ESTIMATOR: Version 1  #
#-------------------------------------#
# Additive linear model

out_lin1 <- linpath(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  C = C,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------------------------#
#  LINEAR MODEL ESTIMATOR: Version 2  #
#-------------------------------------#
# Linear model with D x M interaction

out_lin2 <- linpath(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  C = C,
  interaction_DM = TRUE,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-----------------#
#  IPW ESTIMATOR  #
#-----------------#
# Additive logit models

out_ipw <- ipwpath(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  C = C,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------#
#  COLLATE RESULTS  #
#-------------------#
master <- data.frame(
  param = c("ATE(1,0)", paste0("PSE_{", names(out_lin1$PSE), "}(1,0)")),
  
  # linear models: version 1
  lin1_est = c(
    out_lin1$ATE,
    out_lin1$PSE
  ),
  lin1_ci_low = c(
    out_lin1$ci_ATE[1],
    out_lin1$ci_PSE[,1]
  ),
  lin1_ci_high = c(
    out_lin1$ci_ATE[2],
    out_lin1$ci_PSE[,2]
  ),
  
  # linear models: version 2
  lin2_est = c(
    out_lin2$ATE,
    out_lin2$PSE
  ),
  lin2_ci_low = c(
    out_lin2$ci_ATE[1],
    out_lin2$ci_PSE[,1]
  ),
  lin2_ci_high = c(
    out_lin2$ci_ATE[2],
    out_lin2$ci_PSE[,2]
  ),
  
  # IPW
  ipw_est = c(
    out_ipw$ATE,
    out_ipw$PSE
  ),
  ipw_ci_low = c(
    out_ipw$ci_ATE[1],
    out_ipw$ci_PSE[,1]
  ),
  ipw_ci_high = c(
    out_ipw$ci_ATE[2],
    out_ipw$ci_PSE[,2]
  )
)
rownames(master) <- NULL

width_curr <- getOption("width")
options(width = 500)
master |>
  mutate(
    across(
      .cols = !param,
      .fns = \(x) round(x, 3)
    )
  )
options(width = width_curr)

# Close log
sink()
