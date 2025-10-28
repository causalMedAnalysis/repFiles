# Preliminaries
chapter <- "ch5"
title <- "table_5-1"
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

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
create_dir_if_missing(dir_fig)

# Open log
sink(log_path, split = TRUE)

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch5/table_5-1.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch5/_LOGS/table_5-1_log.txt

# Description: Replicates Chapter 5, Table 5.1: Total, Direct, and Indirect 
#              Effects of College Attendance on CES-D Scores as Estimated from 
#              the One-Mediator-at-a-Time Approach Applied to the NLSY.
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
  #M, # will handle specially below
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

# create two samples for the analysis data, taking complete cases of the 
# relevant variables for each analysis
nlsy_m1 <- nlsy_raw[complete.cases(nlsy_raw[,c(key_vars,M[1])]),] |>
  mutate(
    std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
  )

nlsy_m2 <- nlsy_raw[complete.cases(nlsy_raw[,c(key_vars,M[2])]),] |>
  mutate(
    std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
  )

#-------------------------#
#  MEDIATOR 1, VERSION 1  #
#-------------------------#
# Additive linear model

out_m1v1 <- linmed(
  data = nlsy_m1,
  D = D,
  M = M[1],
  Y = Y,
  C = C,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------------#
#  MEDIATOR 1, VERSION 2  #
#-------------------------#
# Linear model with D x M interaction

out_m1v2 <- linmed(
  data = nlsy_m1,
  D = D,
  M = M[1],
  Y = Y,
  C = C,
  interaction_DM = TRUE,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------------#
#  MEDIATOR 2, VERSION 1  #
#-------------------------#
# Additive linear model

out_m2v1 <- linmed(
  data = nlsy_m2,
  D = D,
  M = M[2],
  Y = Y,
  C = C,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------------#
#  MEDIATOR 2, VERSION 2  #
#-------------------------#
# Linear model with D x M interaction

out_m2v2 <- linmed(
  data = nlsy_m2,
  D = D,
  M = M[2],
  Y = Y,
  C = C,
  interaction_DM = TRUE,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------#
#  COLLATE RESULTS  #
#-------------------#
master <- data.frame(
  mediator = c("M1 (Unemployment)", "", "", "M2 (Household Income)", "", ""),
  param = c("ATE(1,0)", "NDE(1,0)", "NIE(1,0)"),
  
  # version 1: additive linear model
  v1_est = c(
    ## mediator 1
    out_m1v1$ATE,
    out_m1v1$NDE,
    out_m1v1$NIE,
    ## mediator 2
    out_m2v1$ATE,
    out_m2v1$NDE,
    out_m2v1$NIE
  ),
  v1_ci_low = c(
    ## mediator 1
    out_m1v1$ci_ATE[1],
    out_m1v1$ci_NDE[1],
    out_m1v1$ci_NIE[1],
    ## mediator 2
    out_m2v1$ci_ATE[1],
    out_m2v1$ci_NDE[1],
    out_m2v1$ci_NIE[1]
  ),
  v1_ci_high = c(
    ## mediator 1
    out_m1v1$ci_ATE[2],
    out_m1v1$ci_NDE[2],
    out_m1v1$ci_NIE[2],
    ## mediator 2
    out_m2v1$ci_ATE[2],
    out_m2v1$ci_NDE[2],
    out_m2v1$ci_NIE[2]
  ),
  
  # version 2: linear model with D x M interaction
  v2_est = c(
    ## mediator 1
    out_m1v2$ATE,
    out_m1v2$NDE,
    out_m1v2$NIE,
    ## mediator 2
    out_m2v2$ATE,
    out_m2v2$NDE,
    out_m2v2$NIE
  ),
  v2_ci_low = c(
    ## mediator 1
    out_m1v2$ci_ATE[1],
    out_m1v2$ci_NDE[1],
    out_m1v2$ci_NIE[1],
    ## mediator 2
    out_m2v2$ci_ATE[1],
    out_m2v2$ci_NDE[1],
    out_m2v2$ci_NIE[1]
  ),
  v2_ci_high = c(
    ## mediator 1
    out_m1v2$ci_ATE[2],
    out_m1v2$ci_NDE[2],
    out_m1v2$ci_NIE[2],
    ## mediator 2
    out_m2v2$ci_ATE[2],
    out_m2v2$ci_NDE[2],
    out_m2v2$ci_NIE[2]
  )
)

master |>
  mutate(
    across(
      .cols = !c(mediator,param),
      .fns = \(x) round(x, 3)
    )
  )

# Close log
sink()
