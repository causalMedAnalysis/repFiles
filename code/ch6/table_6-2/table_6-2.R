#----------Preliminaries----------#
rm(list = ls())
chapter <- "ch6"
title <- "table_6-2"

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

# Script:      .../code/ch5/table_6-2.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.RDS
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/mrmed.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/dmlmed.R

# Outputs:     .../code/ch6/_LOGS/table_6-2_log.txt

# Description: Replicates Chapter 6, Table 6.2: Multiply Robust Estimates of the Total, 
#                     Natural Direct, and Natural Indirect Effects of College Attendance 
#                     on CES-D scores, as Mediated by Unemployment, from the NLSY
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL DEPENDENCIES AND LOAD RERUIRED PACKAGES
#-------------------------------------------------#

# The following packages are required for replicate results:
packages <-
  c(
    "survey", 
    "gbm",
    "ranger",
    "glmnet",
    "rsample",
    "caret",
    "rlang",
    "tidyverse",
    "Hmisc",
    "SuperLearner",
    "scales",
    "haven"
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
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/mrmed.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/dmlmed.R")

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
  "cesd_age40", 
  D,
  unlist(M),
  C
)

# number of bootstrap replications
n_reps <- 2000

# set seed
seed <- 02138

#-----------------------------#
#        PREPARE DATA         #
#-----------------------------#

nlsy_raw <- as.data.frame(
  readRDS(
    url("https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.RDS")
    )
  )

df <- nlsy_raw[complete.cases(nlsy_raw[, c(D, "cesd_age40", unlist(M), C)]),] %>%
  mutate(std_cesd_age40 = as.numeric(scale(cesd_age40)))

#==============================================================================#
#                                   TABLE 6.2                                  #
#==============================================================================#

#==============================================================================#
#                                   Method 1                                   #
#==============================================================================#

#--------------------------------------------------#
#  Set Up the Exposure and Outcome Models          #
#--------------------------------------------------#

# Exposure ~ Baseline Confounders
D_C_model <- paste(D, " ~ ", paste(C, collapse= "+"))

# Mediator ~ Treatment + Baseline Confounders
M_DC_model <- paste(M[[1]], " ~ ", paste(c(C,D), collapse= "+"))

# Outcome ~  Treatment + Baseline Confounders + Mediators
Y_DMC_model <- paste(Y, " ~ ", paste(c(C,D,M[[1]]), collapse= "+"))

#--------------------------------------------------#
# Column 1: Parametric MR Estimation               #
#--------------------------------------------------#

mrmed1_rst <-
  mrmed(
    D = D,
    Y = Y,
    M = M[[1]],
    D_C_model = D_C_model,
    Y_DMC_model = Y_DMC_model,
    M_DC_model = M_DC_model,
    data = df,
    d = 1,
    dstar = 0,
    boot = TRUE,
    boot_reps = 2000,
    boot_seed = seed
  )

#--------------------------------------------------#
# Column 2: DML Estimation                         #
#--------------------------------------------------#

dmlmed1_rst <-
  dmlmed(
    D = D,
    Y = Y,
    M = M[[1]],
    D_C_model = D_C_model,
    Y_DMC_model = Y_DMC_model,
    M_DC_model = M_DC_model,
    data = df,
    d = 1,
    dstar = 0,
    seed = seed,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger")
  )

#------------------------------------------------------------------------------#
#                                    Method 2                                  #
#------------------------------------------------------------------------------#

#--------------------------------------------------#
#  Set Up the Exposure and Outcome Models          #
#--------------------------------------------------#

# Exposure ~ Baseline Confounders
D_C_model <- as.formula(paste(D, " ~ ", paste(C, collapse= "+")))

# Exposure ~ Baseline Confounders and Mediator
D_MC_model <- as.formula(paste(D, " ~ ", paste(c(C, M[[1]]), collapse= "+")))

# Outcome ~ Baseline Confounders and Exposure
Y_DC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D), collapse= "+")))

#Outcome ~ Baseline Confounders, Exposure and Mediator
Y_DMC_model <- as.formula(paste(Y, " ~ ", paste(c(C, D, M[[1]]), collapse= "+")))

#--------------------------------------------------#
# Column 3 : Parametric MR Estimation              #
#--------------------------------------------------#

mrmed2_rst <-
  mrmed(
    D = D,
    Y = Y,
    M = M[[1]],
    D_C_model = D_C_model,
    D_MC_model = D_MC_model,
    Y_DC_model = Y_DC_model,
    Y_DMC_model = Y_DMC_model,
    data = df,
    d = 1,
    dstar = 0,
    boot = TRUE,
    boot_reps = 2000,
    boot_seed = seed
  )

#--------------------------------------------------#
# Column 4 : DML Estimation                        #
#--------------------------------------------------#

dmlmed2_rst <-
  dmlmed(
    D = D,
    Y = Y,
    M = M[[1]],
    D_C_model = D_C_model,
    D_MC_model = D_MC_model,
    Y_DC_model = Y_DC_model,
    Y_DMC_model = Y_DMC_model,
    data = df,
    d = 1,
    dstar = 0,
    seed = seed,
    SL.library = c("SL.mean", "SL.glmnet","SL.ranger")
  )

#------------------------------------------------------------------------------#
#                               Collate Results                                #
#------------------------------------------------------------------------------#

make_CI <- function(est, se) {
  paste0(
    round(est, 3), " [",
    round(est - 1.96 * se, 3), ", ",
    round(est + 1.96 * se, 3), "]"
  )
}

final_rst <- 
  tibble::tibble(
  estimand = c("ATE(1,0)", "NDE(1,0)", "NIE(1,0)"),
  mrmed1 = c(
    make_CI(mrmed1_rst$est1$`ATE(1,0)`, mrmed1_rst$boot_rst_lst$est1$sd_ATE),
    make_CI(mrmed1_rst$est1$`NDE(1,0)`, mrmed1_rst$boot_rst_lst$est1$sd_NDE),
    make_CI(mrmed1_rst$est1$`NIE(1,0)`, mrmed1_rst$boot_rst_lst$est1$sd_NIE)
  ),
  mrmed2 = c(
    make_CI(mrmed2_rst$est2$`ATE(1,0)`, mrmed2_rst$boot_rst_lst$est2$sd_ATE),
    make_CI(mrmed2_rst$est2$`NDE(1,0)`, mrmed2_rst$boot_rst_lst$est2$sd_NDE),
    make_CI(mrmed2_rst$est2$`NIE(1,0)`, mrmed2_rst$boot_rst_lst$est2$sd_NIE)
  )
 ) %>%
 left_join(
   dmlmed1_rst$est1 %>%
     dplyr::select(
       estimand,
       out
     ) %>%
    rename(
      `DML mrmed1` = out
    ),
   by = "estimand"
 ) %>%
  left_join(
    dmlmed2_rst$est2 %>%
      dplyr::select(
        estimand,
        out
      ) %>%
      rename(
        `DML mrmed2` = out
      ),
    by = "estimand"
  ) %>%
  dplyr::select(
    estimand,
    mrmed1,
    `DML mrmed1`,
    mrmed2,
    `DML mrmed2`
  )

# Open log
sink(log_path, split = TRUE)

# Print table
width_curr <- getOption("width")
options(width = 300)
print(final_rst)

# Close log
sink()
