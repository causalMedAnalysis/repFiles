#----------Preliminaries----------#
rm(list = ls())
chapter <- "ch5"
title <- "table_5-8"

# Specify the root directory:
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text" 

# Define subdirectories for logs and figures:
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

# Script:      .../code/ch5/table_5-8.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/Brader_et_al2008/Brader_et_al2008.RData

# Outputs:     .../code/ch5/_LOGS/table_5-8_log.txt

# Description: Replicates Chapter 5, Table 5.8: Total and Path-Specific Effects 
#              of Negative Media Framing on Support for Immigration as Estimated from 
#              Linear Models and Inverse Probability Weighting.
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#
packages <-
  c(
    "tidyverse", 
    "margins", 
    "mediation", 
    "foreach", 
    "doParallel", 
    "doRNG",
    "paths",
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
Y <- "immigr"

# exposure
D <- "tone_eth"

# mediators
M1 <- "p_harm"
M2 <- "emo"
M <- 
  list(
    M1,
    M2
  )

# baseline confounders
C <- c(
  "ppage", 
  "female", 
  "hs", 
  "sc", 
  "ba", 
  "ppincimp"
  )

# key variables
key_vars <- c(
  "immigr", # unstandardized version of Y
  D,
  unlist(M),
  C
  )

# number of bootstrap replications
nboot <- 2000

# set seed
boot_seed <- 3308004

#-----------------------------#
#        PREPARE DATA         #
#-----------------------------#

# Load the data
temp_file <- tempfile() # define a placeholder to store the data

download.file(
  "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/Brader_et_al2008/Brader_et_al2008.RData", 
  temp_file, 
  mode = "wb")

load(temp_file)

# Process the data
Brader <- 
  Brader %>%
  dplyr::select(
    immigr, 
    emo, 
    p_harm, 
    tone_eth, 
    ppage, 
    ppeducat, 
    ppgender, 
    ppincimp) %>% 
  na.omit() %>%
  mutate(
    immigr = as.numeric(scale(4 - immigr)),
    hs = (ppeducat == "high school"),
    sc = (ppeducat == "some college"),
    ba = (ppeducat == "bachelor's degree or higher"),
    female = (ppgender == "female")) %>%
  mutate_at(
    vars(emo, p_harm, ppage, female, hs, sc, ba, ppincimp), 
    ~ . - mean(., na.rm = TRUE)  
  )

#------------------------------------------------------------------------------#
#                            REPLICATE TABLE 5.8.                              #
#------------------------------------------------------------------------------#

#-------------------------------------------------------#
#     Linear Models Without D x M Interactions:         #
#-------------------------------------------------------#

LinMod <-
  linpath(
    data = Brader,
    D = D,
    M = M,
    Y = Y,
    C = C,
    boot = TRUE,
    boot_reps = 2000,
    boot_parallel = TRUE,
    boot_seed = boot_seed
  )

#---------------------------------------------------#
#     Linear Models With D x M Interactions:        #
#---------------------------------------------------#

LinModX <-
  linpath(
    data = Brader,
    D = D,
    M = M,
    Y = Y,
    C = C,
    interaction_DM = TRUE,
    boot = TRUE,
    boot_reps = 2000,
    boot_parallel = TRUE,
    boot_seed = boot_seed
  )

#--------------------------------------------#
#         Inverse Probability Weighting      #
#--------------------------------------------#

IPW <- 
  ipwpath(
    data = Brader,
    D = D,
    M = M,
    Y = Y,
    C = C,
    boot = TRUE,
    boot_reps = 2000,
    boot_parallel = TRUE,
    boot_seed = boot_seed
  )

#-------------------------------#
#        COLLATE RESULTS        #
#-------------------------------#

# Create the master dataframe
master <- data.frame(
  param = c("ATE(1,0)", paste0("PSE_{", names(LinMod$PSE), "}(1,0)")),
  
  LinMod_est     = c(LinMod$ATE, LinMod$PSE),
  LinMod_ci_low  = c(LinMod$ci_ATE[1], LinMod$ci_PSE[, 1]),
  LinMod_ci_high = c(LinMod$ci_ATE[2], LinMod$ci_PSE[, 2]),
  
  LinModX_est     = c(LinModX$ATE, LinModX$PSE),
  LinModX_ci_low  = c(LinModX$ci_ATE[1], LinModX$ci_PSE[, 1]),
  LinModX_ci_high = c(LinModX$ci_ATE[2], LinModX$ci_PSE[, 2]),
  
  IPW_est     = c(IPW$ATE, IPW$PSE),
  IPW_ci_low  = c(IPW$ci_ATE[1], IPW$ci_PSE[, 1]),
  IPW_ci_high = c(IPW$ci_ATE[2], IPW$ci_PSE[, 2])
  ) %>%
  pivot_longer(
    cols      = -param, 
    names_to  = "model", 
    values_to = "value"
  ) %>%
  separate(
    model, 
    into  = c("prefix", "stat"), 
    sep   = "_", 
    extra = "merge"
  ) %>%
  pivot_wider(
    names_from  = stat, 
    values_from = value
  ) %>%
  mutate(
    final_value = paste0(
      sprintf("%.3f", est), 
      " [", sprintf("%.3f", ci_low), ", ", sprintf("%.3f", ci_high), "]"
    )
  ) %>%
  dplyr::select(param, prefix, final_value) %>%
  pivot_wider(
    names_from  = prefix, 
    values_from = final_value
  )

colnames(master) <- c(
  "Estimand",
  "Additive Linear Model(LinMod)",
  "LinMod with D x M Interactions",
  "Inverse Probability Weighting"
  )

# Open log
sink(log_path, split = TRUE)

# Print table
width_curr <- getOption("width")
options(width = 300)
print(master)

# Close log
sink()
