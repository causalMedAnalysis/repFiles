#----------Preliminaries----------#
rm(list = ls())
chapter <- "ch5"
title <- "table_5-9"

# Specify the root directory:
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text" 

# Define subdirectories for logs and figures:
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

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch5/table_5-9.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/Brader_et_al2008/Brader_et_al2008.RData

# Outputs:     .../code/ch5/_LOGS/table_5-9_log.txt

# Description: Replicates Chapter 5, Table 5.9: Total and Path-Specific Effects 
#              of Negative Media Framing on Support for Immigration as Estimated 
#              from the Regression Imputation Approach.
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#
packages <-
  c(
    "tidyverse", 
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

# set seed:
boot_seed <- 02138

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
#                            REPLICATE TABLE 5.9                               #
#------------------------------------------------------------------------------#

#-------------------------------------------------------------#
#           Regression Imputation Without Interactions        #
#-------------------------------------------------------------#

# Specify the outcome models

# E(Y|D,C)
glm_m0 <- glm(
  immigr ~ ppage + female + hs + sc + ba + ppincimp + tone_eth, 
  data = Brader
)

# E(Y|D,C,M1)
glm_m1 <- glm(
  immigr ~ ppage + female + hs + sc + ba + ppincimp + tone_eth + p_harm, 
  data = Brader
)

# E(Y|D,C,M1,M2)
glm_m2 <- glm(
  immigr ~ ppage + female + hs + sc + ba + ppincimp + tone_eth + p_harm + emo, 
  data = Brader
)

# Combine outcome models into a list
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

# Compute effect estimates
Paths_NoInteraction <-
  pathimp(
    D = D,
    Y = Y,
    M = M,
    Y_models = glm_ymodels,
    D_model = NULL,
    data = Brader,
    boot_reps = 2000,
    boot_seed = boot_seed,
    out_ipw = FALSE
  )$summary_df

#---------------------------------------------------------------#
#       Regression Imputation With D x M Interactions           #
#---------------------------------------------------------------#

# Specify outcome models

# E(Y|D,C)
glm_m0 <- glm(
  immigr ~ ppage + female + hs + sc + ba + ppincimp + tone_eth, 
  data = Brader
)

# E(Y|D,C,M1)
glm_m1 <- glm(
  immigr ~ ppage + female + hs + sc + ba + ppincimp + tone_eth +
    p_harm + tone_eth * p_harm, 
  data = Brader
)

# E(Y|D,C,M1,M2)
glm_m2 <- glm(
  immigr ~ ppage + female + hs + sc + ba + ppincimp + tone_eth +
    p_harm + emo + tone_eth * p_harm + tone_eth * emo, 
  data = Brader
)

# Combine outcome models into a list
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

# Compute effect estimates
Paths_DMInteraction <-
  pathimp(
    D = D,
    Y = Y,
    M = M,
    Y_models = glm_ymodels,
    D_model = NULL,
    data = Brader,
    boot_reps = 2000,
    round_decimal = 3,
    boot_seed = boot_seed,
    out_ipw = FALSE
  )$summary_df

#----------------------------------------------------------------------#
#       Regression Imputation With D x {M, C} Interactions             #
#----------------------------------------------------------------------#

# Specify outcome models

# E(Y|D,C)
glm_m0 <- glm(
  immigr ~  ppage + female + hs + sc + ba + ppincimp + tone_eth + 
    tone_eth * ppage + tone_eth * female + tone_eth * hs + 
    tone_eth * sc + tone_eth * ba + tone_eth * ppincimp, 
  data = Brader
)

# E(Y|D,C,M1)
glm_m1 <- glm(
  immigr ~ ppage + female + hs + sc + ba + ppincimp + tone_eth + 
    p_harm + tone_eth * ppage + tone_eth * female + 
    tone_eth * hs + tone_eth * sc + tone_eth * ba + 
    tone_eth * ppincimp + tone_eth * p_harm, 
  data = Brader
)

# E(Y|D,C,M1,M2)
glm_m2 <- glm(
  immigr ~ ppage + female + hs + sc + ba + ppincimp + tone_eth + 
    p_harm + emo + tone_eth * ppage + tone_eth * female + 
    tone_eth * hs + tone_eth * sc + tone_eth * ba + 
    tone_eth * ppincimp + tone_eth * p_harm + tone_eth * emo,
  data = Brader
)

# Combine outcome models into a list
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

# Compute effect estimates
Paths_DMCInteraction <-
  pathimp(
    D = D,
    Y = Y,
    M = M,
    Y_models = glm_ymodels,
    D_model = NULL,
    data = Brader,
    boot_reps = 2000,
    boot_seed = boot_seed,
    out_ipw = FALSE
  )$summary_df

#------------------------#
#   COLLATE RESULTS      #
#------------------------#

master <-
  reduce(
    list(
      Paths_NoInteraction %>% dplyr::select(-estimator),
      Paths_DMInteraction %>% dplyr::select(-estimator),
      Paths_DMCInteraction %>% dplyr::select(-estimator)
    ),
    left_join,
    by = "estimand"
  )

colnames(master) <- c(
  "estimand",
  "Outcome Models without Interactions",
  "Outcome Models with D x M Interactions",
  "Outcome Models with D x {M, C} Interactions"
)

# Open log
sink(log_path, split = TRUE)

print(master)

# Close log
sink()
