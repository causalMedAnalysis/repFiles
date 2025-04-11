#----------Preliminaries----------#
rm(list = ls())
chapter <- "ch5"
title <- "table_5-9"

# Specify the root directory:
dir_root <- "Please Change to Your Local Directory" 

# Define subdirectories for logs and figures:
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)
dir_tab <- paste0(dir_root, "/table/", chapter)

# Ensure all necessary directories exist under your root folder，
# if not, the function will create folders for you:

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
create_dir_if_missing(dir_tab)

# Open log
sink(log_path, split = TRUE)

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch5/table_5-8.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/Brader_et_al2008/Brader_et_al2008.RData
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/linmed.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwmed.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/linpath.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwpath.R

# Outputs:     .../code/ch5/_LOGS/table_5-6_log.txt

# Description: Replicates Chapter 5, Table 5.6: Total and Path-Specific Effects 
#              of College Attendance on CES-D Scores as Estimated using Regression 
#.             Imputation with NLSY.
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL DEPENDENCIES and LOAD RERUIRED PACKAGES
#------------------------------------------------#

# The following packages are required for replicate results:
packages <-
  c(
    "tidyverse", 
    "paths"
  )

# Below function will automatically download the package you need,
# otherwise simply load the package:
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

# helper functions:
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/linmed.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwmed.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/linpath.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwpath.R")
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/pathimp.R")


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

# baseline confounder(s)
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

# Load the data:

temp_file <- tempfile() # define a placeholder to store the data

download.file(
  "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/Brader_et_al2008/Brader_et_al2008.RData", 
  temp_file, 
  mode = "wb")

load(temp_file)

# Process the data:
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
#                            REPLICATE TABLE 5.9                              #
#------------------------------------------------------------------------------#

#-----------------------------------------------------------------------#
#           Example 1: Regression Imputation Without Interaction:      #
#---------------------------------------------------------------------#

# Specify outcome models:
# E(Y|D,C):
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

# Combine all models into a list
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

# Fit the Paths Model:
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
    boot_parallel = "no",
    out_ipw = FALSE
  )$summary_df


#-----------------------------------------------------------------#
#       Example 2: Regression Imputation With D x M Interaction   #
#-----------------------------------------------------------------#

# Specify outcome models:
# E(Y|D,C):
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

# Combine all models into a list
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

# Fit the Paths Model:
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
    boot_parallel = "no",
    boot_seed = boot_seed,
    out_ipw = FALSE
  )$summary_df

#----------------------------------------------------------------------#
#       Example 3: Regression Imputation With D x {M, C} Interaction   #
#----------------------------------------------------------------------#

# Specify outcome models:
# E(Y|D,C):
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

# Combine all models into a list
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

# Fit the Paths Model:
Paths_DMCInteraction <-
  pathimp(
    D = D,
    Y = Y,
    M = M,
    Y_models = glm_ymodels,
    D_model = NULL,
    data = Brader,
    boot_reps = 2000,
    boot_parallel = "no",
    boot_seed = boot_seed,
    out_ipw = FALSE
  )$summary_df

#-----------------------------#
#   GENERATE FINAL TABLE      #
#-----------------------------#

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
write_csv(master, paste0(dir_tab,"/table5_9.csv"))












