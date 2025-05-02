#----------Preliminaries---------#
rm(list=ls())
chapter <- "ch5"
title <- "figure_5-6"

# Specify the root directory:
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text" 

# Define subdirectories for logs and figures:
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

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch5/table_5-9.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/pathimp.R

# Outputs:     .../code/ch5/_LOGS/figure_5-6_log.txt
#              .../figures/ch5/figure_5-6.png

# Description: Replicates Chapter 5, Figure 5.6: Bias-adjusted Estimates of the 
#              Direct Effect of Issue Framing on Support for Immigration
#-------------------------------------------------------------------------------


#----------------------------------------------------#
#  INSTALL DEPENDENCIES and LOAD RERUIRED PACKAGES   #
#----------------------------------------------------#

# The following packages are required to replicate results
packages <-
  c(
    "tidyverse", 
    "paths",
    "latex2exp"
  )

# Function below will automatically download the packages you need
# Otherwise simply load the packages

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
  "ppincimp",
  "female_sens",
  "ba_sens"
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
  ) %>%
  mutate(
    female_sens = as.double(ppgender == "female"),
    ba_sens = as.double(ppeducat == "bachelor's degree or higher")
  )


#------------------------------------------------------------------------------#
#                            REPLICATE FIGURE 5.6                              #
#------------------------------------------------------------------------------#

#-------------------------------------------------------------------#
#   Step 1: Compute Effect Estimates using Regression Imputation    #
#-------------------------------------------------------------------#

# Specify outcome models

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

# Compute estimates
Paths_Model <-
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
  )

#----------------------------------------------#
#   Step 2: Perform Sensitivity Analysis       #
#----------------------------------------------#

Paths_Sens_Model <- sens(
  Paths_Model$org_obj, 
  confounded = "M1", 
  estimand = "direct",
  gamma_values = seq(0, 2, 0.005),
  eta_values = seq(-0.5, 0.5, 0.005)
  )

#-----------------------------------#
#   Step 3: Add Reference Points    #
#-----------------------------------#

# Model_UY: E(Y | C, D, M1, M2, U)
Model_UY <- lm(
  immigr ~ ppage + female_sens + hs + sc + 
    ba_sens + ppincimp + tone_eth + p_harm + emo, 
  data = Brader
  )

# Model_UD_fem: E(U | D, M1, M2, U1)
Model_UD_fem <- lm(
  female_sens ~ ppage + hs + sc + ba_sens + 
    ppincimp + tone_eth + p_harm + emo, 
  data = Brader
  )

# Model_UD_ba: E(U | D, M1, M2, U2)
Model_UD_ba <- lm(
  ba_sens ~ ppage + hs + sc + female_sens + 
    ppincimp + tone_eth + p_harm + emo, 
  data = Brader
  )

delta_UY_fem <- Model_UY$coefficients["female_sens"] # effect of U on Y
delta_DU_fem <- Model_UD_fem$coefficients["tone_eth"] # effect of U on D

delta_UY_ba <- Model_UY$coefficients["ba_sens"] # effect of U on Y
delta_DU_ba <- Model_UD_ba$coefficients["tone_eth"] # effect of U on D

#--------------------------------------#
#   Step 4: Generate Contour Plot      #
#-------------------------------- -----#

plot(
  Paths_Sens_Model, 
  outcome_name = "Support for Immigration") +
  xlab(bquote("Conditional Difference in the Prevalence of U across levels of D"~(delta[DU]))) +
  ylab(bquote("Conditional Mean Difference in Y across levels of U"~(delta[UY]))) +
  annotate("point", x = delta_DU_fem, y = delta_UY_fem) +
  annotate("text", x = delta_DU_fem + 0.05, y = delta_UY_fem, label = "female", size = 5) +
  annotate("point", x = delta_DU_ba, y = delta_UY_ba) +
  annotate("text", x = delta_DU_ba - 0.05, y = delta_UY_ba, label = "college\n graduate", size = 5) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c(NA, "grey70"),  na.value = NA)

ggsave(paste0(dir_fig,"/figure_5-6.png"), width = 10, height = 7)

