# Preliminaries
chapter <- "ch4"
title <- "table_4-5"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch4/table_4-5.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwvent.R

# Outputs:     .../code/ch4/_LOGS/table_4-5_log.txt

# Description: Replicates Chapter 4, Table 4.5: Bootstrap Inferential Statistics 
#              for the Interventional Effects of College Attendance on CES-D 
#              Scores Computed from the NLSY.
#-------------------------------------------------------------------------------


#-------------#
#  LIBRARIES  #
#-------------#
library(tidyverse)
library(haven)




#-----------------------------#
#  LOAD CAUSAL MED FUNCTIONS  #
#-----------------------------#
# utilities
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/utils_bare.R")
# RWR estimator
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/rwrlite.R")
# simulation estimator
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/medsim.R")
# IPW estimator
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwvent.R")
#source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/ipwvent.R")




#------------------#
#  SPECIFICATIONS  #
#------------------#
# outcome
Y <- "std_cesd_age40"

# exposure
D <- "att22"

# mediator
M <- "log_faminc_adj_age3539"

# exposure-induced confounder
L <- "ever_unemp_age3539"

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
  L,
  C
)

# mediator value for CDE
m <- log(5e4)

# number of simulations for simulation estimator
n_sims <- 2000

# number of bootstrap replications
#n_reps <- 2000
n_reps <- 5




#----------------#
#  PREPARE DATA  #
#----------------#
nlsy_raw <- read_stata(
  #file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
  file = "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Data/NLSY79/nlsy79BK_ed2.dta"
)

nlsy <- nlsy_raw[complete.cases(nlsy_raw[,key_vars]),] |>
  mutate(
    std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
  )




#------------------#
#  MODEL FORMULAE  #
#------------------#
# L model: With D x C interactions
# M model: With D x C interactions
# Y model: With D x M, D x C, M x C, M x L interactions
# Note that M here is the log form

# L and M model formulae
## main effects
predictors_LM <- paste(c(D,C), collapse = " + ")
## D x C interactions
predictors_LM <- paste(
  predictors_LM,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## full formula
(formula_L_string <- paste(L, "~", predictors_LM))
(formula_M_string <- paste(M, "~", predictors_LM))
formula_L <- as.formula(formula_L_string)
formula_M <- as.formula(formula_M_string)

# Y model formula
## main effects
predictors_Y <- paste(c(D,M,L,C), collapse = " + ")
## D x M interaction
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(D, M, sep = ":", collapse = " + ")
)
## D x C interactions
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## M x C interactions
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(M, C, sep = ":", collapse = " + ")
)
## M x L interaction
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(M, L, sep = ":", collapse = " + ")
)
## full formula
(formula_Y_string <- paste(Y, "~", predictors_Y))
formula_Y <- as.formula(formula_Y_string)




#-----------------#
#  RWR ESTIMATOR  #
#-----------------#
out_rwr <- rwrlite(
  data = nlsy,
  D = D,
  C = C,
  m = m,
  Y_formula = formula_Y,
  M_formula = formula_M,
  L_formula_list = list(formula_L),
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
  # ^ Note that parallelizing the bootstrap is optional, but requires that you 
  # have installed the following R packages: doParallel, doRNG, foreach.
  # The rwrlite() function also requires that you have installed the rwrmed R 
  # package.
  # (You do not need to load any of these packages beforehand, with the library 
  # function.)
  # If you choose not to parallelize the bootstrap (by setting the boot_parallel 
  # argument to FALSE), the results may differ slightly, due to simulation 
  # variance (even if you specify the same seed).
)




#------------------------#
#  SIMULATION ESTIMATOR  #
#------------------------#
# Define model specifications
out_sim_specs <- list(
  ## L model
  list(
    func = "glm",
    formula = formula_L,
    args = list(family = "binomial")
  ),
  ## M model
  list(
    func = "lm",
    formula = formula_M
  ),
  ## Y model
  list(
    func = "lm",
    formula = formula_Y
  )
)

# Estimate interventional effects
out_sim <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = M,
  model_spec = out_sim_specs,
  boot = TRUE,
  reps = n_reps,
  seed = 3308004
  # ^ Because of its resource intensity, running a non-parallelized bootstrap 
  # of the simulation estimator is not advisable. Therefore, unlike the other 
  # estimator functions, the parallelized bootstrap is the only implemented 
  # bootstrap in the medsim() function. If you request a bootstrap (as in the 
  # code above), you must have installed the following R packages: 
  # doParallel, doRNG, foreach.
  # (You do not need to load those packages beforehand, with the library 
  # function.)
)

# Estimate CDE(1,0,ln(50K))
out_sim_cde <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = paste0(M,"=",m),
  model_spec = out_sim_specs,
  boot = TRUE,
  reps = n_reps,
  seed = 3308004
)




#-----------------#
#  IPW ESTIMATOR  #
#-----------------#
# to do




#-------------------#
#  COLLATE RESULTS  #
#-------------------#
master <- data.frame(
  param = c("OE(1,0)", "IDE(1,0)", "NIE(1,0)", "CDE(1,0,ln(50K))"),
  
  # RWR
  rwr_pvalue = c(
    out_rwr$pvalue_OE,
    out_rwr$pvalue_IDE,
    out_rwr$pvalue_IIE,
    out_rwr$pvalue_CDE
  ),
  rwr_ci_low = c(
    out_rwr$ci_OE[1],
    out_rwr$ci_IDE[1],
    out_rwr$ci_IIE[1],
    out_rwr$ci_CDE[1]
  ),
  rwr_ci_high = c(
    out_rwr$ci_OE[2],
    out_rwr$ci_IDE[2],
    out_rwr$ci_IIE[2],
    out_rwr$ci_CDE[2]
  ),
  
  # simulation
  sim_pvalue = c(
    out_sim$pval[3],
    out_sim$pval[1],
    out_sim$pval[2],
    out_sim_cde$pval[1]
  ),
  sim_ci_low = c(
    out_sim$ll.95ci[3],
    out_sim$ll.95ci[1],
    out_sim$ll.95ci[2],
    out_sim_cde$ll.95ci[1]
  ),
  sim_ci_high = c(
    out_sim$ul.95ci[3],
    out_sim$ul.95ci[1],
    out_sim$ul.95ci[2],
    out_sim_cde$ul.95ci[1]
  ),
  
  # IPW
  ipw_pvalue = c(
    NA_real_,
    NA_real_,
    NA_real_,
    NA_real_
  ),
  ipw_ci_low = c(
    NA_real_,
    NA_real_,
    NA_real_,
    NA_real_
  ),
  ipw_ci_high = c(
    NA_real_,
    NA_real_,
    NA_real_,
    NA_real_
  )
)

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

