# Preliminaries
chapter <- "ch3"
title <- "figure_3-8"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/figure_3-8.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwmed.R

# Outputs:     .../code/ch3/_LOGS/figure_3-8_log.txt
#              .../figures/ch3/figure_3-8.png

# Description: Replicates Chapter 3, Figure 3-8: Kernel Density Plot of 
#              Bootstrap Estimates for NIE-hat(1,0)^ipw based on the NLSY.
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
# IPW estimator
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwmed.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/ipwmed.R")




#------------------#
#  SPECIFICATIONS  #
#------------------#
# outcome
Y <- "std_cesd_age40"

# exposure
D <- "att22"

# mediator
M <- "ever_unemp_age3539"

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
  C
)

# number of bootstrap replications
n_reps <- 2000




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
  # ^ Note that parallelizing the bootstrap is optional, but requires that you 
  # have installed the following R packages: doParallel, doRNG, foreach.
  # (You do not need to load those packages beforehand, with the library 
  # function.)
  # If you choose not to parallelize the bootstrap (by setting the boot_parallel 
  # argument to FALSE), the results may differ slightly, due to simulation 
  # variance (even if you specify the same seed).
)




#-----------------#
#  CREATE FIGURE  #
#-----------------#
dens <- data.frame(NIE = out1$boot_NIE) |>
  ggplot(aes(x=NIE)) +
  geom_density(adjust=0.9, alpha=0.1, fill="black") +
  theme_bw() +
  scale_y_continuous(
    name = "Density",
    limits = c(0, 60),
    breaks = seq(0, 60, 10)
  ) +
  scale_x_continuous(
    name = expression(hat(theta)[b]),
    limits = c(-0.05, 0.02),
    breaks = round(seq(-0.05, 0.02, 0.01), 2)
  )

dens_data <- ggplot_build(dens)$data[[1]]

dens <- dens + 
  geom_area(
    data = subset(dens_data, x>0),
    aes(x=x, y=y),
    fill = "grey30"
  )

print(dens)

# Save the figure
ggsave(
  paste0(dir_fig, "/", title, ".png"),
  height = 4.5,
  width = 4.5,
  units = "in",
  dpi = 600
)


# Close log
sink()

