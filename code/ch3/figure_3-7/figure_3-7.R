# Preliminaries
chapter <- "ch3"
title <- "figure_3-7"

#dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"
dir_root <- "C:/Users/Geoffrey Wodtke/Desktop/repFiles-Dev"

dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

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
create_dir_if_missing(dir_fig)

# Open log
sink(log_path, split = TRUE)

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/figure_3-7.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch3/_LOGS/figure_3-7_log.txt
#              .../figures/ch3/figure_3-7.png

# Description: Replicates Chapter 3, Figure 3-7: Histogram of Bootstrap 
#              Estimates for NIE-hat(1,0)^ipw based on the NLSY.
#-------------------------------------------------------------------------------

#---------------------------------------------------------#
#  INSTALL DEPENDENCIES AND LOAD CAUSAL MED FUNCTIONS     #
#---------------------------------------------------------#
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

#install_github("causalMedAnalysis/causalMedR-Dev")
#library(causalMedR)

install.packages("C:/Users/Geoffrey Wodtke/Desktop/cmedR_0.1.0.tar.gz", repos = NULL, type = "source")
library(cmedR)

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
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
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
formula1_D <- as.formula(paste(D, "~", predictors1_D))

# D model 2 formula: s(D|C,M)
predictors2_D <- paste(c(M,C), collapse = " + ")
formula2_D <- as.formula(paste(D, "~", predictors2_D))

# M model formula: g(M|C,D)
predictors_M <- paste(c(D,C), collapse = " + ")
formula_M <- as.formula(paste(M, "~", predictors_M))

# Estimate ATE(1,0), NDE(1,0), NIE(1,0)
out1 <- ipwmed(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  D_C_model = formula1_D,
  D_CM_model = formula2_D,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-----------------#
#  CREATE FIGURE  #
#-----------------#
data.frame(NIE = out1$boot_NIE) |>
  ggplot(aes(x=NIE)) +
  geom_histogram(binwidth=0.00125, color="black", fill="grey90") +
  #geom_density(adjust=0.8, alpha=0.1, fill="black") +
  geom_vline(xintercept=out1$ci_NIE[1], linetype=2, linewidth=0.5) +
  geom_vline(xintercept=out1$ci_NIE[2], linetype=2, linewidth=0.5) +
  theme_bw() +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 160),
    breaks = seq(0, 160, 20)
  ) +
  scale_x_continuous(
    name = expression(hat(theta)[b]),
    limits = c(-0.05, 0.02),
    breaks = round(seq(-0.05, 0.02, 0.01), 2)
  )

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
