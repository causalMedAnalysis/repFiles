# Preliminaries
chapter <- "ch4"
title <- "figure_4-10"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch4/figure_4-10.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/plowUse/plowUse.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R

# Outputs:     .../code/ch4/_LOGS/figure_4-10_log.txt

# Description: Replicates Chapter 4, Figure 4.10: Estimates of the 
#              Interventional Direct Effect and Overall Effect Adjusted for Bias 
#              Due to Unobserved Exposure-Outcome Confounding.
#-------------------------------------------------------------------------------


#------------------------#
#  INSTALL DEPENDENCIES  #
#------------------------#
# First, install dependencies available on CRAN.
# The following packages are used to create Figure 4.10.
dependencies_cran <- c("gridExtra", "metR")

#install.packages(dependencies_cran)
# ^ Uncomment this line above to install these packages.


# Second, install the rwrmed R package, which is available to install from 
# GitHub.
# To install the package directly, you must first have installed the devtools 
# package (which is available on CRAN).

#install.packages("devtools")
# ^ Uncomment this line above to install the devtools package, if you have not 
# already done so.

#devtools::install_github("xiangzhou09/rwrmed")
# ^ Uncomment this line above to install the rwrmed package from GitHub.




#-------------#
#  LIBRARIES  #
#-------------#
library(gridExtra)
library(metR)
library(tidyverse)
library(haven)




#-----------------------------#
#  LOAD CAUSAL MED FUNCTIONS  #
#-----------------------------#
# utilities
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")
# RWR estimator
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R")
# ^ Note that rwrlite() is a wrapper for two functions from the rwrmed package. 
# It requires that you have installed rwrmed. (But you do not need to load the 
# rwrmed package beforehand, with the library function.)




#------------------#
#  SPECIFICATIONS  #
#------------------#
# outcome
Y <- "women_politics"

# exposure
D <- "plow"

# mediator
M <- "ln_income"

# exposure-induced confounder
L <- "authGovCat"

# baseline confounder(s)
C <- c(
  "agricultural_suitability",
  "tropical_climate",
  "large_animals",
  "rugged"
)

# key variables
key_vars <- c(
  Y,
  D,
  M,
  "polity2_2000", # source variable for L
  C,
  "isocode" # observation/country identifier
)

# mediator value for CDE
m <- 7.5 # roughly equal to log(1800)




#----------------#
#  PREPARE DATA  #
#----------------#
plow_raw <- read_stata(
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/plowUse/plowUse.dta"
)

plow <- na.omit(plow_raw[,key_vars]) |>
  mutate(
    plow = round(plow),
    women_politics = women_politics / 100,
    authGovCat = case_when(
      polity2_2000>=6   ~ 3,
      polity2_2000>=-5  ~ 2,
      polity2_2000>=-10 ~ 1
    )
  )




#------------------#
#  MODEL FORMULAE  #
#------------------#
# L model: Additive
# M model: Additive
# Y model: With D x M interaction

# L and M model formulae
predictors_LM <- paste(c(D,C), collapse = " + ")
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
## full formula
(formula_Y_string <- paste(Y, "~", predictors_Y))
formula_Y <- as.formula(formula_Y_string)




#--------------------#
#  ESTIMATE EFFECTS  #
#--------------------#
out <- rwrlite(
  data = plow,
  D = D,
  C = C,
  m = m,
  Y_formula = formula_Y,
  M_formula = formula_M,
  L_formula_list = list(formula_L)
)




#---------------------------#
#  BIAS-ADJUSTED ESTIMATES  #
#---------------------------#
# Specify range for sensitivity parameters under D-Y confounding
sens_grid <- expand.grid(
  delta_UYgivCDLM = seq(-0.1, 0.1, 0.01),
  delta_DUgivC = seq(-0.1, 0.1, 0.01)
)

# Compute bias-adjusted estimates
adj_grid <- cbind(
  sens_grid,
  ide_adj = out$IDE - (sens_grid$delta_UYgivCDLM * sens_grid$delta_DUgivC),
  oe_adj = out$OE - (sens_grid$delta_UYgivCDLM * sens_grid$delta_DUgivC)
)




#------------------------#
#  CREATE CONTOUR PLOTS  #
#------------------------#
# Create contour plots of bias-adjusted estimates against sensitivity parameters

# Common specifications
common_specs <- list(
  scale_colour_distiller(palette="Greys", direction=1),
  xlab(expression(delta["UY|C,D,L,M"])),
  ylab(expression(delta["DU|C"])),
  scale_x_continuous(breaks=seq(-0.1, 0.1, 0.02)),
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.02)),
  theme_bw(base_size=11),
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
)

# IDE plot
plot_ide <- ggplot(
  adj_grid,
  aes(x=delta_UYgivCDLM, y=delta_DUgivC, z=ide_adj, color=after_stat(level))
) +  
  geom_contour(
    breaks = seq(
      round(min(adj_grid$ide_adj), 3),
      round(max(adj_grid$ide_adj), 3),
      0.002
    ),
    show.legend = FALSE
  ) +
  ggtitle("A. Bias-adjusted IDE(1,0) Estimates") +
  geom_text_contour(
    breaks = seq(
      round(min(adj_grid$ide_adj), 3),
      round(max(adj_grid$ide_adj), 3),
      0.002
    ),
    stroke = 0.3,
    size = 3,
    skip = 0,
    color = "black"
  ) +
  common_specs

# OE plot
plot_oe <- ggplot(
  adj_grid,
  aes(x=delta_UYgivCDLM, y=delta_DUgivC, z=oe_adj, color=after_stat(level))
) +  
  geom_contour(
    breaks = seq(
      round(min(adj_grid$oe_adj), 3),
      round(max(adj_grid$oe_adj), 3),
      0.002
    ),
    show.legend = FALSE
  ) +
  ggtitle("B. Bias-adjusted OE(1,0) Estimates") +
  geom_text_contour(
    breaks = seq(
      round(min(adj_grid$oe_adj), 3),
      round(max(adj_grid$oe_adj), 3),
      0.002
    ),
    stroke = 0.3,
    size = 3,
    skip = 0,
    color = "black"
  ) +
  common_specs

# Combine the plots
plot_comb <- grid.arrange(plot_ide, plot_oe, ncol=1)

# Save the combined plot
ggsave(
  paste0(dir_fig, "/", title, ".png"),
  plot = plot_comb,
  height = 8,
  width = 4.5,
  units = "in",
  dpi = 600
)


# Close log
sink()

