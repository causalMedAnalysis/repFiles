# Preliminaries
chapter <- "ch3"
title <- "figure_3-10"
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"
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

# Open log
sink(log_path, split = TRUE)

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/figure_3-10.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/JOBSII/Jobs-NoMiss-Binary.dta

# Outputs:     .../code/ch3/_LOGS/figure_3-10_log.txt

# Description: Replicates Chapter 3, Figure 3-10: Estimates of Natural Direct 
#              and Indirect Effects from JOBSII Adjusted for Bias due to 
#              Unobserved Mediator-Outcome Confounding.
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#
packages <-
  c(
    "tidyverse",
    "haven",
    "gridExtra",
    "metR", 
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
Y <- "work1"

# exposure
D <- "treat"

# mediator
M <- "job_seek"

# baseline confounders
C <- c(
  "econ_hard",
  "sex",
  "age",
  "nonwhite",
  "educ",
  "income"
)

# mediator value for CDE
m <- 4

#----------------#
#  PREPARE DATA  #
#----------------#
jobs_raw <- read_stata(
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/JOBSII/Jobs-NoMiss-Binary.dta"
)

jobs <- jobs_raw |>
  zap_labels()

#--------------------#
#  ESTIMATE EFFECTS  #
#--------------------#
# Linear model with D x M interaction

out <- linmed(
  data = jobs,
  D = D,
  M = M,
  Y = Y,
  C = C,
  m = m,
  interaction_DM = TRUE
)

#---------------------------#
#  BIAS-ADJUSTED ESTIMATES  #
#---------------------------#
# Specify range for sensitivity parameters under M-Y confounding
sens_grid <- expand.grid(
  delta_UYgivCDM = seq(-0.1, 0.1, 0.01),
  delta_DUgivCM = seq(-0.1, 0.1, 0.01)
)

# Compute bias-adjusted estimates
adj_grid <- cbind(
  sens_grid,
  nde_adj = out$NDE - (sens_grid$delta_UYgivCDM * sens_grid$delta_DUgivCM),
  nie_adj = out$NIE + (sens_grid$delta_UYgivCDM * sens_grid$delta_DUgivCM)
)

#------------------------#
#  CREATE CONTOUR PLOTS  #
#------------------------#
# Create contour plots of bias-adjusted estimates against sensitivity parameters
# Common specifications
common_specs <- list(
  scale_colour_distiller(palette="Greys", direction=1),
  xlab(expression(delta["UY|C,D,M"])),
  ylab(expression(delta["DU|C,M"])),
  scale_x_continuous(breaks=seq(-0.1, 0.1, 0.02)),
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.02)),
  theme_bw(base_size=11),
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
)

# NDE plot
plot_nde <- ggplot(
  adj_grid,
  aes(x=delta_UYgivCDM, y=delta_DUgivCM, z=nde_adj, color=after_stat(level))
) +  
  geom_contour(
    breaks = seq(
      round(min(adj_grid$nde_adj), 3),
      round(max(adj_grid$nde_adj), 3),
      0.002
    ),
    show.legend = FALSE
  ) +
  ggtitle("A. Bias-adjusted NDE(1,0) Estimates") +
  geom_text_contour(
    breaks = seq(
      round(min(adj_grid$nde_adj), 3),
      round(max(adj_grid$nde_adj), 3),
      0.002
    ),
    stroke = 0.3,
    size = 3,
    skip = 0,
    color = "black"
  ) +
  common_specs

# NIE plot
plot_nie <- ggplot(
  adj_grid,
  aes(x=delta_UYgivCDM, y=delta_DUgivCM, z=nie_adj, color=after_stat(level))
) +  
  geom_contour(
    breaks = seq(
      round(min(adj_grid$nie_adj), 3),
      round(max(adj_grid$nie_adj), 3),
      0.002
    ),
    show.legend = FALSE
  ) +
  ggtitle("B. Bias-adjusted NIE(1,0) Estimates") +
  geom_text_contour(
    breaks = seq(
      round(min(adj_grid$nie_adj), 3),
      round(max(adj_grid$nie_adj), 3),
      0.002
    ),
    stroke = 0.3,
    size = 3,
    skip = 0,
    color = "black"
  ) +
  common_specs

# Combine the plots
plot_comb <- grid.arrange(plot_nde, plot_nie, ncol=1)

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
