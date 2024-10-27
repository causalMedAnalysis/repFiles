# Preliminaries
chapter <- "ch3"
title <- "table_3-1"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/table_3-1.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/impcde.R

# Outputs:     .../code/ch3/_LOGS/table_3-1_log.txt

# Description: Replicates Chapter 3, Table 3-1: Case Counts and Sample Means for 
#              CES-D Scores (Y) by College Attendance (D), Unemployment Status 
#              (M), and Maternal Education (C), NLSY79.
#              Also replicates nonparametric estimates of the ATE, CDE, NDE, and 
#              NIE, which are reported in the text following Table 3-1.
#-------------------------------------------------------------------------------


#-------------#
#  LIBRARIES  #
#-------------#
library(margins)
library(tidyverse)
library(haven)




#-----------------------------#
#  LOAD CAUSAL MED FUNCTIONS  #
#-----------------------------#
# regression imputation CDE estimator
#source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/impcde.R")
source("C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/test project/R/impcde.R")




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
C <- "momcol"

# key variables
key_vars <- c(
  "cesd_age40", # unstandardized version of Y
  D,
  M,
  "momedu" # source variable for C
)

# mediator value for CDE
m <- 0




#----------------#
#  PREPARE DATA  #
#----------------#
nlsy_raw <- read_stata(
  #file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
  file = "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Data/NLSY79/nlsy79BK_ed2.dta"
)

nlsy <- nlsy_raw[complete.cases(nlsy_raw[,key_vars]),] |>
  mutate(
    momcol = as.numeric(momedu>12),
    std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
  )




#-------------------------------#
#  CASE COUNTS & OUTCOME MEANS  #
#-------------------------------#
nlsy |>
  # duplicate the dataframe to create a sub-total across the mediator
  mutate(
    ever_unemp_age3539 = 9
  ) |>
  bind_rows(nlsy) |>
  mutate(
    ever_unemp_age3539_fct = factor(
      ever_unemp_age3539,
      levels = c(0,1,9),
      labels = c("0","1","Total")
    )
  ) |>
  # get case counts and outcome means, by D, M, and C
  group_by(att22, ever_unemp_age3539_fct, momcol) |>
  summarize(.groups = "drop",
    mean = round(mean(std_cesd_age40),2),
    n = n()
  ) |>
  pivot_wider(
    names_from = momcol,
    names_prefix = "C",
    names_vary = "slowest",
    values_from = c(mean, n)
  )




#---------------------------#
#  NONPARAMETRIC ESTIMATES  #
#---------------------------#
# The following estimates are reported in the text after Table 3-1.

# Estimate ATE(1,0)
m1 <- lm(
  std_cesd_age40 ~ factor(att22)*momcol,
  data = nlsy
)

ATEhat <- margins_summary(m1) |>
  filter(factor==paste0(D,"1")) |>
  pull(AME)


# Estimate CDE(1,0,0)
m2 <- lm(
  std_cesd_age40 ~ att22*ever_unemp_age3539*momcol,
  data = nlsy
)

CDE0hat <- impcde(
  data = nlsy,
  model_y = m2,
  D = D,
  M = M,
  m = m
)


# Estimate NDE(1,0) and NIE(1,0)
gdata <- nlsy

gdata[[D]] <- 0
gdata[[M]] <- 0
EhatY_D0M0C <- predict(m2, gdata)

gdata[[M]] <- 1
EhatY_D0M1C <- predict(m2, gdata)

gdata[[D]] <- 1
gdata[[M]] <- 0
EhatY_D1M0C <- predict(m2, gdata)

gdata[[M]] <- 1
EhatY_D1M1C <- predict(m2, gdata)

m3 <- lm(
  ever_unemp_age3539 ~ att22*momcol,
  data = nlsy
)
gdata <- nlsy

gdata[[D]] <- 0
PhatM_D0C <- predict(m3, gdata)

gdata[[D]] <- 1
PhatM_D1C <- predict(m3, gdata)

NDEhat <- mean(
  (EhatY_D1M0C-EhatY_D0M0C) * (1-PhatM_D0C) + 
  (EhatY_D1M1C-EhatY_D0M1C) * PhatM_D0C
)
NIEhat <- mean(
  (EhatY_D1M0C) * ((1-PhatM_D1C) - (1-PhatM_D0C)) +
  (EhatY_D1M1C) * (PhatM_D1C - PhatM_D0C)
)


# Collate estimates
npest <- data.frame(
  param = c("ATE(1,0)", "CDE(1,0,0)", "NDE(1,0)", "NIE(1,0)"),
  est = c(ATEhat, CDE0hat, NDEhat, NIEhat)
)

npest |>
  mutate(
    est = round(est, 2)
  )


# Close log
sink()

