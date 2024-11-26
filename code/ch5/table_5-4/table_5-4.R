# Preliminaries
chapter <- "ch5"
title <- "table_5-4"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch5/table_5-4.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch5/_LOGS/table_5-4_log.txt

# Description: Replicates Chapter 5, Table 5.4: Case Counts and Sample Means for 
#              CES-D Scores (Y) by College Attendance (D), Unemployment Status 
#              (M1), Household Income (M2), and Maternal Education (C), NLSY. 
#              Also replicates nonparametric estimates of path-specific effects, 
#              which are reported in the text preceding Table 5.4.
#-------------------------------------------------------------------------------


#-------------#
#  LIBRARIES  #
#-------------#
library(tidyverse)
library(haven)




#------------------#
#  SPECIFICATIONS  #
#------------------#
# outcome
Y <- "std_cesd_age40"

# exposure
D <- "att22"

# mediators
M <- c(
  "ever_unemp_age3539",
  "incgt50k"
)

# baseline confounder(s)
C <- "momcol"

# key variables
key_vars <- c(
  "cesd_age40", # unstandardized version of Y
  D,
  M[1],
  "faminc_adj_age3539", # source variable for M[2]
  "momedu" # source variable for C
)




#----------------#
#  PREPARE DATA  #
#----------------#
nlsy_raw <- read_stata(
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
)

nlsy <- nlsy_raw[complete.cases(nlsy_raw[,key_vars]),] |>
  mutate(
    momcol = as.numeric(momedu>12),
    incgt50k = as.numeric(faminc_adj_age3539>=50000),
    std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40)
  )




#-------------------------------#
#  CASE COUNTS & OUTCOME MEANS  #
#-------------------------------#
agg <- nlsy |>
  group_by(att22, ever_unemp_age3539, incgt50k, momcol) |>
  summarize(.groups = "drop",
    mean = mean(std_cesd_age40),
    n = n()
  )
agg |>
  mutate(
    mean = round(mean,2)
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
# The following estimates are reported in the text before Table 5.4.

# First, add counts and proportions aggregated to various levels
agg2 <- agg |>
  group_by(pick(all_of(C))) |>
  mutate(n_C = sum(n)) |>
  group_by(pick(all_of(c(C,D,M[1])))) |>
  mutate(n_CDM1 = sum(n)) |>
  group_by(pick(all_of(c(C,D)))) |>
  mutate(n_CD = sum(n)) |>
  ungroup() |>
  mutate(
    n_tot = sum(n),
    prop_C = n_C / n_tot,
    prop_M1_CD = n_CDM1 / n_CD,
    prop_M2_CDM1 = n / n_CDM1,
    sign = ifelse(as.logical(.data[[D]]), 1, -1)
  )


# Estimate ATE(1,0)
ATEhat_npl <- agg2 |>
  summarize(
    est = sum(sign * mean * prop_M2_CDM1 * prop_M1_CD * prop_C)
  ) |>
  pull(est)


# Estimate PSE_{D->Y}(1,0)
PSE1hat_npl <- agg2 |>
  left_join(
    agg2 |>
      filter(.data[[D]]==0) |>
      select(all_of(c(C,M)), prop_M1_CD, prop_M2_CDM1) |>
      rename(
        prop_M1_Cdstar = prop_M1_CD,
        prop_M2_CdstarM1 = prop_M2_CDM1
      )
    ,
    by = c(C,M)
  ) |>
  summarize(
    est = sum(sign * mean * prop_M2_CdstarM1 * prop_M1_Cdstar * prop_C)
  ) |>
  pull(est)


# Estimate PSE_{D->M2->Y}(1,0)
PSE2hat_npl <- agg2 |>
  left_join(
    agg2 |>
      filter(.data[[D]]==0) |>
      select(all_of(c(C,M)), prop_M1_CD, prop_M2_CDM1) |>
      rename(
        prop_M1_Cdstar = prop_M1_CD,
        prop_M2_CdstarM1 = prop_M2_CDM1
      )
    ,
    by = c(C,M)
  ) |>
  filter(.data[[D]]==1) |>
  summarize(
    est = sum((prop_M2_CDM1 - prop_M2_CdstarM1) * mean * prop_M1_Cdstar * prop_C)
    # ^ since we have filtered to D==1, prop_M2_CDM1 is equal to prop_M2_CdM1 
    # (i.e., the estimate of P(M2|C,d,M1))
  ) |>
  pull(est)


# Estimate PSE_{D->M1~>Y}(1,0)
PSE3hat_npl <- agg2 |>
  left_join(
    agg2 |>
      filter(.data[[D]]==0) |>
      select(all_of(c(C,M)), prop_M1_CD) |>
      rename(prop_M1_Cdstar = prop_M1_CD)
    ,
    by = c(C,M)
  ) |>
  filter(.data[[D]]==1) |>
  summarize(
    est = sum((prop_M1_CD - prop_M1_Cdstar) * mean * prop_M2_CDM1 * prop_C)
    # ^ since we have filtered to D==1, prop_M1_CD is equal to prop_M1_Cd 
    # (i.e., the estimate of P(M1|C,d))
  ) |>
  pull(est)


# Collate estimates
npest <- data.frame(
  param = c("ATE(1,0)", "PSE_{D->Y}(1,0)", "PSE_{D->M2->Y}(1,0)", "PSE_{D->M1~>Y}(1,0)"),
  est = c(ATEhat_npl, PSE1hat_npl, PSE2hat_npl, PSE3hat_npl)
)

npest |>
  mutate(
    est = round(est, 2)
  )


# Close log
sink()

