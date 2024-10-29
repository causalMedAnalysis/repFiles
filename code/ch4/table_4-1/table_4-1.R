# Preliminaries
chapter <- "ch4"
title <- "table_4-1"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch4/table_4-1.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch4/_LOGS/table_4-1_log.txt

# Description: Replicates Chapter 4, Table 4.1: Case Counts and Sample Means for 
#              CES-D Scores (Y) by College Attendance (D), Unemployment Status 
#              (L), Household Income (M), and Maternal Education (C), NLSY.
#              Also replicates nonparametric estimates of the CDE, IDE, IIE, and 
#              OE, which are reported in the text following Table 4-1.
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

# mediator
M <- "incgt50k"

# exposure-induced confounder
L <- "ever_unemp_age3539"

# baseline confounder(s)
C <- "momcol"

# key variables
key_vars <- c(
  "cesd_age40", # unstandardized version of Y
  D,
  "faminc_adj_age3539", # source variable for C
  L,
  "momedu" # source variable for C
)

# mediator value for CDE
m <- 1




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
# The following estimates are reported in the text after Table 4-1.

# First, add counts and proportions aggregated to various levels
agg2 <- agg |>
  group_by(pick(all_of(C))) |>
  mutate(n_C = sum(n)) |>
  group_by(pick(all_of(c(L,C,D)))) |>
  mutate(n_LCD = sum(n)) |>
  group_by(pick(all_of(c(M,C,D)))) |>
  mutate(n_MCD = sum(n)) |>
  group_by(pick(all_of(c(C,D)))) |>
  mutate(n_CD = sum(n)) |>
  ungroup() |>
  mutate(
    n_tot = sum(n),
    prop_C = n_C / n_tot,
    prop_L_CD = n_LCD / n_CD,
    prop_M_CD = n_MCD / n_CD,
    sign = ifelse(as.logical(.data[[D]]), 1, -1)
  )


# Estimate CDE(1,0,50K+)
CDEhat_npl <- agg2 |>
  filter(.data[[M]]==m) |>
  mutate(
    product = sign * mean * prop_L_CD * prop_C
  ) |>
  summarize(
    est = sum(product)
  ) |>
  pull(est)


# Estimate IDE(1,0)
IDEhat_npl <- agg2 |>
  left_join(
    agg2 |>
      filter(.data[[D]]==0) |>
      select(all_of(c(M,L,C)), prop_M_CD) |>
      rename(prop_M_Cdstar = prop_M_CD)
    ,
    by = c(M,L,C)
  ) |>
  mutate(
    product = sign * mean * prop_L_CD * prop_M_Cdstar * prop_C
  ) |>
  summarize(
    est = sum(product)
  ) |>
  pull(est)


# Estimate IIE(1,0)
IIEhat_npl <- agg2 |>
  left_join(
    agg2 |>
      filter(.data[[D]]==0) |>
      select(all_of(c(M,L,C)), prop_M_CD) |>
      rename(prop_M_Cdstar = prop_M_CD)
    ,
    by = c(M,L,C)
  ) |>
  filter(.data[[D]]==1) |>
  mutate(
    product = (prop_M_CD - prop_M_Cdstar) * mean * prop_L_CD * prop_C
    # ^ since we have filtered to D==1, prop_M_CD is equal to prop_M_Cd 
    # (i.e., the estimate of P(M|C,d))
  ) |>
  summarize(
    est = sum(product)
  ) |>
  pull(est)


# Estimate OE(1,0)
OEhat_npl <- IDEhat_npl + IIEhat_npl


# Collate estimates
npest <- data.frame(
  param = c("CDE(1,0,50K+)", "IDE(1,0)", "IIE(1,0)", "OE(1,0)"),
  est = c(CDEhat_npl, IDEhat_npl, IIEhat_npl, OEhat_npl)
)

npest |>
  mutate(
    est = round(est, 2)
  )


# Close log
sink()

