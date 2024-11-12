# Preliminaries
chapter <- "ch4"
title <- "table_4-4"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch4/table_4-4.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch4/_LOGS/table_4-4_log.txt

# Description: Replicates Chapter 4, Table 4.4: Interventional Effects of 
#              College Attendance on CES-D Scores as Estimated from the NLSY 
#              Using Inverse Probability Weighting.
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
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")




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




#------------------------------#
#  DEFINE CUSTOM IPW FUNCTION  #
#------------------------------#
custom_ipwvent <- function(
  data,
  D,
  M,
  Y,
  L,
  m = 0,
  D_formula,
  L_formula,
  M_formula,
  stabilize = TRUE,
  censor = TRUE,
  censor_low = 0.01,
  censor_high = 0.99
) {
  # fit specified models
  ## estimate f(D|C)
  D_model <- glm(
    D_formula,
    data = data,
    family = binomial(link = "logit")
  )
  ## estimate q(L|C,D)
  L_model <- glm(
    L_formula,
    data = data,
    family = binomial(link = "logit")
  )
  ## estimate g(M|C,D,L)
  M_model <- lm(
    M_formula,
    data = data
  )
  
  # D and L values
  D_values <- unique(data[[D]])
  L_values <- unique(data[[L]])
  cross_DL <- expand.grid(D_val = D_values, L_val = L_values)
  
  # estimate p(D|C) for each value of D
  pDd_C <- lapply(
    D_values,
    \(d) dbinom(
      d,
      size = 1,
      prob = predict(D_model, newdata = data, type = "response")
    )
  )
  names(pDd_C) <- paste0("pD", D_values, "_C")
  list2env(pDd_C, envir = environment())
  
  # estimate p(D|C) for the observed D
  pD_C <- dbinom(
    data[[D]],
    size = 1,
    prob = predict(D_model, newdata = data, type = "response")
  )
  
  # estimate p(D) for each value of D
  pD1 <- mean(data[[D]])
  pD0 <- 1 - pD1
  
  # estimate p(D) for the observed D
  pD <- ifelse(as.logical(data[[D]]), pD1, pD0)
  
  # estimate p(L|D,C) for each value of L and D
  pLl_DdC <- mapply(
    function(d, l) {
      temp_data <- data
      temp_data[[D]] <- d
      dbinom(
        l,
        size = 1,
        prob = predict(L_model, newdata = temp_data, type = "response")
      )
    },
    cross_DL$D_val,
    cross_DL$L_val,
    SIMPLIFY = FALSE
  )
  names(pLl_DdC) <- paste0("pL", cross_DL$L_val, "_D", cross_DL$D_val, "C")
  list2env(pLl_DdC, envir = environment())
  
  # estimate p(M|D,L,C) for the observed M and each value of D and L
  pM_DdLlC <- mapply(
    function(d, l) {
      temp_data <- data
      temp_data[[D]] <- d
      temp_data[[L]] <- l
      dnorm(
        temp_data[[M]],
        mean = predict(M_model, newdata = temp_data, type = "response"),
        sd = sigma(M_model)
      )
    },
    cross_DL$D_val,
    cross_DL$L_val,
    SIMPLIFY = FALSE
  )
  names(pM_DdLlC) <- paste0("pM_D", cross_DL$D_val, "L", cross_DL$L_val, "C")
  list2env(pM_DdLlC, envir = environment())
  
  # estimate p(M|D,L,C) for the observed M, observed L, and each value of D
  pM_DdLC <- lapply(
    D_values,
    function(d) {
      temp_data <- data
      temp_data[[D]] <- d
      dnorm(
        temp_data[[M]],
        mean = predict(M_model, newdata = temp_data, type = "response"),
        sd = sigma(M_model)
      )
    }
  )
  names(pM_DdLC) <- paste0("pM_D", D_values, "LC")
  list2env(pM_DdLC, envir = environment())
  
  # estimate p(M|D,L,C) for the observed M, observed D, and observed L
  pM_DLC <- dnorm(
    data[[M]],
    mean = predict(M_model, newdata = data, type = "response"),
    sd = sigma(M_model)
  )
  
  # estimate p(M|D) for the observed M and observed D
  M_mean <- ave(data[[M]], data[[D]])
  M_dev <- data[[M]] - M_mean
  pM_D <- dnorm(
    data[[M]],
    mean = M_mean,
    sd = sd(M_dev)
  )
  
  # for convenience, create logical vectors to identify subgroups
  group_D0 <- data[[D]]==0
  group_D1 <- data[[D]]==1
  
  # create IPWs
  w1 <- ifelse(
    group_D0,
    ((pM_D0L0C * pL0_D0C) + (pM_D0L1C * pL1_D0C)) /
      (pD0_C * pM_D0LC)
    ,
    0
  )
  w2 <- ifelse(
    group_D1,
    ((pM_D1L0C * pL0_D1C) + (pM_D1L1C * pL1_D1C)) /
      (pD1_C * pM_D1LC)
    ,
    0
  )
  w3 <- ifelse(
    group_D1,
    ((pM_D0L0C * pL0_D0C) + (pM_D0L1C * pL1_D0C)) /
      (pD1_C * pM_D1LC)
    ,
    0
  )
  w4 <- 1 / (pM_DLC * pD_C)
  
  # stabilize IPWs
  if (stabilize) {
    w1 <- w1 * pD0
    w2 <- w2 * pD1
    w3 <- w3 * pD1
    w4 <- w4 * pM_D * pD
  }
  
  # censor IPWs (among the appropriate subgroups with non-zero IPWs)
  if (censor) {
    w1[group_D0] <- trimQ(w1[group_D0], low = censor_low, high = censor_high)
    w2[group_D1] <- trimQ(w2[group_D1], low = censor_low, high = censor_high)
    w3[group_D1] <- trimQ(w3[group_D1], low = censor_low, high = censor_high)
    w4 <- trimQ(w4, low = censor_low, high = censor_high)
  }
  
  # estimate OE, IDE, and IIE
  Ehat_Y0M0 <- weighted.mean(data[[Y]], w1)
  Ehat_Y1M1 <- weighted.mean(data[[Y]], w2)
  Ehat_Y1M0 <- weighted.mean(data[[Y]], w3)
  OE  <- Ehat_Y1M1 - Ehat_Y0M0
  IDE <- Ehat_Y1M0 - Ehat_Y0M0
  IIE <- Ehat_Y1M1 - Ehat_Y1M0
  
  # estimate CDE
  Y_model <- lm(
    as.formula(paste0(Y,"~",D,"*",M)),
    data = data,
    weights = w4
  )
  CDE <- 
    Y_model$coefficients[[D]] + 
    Y_model$coefficients[[paste0(D,":",M)]] * m
  
  # compile and output
  out <- list(
    OE = OE,
    IDE = IDE,
    IIE = IIE,
    CDE = CDE,
    weights1 = w1,
    weights2 = w2,
    weights3 = w3,
    weights4 = w4,
    model_D = D_model,
    model_L = L_model,
    model_M = M_model
  )
  return(out)
}




#-------------#
#  VERSION 1  #
#-------------#
# D model: Additive logit model
# L model: Additive logit model
# M model: Additive linear model

# D model formula
predictors1_D <- paste(C, collapse = " + ")
(formula1_D_string <- paste(D, "~", predictors1_D))
formula1_D <- as.formula(formula1_D_string)

# L model formula
predictors1_L <- paste(c(D,C), collapse = " + ")
(formula1_L_string <- paste(L, "~", predictors1_L))
formula1_L <- as.formula(formula1_L_string)

# M model formula
predictors1_M <- paste(c(D,C,L), collapse = " + ")
(formula1_M_string <- paste(M, "~", predictors1_M))
formula1_M <- as.formula(formula1_M_string)

# Estimate effects
out1 <- custom_ipwvent(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  L = L,
  m = m,
  D_formula = formula1_D,
  L_formula = formula1_L,
  M_formula = formula1_M
)




#-------------#
#  VERSION 2  #
#-------------#
# D model: Additive logit model
# L model: Logit model with D x C interactions
# M model: Linear model with D x C, D x L interactions

# D model formula
formula2_D <- formula1_D

# L model formula
## main effects
predictors2_L <- paste(c(D,C), collapse = " + ")
## D x C interactions
predictors2_L <- paste(
  predictors2_L,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## full formula
(formula2_L_string <- paste(L, "~", predictors2_L))
formula2_L <- as.formula(formula2_L_string)

# M model formula
## main effects
predictors2_M <- paste(c(D,C,L), collapse = " + ")
## D x C interactions
predictors2_M <- paste(
  predictors2_M,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## D x L interaction
predictors2_M <- paste(
  predictors2_M,
  "+",
  paste(D, L, sep = ":", collapse = " + ")
)
## full formula
(formula2_M_string <- paste(M, "~", predictors2_M))
formula2_M <- as.formula(formula2_M_string)

# Estimate effects
out2 <- custom_ipwvent(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  L = L,
  m = m,
  D_formula = formula2_D,
  L_formula = formula2_L,
  M_formula = formula2_M
)




#---------------------#
#  COLLATE ESTIMATES  #
#---------------------#
master <- data.frame(
  param = c("OE(1,0)", "IDE(1,0)", "IIE(1,0)", "CDE(1,0,ln(50K))"),
  est_v1 = c(
    out1$OE,
    out1$IDE,
    out1$IIE,
    out1$CDE
  ),
  est_v2 = c(
    out2$OE,
    out2$IDE,
    out2$IIE,
    out2$CDE
  )
)

master |>
  mutate(
    across(
      .cols = starts_with("est_"),
      .fns = \(x) round(x, 3)
    )
  )


# Close log
sink()

