# Preliminaries
chapter <- "ch4"
title <- "figure_4-9"
dir_root <- "C:/Users/ashiv/OneDrive/Documents/Wodtke/Causal Mediation Analysis Book/Programming/Programs/Replication"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")
dir_fig <- paste0(dir_root, "/figures/", chapter)

# Open log
sink(log_path, split = TRUE)
#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch4/table_4-9.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/plowUse/plowUse.dta
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/medsim.R
#              https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/ipwvent.R

# Outputs:     .../code/ch4/_LOGS/figure_4-9_log.txt
#              .../figures/ch4/figure_4-9.png

# Description: Replicates Chapter 4, Figure 4.9: Interventional and Controlled 
#              Direct Effects of Historical Plow Use on Women's Contemporary 
#              Representation in Government.
#-------------------------------------------------------------------------------


#------------------------#
#  INSTALL DEPENDENCIES  #
#------------------------#
# First, install dependencies available on CRAN.
dependencies_cran <- c(
  # The following package is used to fit an ordered logistic regression model.
  "VGAM",
  # The following three packages are used to parallelize the bootstrap.
  "doParallel",
  "doRNG",
  "foreach"
)

#install.packages(dependencies_cran)
# ^ Uncomment this line above to install these packages.

# (And note that, once you have installed these packages, there is no need for 
# you to load these packages with the library function to run the code in this 
# script.)


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
library(VGAM)
library(tidyverse)
library(haven)




#-----------------------------#
#  LOAD CAUSAL MED FUNCTIONS  #
#-----------------------------#
# utilities
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")
# RWR estimator
source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/rwrlite.R")




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

# number of simulations for simulation estimator
n_sims <- 2000

# number of bootstrap replications
n_reps <- 2000




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
# D model: Additive
# L model: Additive
# M model (for RWR and simulation): Additive
# M model (for IPW): Additive, with L term included
# Y model: With D x M interaction

# D model formula
predictors_D <- paste(C, collapse = " + ")
(formula_D_string <- paste(D, "~", predictors_D))
formula_D <- as.formula(formula_D_string)

# L and M model formulae
predictors_LM <- paste(c(D,C), collapse = " + ")
(formula_L_string <- paste(L, "~", predictors_LM))
(formula_M_string <- paste(M, "~", predictors_LM))
formula_L <- as.formula(formula_L_string)
formula_M <- as.formula(formula_M_string)

# M model formula for IPW
predictors2_M <- paste(c(D,C,L), collapse = " + ")
(formula2_M_string <- paste(M, "~", predictors2_M))
formula2_M <- as.formula(formula2_M_string)

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




#-----------------#
#  RWR ESTIMATOR  #
#-----------------#
out_rwr <- rwrlite(
  data = plow,
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
# We will fit the following models:
## L model: ordinal logit
## M model: linear
## Y model: log-binomial

# However, the medsim function does not currently support log-binomial models, 
# so we will create a custom function for this script.

# Define inner custom simulation function
# ----------------------------------------
custom_sim_inner <- function(
  data,
  D,
  M,
  Y,
  L,
  m = 0,
  L_formula,
  M_formula,
  Y_formula,
  n_sims = 1000,
  minimal = FALSE
) {
  # fit specified models
  ## estimate q(L|C,D)
  L_model <- vgam(
    L_formula,
    family = cumulative(parallel = TRUE, reverse = FALSE),
    data = data
  )
  ## estimate g(M|C,D)
  M_model <- lm(
    M_formula,
    data = data
  )
  ## estimate h(Y|C,D,L,M)
  Y_model <- glm(
    Y_formula,
    data = data,
    family = quasibinomial(link = "logit")
  )
  
  # create copies of the data
  sim_data <- sim_data0 <- sim_data1 <- data
  sim_data0[[D]] <- 0
  sim_data1[[D]] <- 1
  
  # predict p from q-hat(L|C,D)
  ## D=0
  phat_L0 <- predict(L_model, newdata = sim_data0, type = "response")
  phat1_L0 <- phat_L0[,1]
  phat2_L0 <- phat_L0[,2]
  phat3_L0 <- phat_L0[,3]
  ## D=1
  phat_L1 <- predict(L_model, newdata = sim_data1, type = "response")
  phat1_L1 <- phat_L1[,1]
  phat2_L1 <- phat_L1[,2]
  phat3_L1 <- phat_L1[,3]
  
  # predict M from g-hat(M|C,D)
  ## D=0
  ehat_M0 <- predict(M_model, newdata = sim_data0, type = "response")
  ## D=1
  ehat_M1 <- predict(M_model, newdata = sim_data1, type = "response")
  
  # initialize potential outcome vectors
  Y0L0M0 <- Y1L1M1 <- Y1L1M0 <- Y0L0m <- Y1L1m <- rep(NA_real_, n_sims)
  
  
  # loop over simulations
  for (sim in seq_len(n_sims)) {
    ## simulate L values for D=0
    runif <- runif(nrow(sim_data))
    L0 <- ifelse(
      runif <= phat1_L0,
      1,
      ifelse(
        runif > phat1_L0 & runif <= (phat1_L0 + phat2_L0),
        2,
        ifelse(
          runif > (phat1_L0 + phat2_L0),
          3,
          NA_real_
        )
      )
    )
    
    ## simulate M values for D=0
    M0 <- rnorm(nrow(sim_data), mean = ehat_M0, sd = sigma(M_model))
    
    ## simulate L values for D=1
    runif <- runif(nrow(sim_data))
    L1 <- ifelse(
      runif <= phat1_L1,
      1,
      ifelse(
        runif > phat1_L1 & runif <= (phat1_L1 + phat2_L1),
        2,
        ifelse(
          runif > (phat1_L1 + phat2_L1),
          3,
          NA_real_
        )
      )
    )
    
    ## simulate M values for D=1
    M1 <- rnorm(nrow(sim_data), mean = ehat_M1, sd = sigma(M_model))
    
    ## simulate potential outcomes Y(1, M-script(1|C))
    temp_sim_data <- sim_data
    temp_sim_data[[D]] <- 1
    temp_sim_data[[L]] <- L1
    temp_sim_data[[M]] <- M1
    phat_Y1L1M1 <- predict(Y_model, newdata = temp_sim_data, type = "response")
    Y1L1M1[sim] <- mean(
      rbinom(n = nrow(temp_sim_data), size = 100, prob = phat_Y1L1M1)
      /
      100
    )
    
    ## simulate potential outcomes Y(1, M-script(0|C))
    temp_sim_data <- sim_data
    temp_sim_data[[D]] <- 1
    temp_sim_data[[L]] <- L1
    temp_sim_data[[M]] <- M0
    phat_Y1L1M0 <- predict(Y_model, newdata = temp_sim_data, type = "response")
    Y1L1M0[sim] <- mean(
      rbinom(n = nrow(temp_sim_data), size = 100, prob = phat_Y1L1M0)
      /
      100
    )
    
    ## simulate potential outcomes Y(0, M-script(0|C))
    temp_sim_data <- sim_data
    temp_sim_data[[D]] <- 0
    temp_sim_data[[L]] <- L0
    temp_sim_data[[M]] <- M0
    phat_Y0L0M0 <- predict(Y_model, newdata = temp_sim_data, type = "response")
    Y0L0M0[sim] <- mean(
      rbinom(n = nrow(temp_sim_data), size = 100, prob = phat_Y0L0M0)
      /
      100
    )
    
    ## simulate potential outcomes Y(0, m)
    temp_sim_data <- sim_data
    temp_sim_data[[D]] <- 0
    temp_sim_data[[L]] <- L0
    temp_sim_data[[M]] <- m
    phat_Y0L0m <- predict(Y_model, newdata = temp_sim_data, type = "response")
    Y0L0m[sim] <- mean(
      rbinom(n = nrow(temp_sim_data), size = 100, prob = phat_Y0L0m)
      /
      100
    )
    
    ## simulate potential outcomes Y(1, m)
    temp_sim_data <- sim_data
    temp_sim_data[[D]] <- 1
    temp_sim_data[[L]] <- L1
    temp_sim_data[[M]] <- m
    phat_Y1L1m <- predict(Y_model, newdata = temp_sim_data, type = "response")
    Y1L1m[sim] <- mean(
      rbinom(n = nrow(temp_sim_data), size = 100, prob = phat_Y1L1m)
      /
      100
    )
  }
  
  
  # estimate effects
  IDE <- mean(Y1L1M0) - mean(Y0L0M0)  
  IIE <- mean(Y1L1M1) - mean(Y1L1M0)
  OE <- mean(Y1L1M1) - mean(Y0L0M0)
  CDE <- mean(Y1L1m) - mean(Y0L0m)
  
  # compile and output
  if (minimal) {
    out <- list(
      OE = OE,
      IDE = IDE,
      IIE = IIE,
      CDE = CDE
    )
  }
  else {
    out <- list(
      OE = OE,
      IDE = IDE,
      IIE = IIE,
      CDE = CDE,
      model_L = L_model,
      model_M = M_model,
      model_Y = Y_model
    )
  }
  return(out)
}


# Define outer custom simulation function (bootstrapping the inner function)
# ----------------------------------------
custom_sim <- function(
  data,
  D,
  M,
  Y,
  L,
  m = 0,
  L_formula,
  M_formula,
  Y_formula,
  n_sims = 1000,
  boot = FALSE,
  boot_reps = 1000,
  boot_conf_level = 0.95,
  boot_seed = NULL,
  boot_parallel = FALSE,
  boot_cores = max(c(parallel::detectCores()-2,1))
) {
  # load data
  data_outer <- data
  
  
  # create adjusted boot_parallel logical
  boot_parallel_rev <- ifelse(boot_cores>1, boot_parallel, FALSE)
  
  
  # preliminary error/warning checks for the bootstrap
  if (boot) {
    if (boot_parallel & boot_cores==1) {
      warning(paste(strwrap("Warning: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but you do not have enough cores available for parallelization. The bootstrap will proceed without parallelization."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("doParallel", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'doParallel' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("doRNG", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'doRNG' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("foreach", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'foreach' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
  }
  
  
  # compute point estimates
  est <- custom_sim_inner(
    data = data_outer,
    D = D,
    M = M,
    Y = Y,
    L = L,
    m = m,
    L_formula = L_formula,
    M_formula = M_formula,
    Y_formula = Y_formula,
    n_sims = n_sims,
    minimal = FALSE
  )
  
  
  # bootstrap, if requested
  if (boot) {
    # bootstrap function
    boot_fnc <- function() {
      # sample from the data with replacement
      boot_data <- data_outer[sample(nrow(data_outer), size = nrow(data_outer), replace = TRUE), ]
      
      # compute point estimates in the replicate sample
      custom_sim_inner(
        data = boot_data,
        D = D,
        M = M,
        Y = Y,
        L = L,
        m = m,
        L_formula = L_formula,
        M_formula = M_formula,
        Y_formula = Y_formula,
        n_sims = n_sims,
        minimal = TRUE
      )
    }
    
    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type="PSOCK")
      doParallel::registerDoParallel(cl=x_cluster)
      parallel::clusterExport(
        cl = x_cluster, 
        varlist = c("custom_sim_inner", "vgam", "cumulative"),
        envir = environment()
      )
      `%dopar%` <- foreach::`%dopar%`
    }
    
    # set seed
    if (!is.null(boot_seed)) {
      set.seed(boot_seed)
      if (boot_parallel) {
        doRNG::registerDoRNG(boot_seed)
      }
    }
    
    # compute estimates for each replicate sample
    if (boot_parallel_rev) {
      boot_res <- foreach::foreach(i = 1:boot_reps, .combine = comb_list_vec) %dopar% {
        boot_fnc()
      }
      boot_OE <- boot_res$OE
      boot_IDE <- boot_res$IDE
      boot_IIE <- boot_res$IIE
      boot_CDE <- boot_res$CDE
    }
    else {
      boot_OE <- rep(NA_real_, boot_reps)
      boot_IDE <- rep(NA_real_, boot_reps)
      boot_IIE <- rep(NA_real_, boot_reps)
      boot_CDE <- rep(NA_real_, boot_reps)
      for (i in seq_len(boot_reps)) {
        boot_iter <- boot_fnc()
        boot_OE[i] <- boot_iter$OE
        boot_IDE[i] <- boot_iter$IDE
        boot_IIE[i] <- boot_iter$IIE
        boot_CDE[i] <- boot_iter$CDE
      }
    }
    
    # clean up
    if (boot_parallel_rev) {
      parallel::stopCluster(x_cluster)
      rm(x_cluster)
    }
    
    # compute bootstrap confidence intervals 
    # from percentiles of the bootstrap distributions
    boot_alpha <- 1 - boot_conf_level
    boot_ci_probs <- c(
      boot_alpha/2,
      1 - boot_alpha/2
    )
    boot_ci <- function(x) {
      quantile(x, probs=boot_ci_probs)
    }
    ci_OE <- boot_ci(boot_OE)
    ci_IDE <- boot_ci(boot_IDE)
    ci_IIE <- boot_ci(boot_IIE)
    ci_CDE <- boot_ci(boot_CDE)
    
    # compute two-tailed bootstrap p-values
    boot_pval <- function(x) {
      2 * min(
        mean(x < 0),
        mean(x > 0)
      )
    }
    pvalue_OE <- boot_pval(boot_OE)
    pvalue_IDE <- boot_pval(boot_IDE)
    pvalue_IIE <- boot_pval(boot_IIE)
    pvalue_CDE <- boot_pval(boot_CDE)
    
    # compile bootstrap results
    boot_out <- list(
      ci_OE = ci_OE,
      ci_IDE = ci_IDE,
      ci_IIE = ci_IIE,
      ci_CDE = ci_CDE,
      pvalue_OE = pvalue_OE,
      pvalue_IDE = pvalue_IDE,
      pvalue_IIE = pvalue_IIE,
      pvalue_CDE = pvalue_CDE,
      boot_OE = boot_OE,
      boot_IDE = boot_IDE,
      boot_IIE = boot_IIE,
      boot_CDE = boot_CDE
    )
  }
  
  
  # final output
  out <- est
  if (boot) {
    out <- append(out, boot_out)
  }
  return(out)
}


# Run custom simulation function
# ----------------------------------------
set.seed(60657)
out_sim <- custom_sim(
  data = plow,
  D = D,
  M = M,
  Y = Y,
  L = L,
  m = m,
  L_formula = formula_L,
  M_formula = formula_M,
  Y_formula = formula_Y,
  n_sims = n_sims,
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
#  IPW ESTIMATOR  #
#-----------------#
# Define inner custom IPW function
# ----------------------------------------
custom_ipwvent_inner <- function(
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
    censor_high = 0.99,
    minimal = FALSE
) {
  # fit specified models
  ## estimate f(D|C)
  D_model <- glm(
    D_formula,
    data = data,
    family = binomial(link = "logit")
  )
  ## estimate q(L|C,D)
  L_model <- vgam(
    L_formula,
    family = cumulative(parallel = TRUE, reverse = FALSE),
    data = data
  )
  ## estimate g(M|C,D,L)
  M_model <- lm(
    M_formula,
    data = data
  )
  
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
      if (sum(data[[L]]==l)==0) {
        p <- rep(0, nrow(temp_data))
      }
      else {
        p <- predict(L_model, newdata = temp_data, type = "response")[,as.character(l)]
      }
      return(p)
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
    ((pM_D0L1C * pL1_D0C) + (pM_D0L2C * pL2_D0C) + (pM_D0L3C * pL3_D0C)) /
      (pD0_C * pM_D0LC)
    ,
    0
  )
  w2 <- ifelse(
    group_D1,
    ((pM_D1L1C * pL1_D1C) + (pM_D1L2C * pL2_D1C) + (pM_D1L3C * pL3_D1C)) /
      (pD1_C * pM_D1LC)
    ,
    0
  )
  w3 <- ifelse(
    group_D1,
    ((pM_D0L1C * pL1_D0C) + (pM_D0L2C * pL2_D0C) + (pM_D0L3C * pL3_D0C)) /
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
  if (minimal) {
    out <- list(
      OE = OE,
      IDE = IDE,
      IIE = IIE,
      CDE = CDE
    )
  }
  else {
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
  }
  return(out)
}


# Define outer custom IPW function (bootstrapping the inner function)
# ----------------------------------------
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
    censor_high = 0.99,
    boot = FALSE,
    boot_reps = 1000,
    boot_conf_level = 0.95,
    boot_seed = NULL,
    boot_parallel = FALSE,
    boot_cores = max(c(parallel::detectCores()-2,1))
) {
  # load data
  data_outer <- data
  
  
  # create adjusted boot_parallel logical
  boot_parallel_rev <- ifelse(boot_cores>1, boot_parallel, FALSE)
  
  
  # preliminary error/warning checks for the bootstrap
  if (boot) {
    if (boot_parallel & boot_cores==1) {
      warning(paste(strwrap("Warning: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but you do not have enough cores available for parallelization. The bootstrap will proceed without parallelization."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("doParallel", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'doParallel' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("doRNG", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'doRNG' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
    if (boot_parallel_rev & !requireNamespace("foreach", quietly = TRUE)) {
      stop(paste(strwrap("Error: You requested a parallelized bootstrap (boot=TRUE and boot_parallel=TRUE), but the required package 'foreach' has not been installed. Please install this package if you wish to run a parallelized bootstrap."), collapse = "\n"))
    }
  }
  
  
  # D and L values
  # (Defining these in the outer function in case the bootstrap resampling 
  # fails to sample one of the D or L values.)
  D_values <- unique(data_outer[[D]])
  L_values <- unique(data_outer[[L]])
  cross_DL <- expand.grid(D_val = D_values, L_val = L_values)
  
  
  # reset the environment of the inner function, so that it can access the 
  # D_values, L_values, and cross_DL objects
  environment(custom_ipwvent_inner) <- environment()
  
  
  # compute point estimates
  est <- custom_ipwvent_inner(
    data = data_outer,
    D = D,
    M = M,
    Y = Y,
    L = L,
    m = m,
    D_formula = D_formula,
    L_formula = L_formula,
    M_formula = M_formula,
    stabilize = stabilize,
    censor = censor,
    censor_low = censor_low,
    censor_high = censor_high,
    minimal = FALSE
  )
  
  
  # bootstrap, if requested
  if (boot) {
    # bootstrap function
    boot_fnc <- function() {
      # sample from the data with replacement
      boot_data <- data_outer[sample(nrow(data_outer), size = nrow(data_outer), replace = TRUE), ]
      
      # compute point estimates in the replicate sample
      custom_ipwvent_inner(
        data = boot_data,
        D = D,
        M = M,
        Y = Y,
        L = L,
        m = m,
        D_formula = D_formula,
        L_formula = L_formula,
        M_formula = M_formula,
        stabilize = stabilize,
        censor = censor,
        censor_low = censor_low,
        censor_high = censor_high,
        minimal = TRUE
      )
    }
    
    # parallelization prep, if parallelization requested
    if (boot_parallel_rev) {
      x_cluster <- parallel::makeCluster(boot_cores, type="PSOCK")
      doParallel::registerDoParallel(cl=x_cluster)
      parallel::clusterExport(
        cl = x_cluster, 
        varlist = c("custom_ipwvent_inner", "trimQ", "vgam", "cumulative"),
        envir = environment()
      )
      `%dopar%` <- foreach::`%dopar%`
    }
    
    # set seed
    if (!is.null(boot_seed)) {
      set.seed(boot_seed)
      if (boot_parallel) {
        doRNG::registerDoRNG(boot_seed)
      }
    }
    
    # compute estimates for each replicate sample
    if (boot_parallel_rev) {
      boot_res <- foreach::foreach(i = 1:boot_reps, .combine = comb_list_vec) %dopar% {
        boot_fnc()
      }
      boot_OE <- boot_res$OE
      boot_IDE <- boot_res$IDE
      boot_IIE <- boot_res$IIE
      boot_CDE <- boot_res$CDE
    }
    else {
      boot_OE <- rep(NA_real_, boot_reps)
      boot_IDE <- rep(NA_real_, boot_reps)
      boot_IIE <- rep(NA_real_, boot_reps)
      boot_CDE <- rep(NA_real_, boot_reps)
      for (i in seq_len(boot_reps)) {
        boot_iter <- boot_fnc()
        boot_OE[i] <- boot_iter$OE
        boot_IDE[i] <- boot_iter$IDE
        boot_IIE[i] <- boot_iter$IIE
        boot_CDE[i] <- boot_iter$CDE
      }
    }
    
    # clean up
    if (boot_parallel_rev) {
      parallel::stopCluster(x_cluster)
      rm(x_cluster)
    }
    
    # warning messages for missing estimates
    n_miss_OE <- sum(is.na(boot_OE))
    n_miss_IDE <- sum(is.na(boot_IDE))
    n_miss_IIE <- sum(is.na(boot_IIE))
    n_miss_CDE <- sum(is.na(boot_CDE))
    if (n_miss_OE>0) {
      warning(paste(strwrap(paste0("Warning: ", n_miss_OE, " of ", boot_reps, " bootstrap replicate samples are missing an estimate for the OE.")), collapse = "\n"))
    }
    if (n_miss_IDE>0) {
      warning(paste(strwrap(paste0("Warning: ", n_miss_IDE, " of ", boot_reps, " bootstrap replicate samples are missing an estimate for the IDE.")), collapse = "\n"))
    }
    if (n_miss_IIE>0) {
      warning(paste(strwrap(paste0("Warning: ", n_miss_IIE, " of ", boot_reps, " bootstrap replicate samples are missing an estimate for the IIE.")), collapse = "\n"))
    }
    if (n_miss_CDE>0) {
      warning(paste(strwrap(paste0("Warning: ", n_miss_CDE, " of ", boot_reps, " bootstrap replicate samples are missing an estimate for the CDE.")), collapse = "\n"))
    }
    
    # compute bootstrap confidence intervals 
    # from percentiles of the bootstrap distributions
    boot_alpha <- 1 - boot_conf_level
    boot_ci_probs <- c(
      boot_alpha/2,
      1 - boot_alpha/2
    )
    boot_ci <- function(x) {
      quantile(x, probs=boot_ci_probs, na.rm=TRUE)
    }
    ci_OE <- boot_ci(boot_OE)
    ci_IDE <- boot_ci(boot_IDE)
    ci_IIE <- boot_ci(boot_IIE)
    ci_CDE <- boot_ci(boot_CDE)
    
    # compute two-tailed bootstrap p-values
    boot_pval <- function(x) {
      2 * min(
        mean(x < 0, na.rm = TRUE),
        mean(x > 0, na.rm = TRUE)
      )
    }
    pvalue_OE <- boot_pval(boot_OE)
    pvalue_IDE <- boot_pval(boot_IDE)
    pvalue_IIE <- boot_pval(boot_IIE)
    pvalue_CDE <- boot_pval(boot_CDE)
    
    # compile bootstrap results
    boot_out <- list(
      ci_OE = ci_OE,
      ci_IDE = ci_IDE,
      ci_IIE = ci_IIE,
      ci_CDE = ci_CDE,
      pvalue_OE = pvalue_OE,
      pvalue_IDE = pvalue_IDE,
      pvalue_IIE = pvalue_IIE,
      pvalue_CDE = pvalue_CDE,
      boot_OE = boot_OE,
      boot_IDE = boot_IDE,
      boot_IIE = boot_IIE,
      boot_CDE = boot_CDE
    )
  }
  
  
  # final output
  out <- est
  if (boot) {
    out <- append(out, boot_out)
  }
  return(out)
}


# Run custom IPW function
# ----------------------------------------
out_ipw <- custom_ipwvent(
  data = plow,
  D = D,
  M = M,
  Y = Y,
  L = L,
  m = m,
  D_formula = formula_D,
  L_formula = formula_L,
  M_formula = formula2_M,
  censor_low = 0.02,
  censor_high = 0.98,
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




#-------------------#
#  COLLATE RESULTS  #
#-------------------#
master <- data.frame(
  param = c("OE(1,0)", "IDE(1,0)", "NIE(1,0)", "CDE(1,0,1.8K)"),
  
  # RWR
  rwr_est = c(
    out_rwr$OE,
    out_rwr$IDE,
    out_rwr$IIE,
    out_rwr$CDE
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
  sim_est = c(
    out_sim$OE,
    out_sim$IDE,
    out_sim$IIE,
    out_sim$CDE
  ),
  sim_ci_low = c(
    out_sim$ci_OE[1],
    out_sim$ci_IDE[1],
    out_sim$ci_IIE[1],
    out_sim$ci_CDE[1]
  ),
  sim_ci_high = c(
    out_sim$ci_OE[2],
    out_sim$ci_IDE[2],
    out_sim$ci_IIE[2],
    out_sim$ci_CDE[2]
  ),
  
  # IPW
  ipw_est = c(
    out_ipw$OE,
    out_ipw$IDE,
    out_ipw$IIE,
    out_ipw$CDE
  ),
  ipw_ci_low = c(
    out_ipw$ci_OE[1],
    out_ipw$ci_IDE[1],
    out_ipw$ci_IIE[1],
    out_ipw$ci_CDE[1]
  ),
  ipw_ci_high = c(
    out_ipw$ci_OE[2],
    out_ipw$ci_IDE[2],
    out_ipw$ci_IIE[2],
    out_ipw$ci_CDE[2]
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




#-----------------#
#  CREATE FIGURE  #
#-----------------#
master |>
  pivot_longer(
    cols = !param,
    names_to = c("stub", ".value"),
    names_pattern = "([^_]+)_(.+)"
  ) |>
  mutate(
    method = factor(
      case_when(
        stub=="rwr" ~ "RWR Estimates",
        stub=="sim" ~ "Simulation Estimates",
        stub=="ipw" ~ "IPW Estimates"
      ),
      levels = c("RWR Estimates", "Simulation Estimates", "IPW Estimates")
    )
  ) |>
  ggplot(aes(x = est, y = param)) +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.1) +
  facet_wrap(vars(method), ncol = 3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(-0.12, 0.08, by = 0.02)) +
  xlab("Difference in Proportions") +
  ylab("Estimand") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(panel.spacing = grid::unit(1, "lines"))

# Save the figure
ggsave(
  paste0(dir_fig, "/", title, ".png"),
  height = 4.5,
  width = 9,
  units = "in",
  dpi = 600
)


# Close log
sink()

