# Preliminaries
chapter <- "ch3"
title <- "table_3-7"
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")

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

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch3/table_3-7.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/JOBSII/Jobs-NoMiss-Binary.dta

# Outputs:     .../code/ch3/_LOGS/table_3-7_log.txt

# Description: Replicates Chapter 3, Table 3-7: Total, Direct, and Indirect
#              Effects of Job Training on Employment as Estimated from JOBSII Study.
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#
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

# number of simulations for simulation estimator
n_sims <- 1000

# number of bootstrap replications
n_reps <- 2000

#----------------#
#  PREPARE DATA  #
#----------------#
jobs_raw <- read_stata(
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/JOBSII/Jobs-NoMiss-Binary.dta"
)

jobs <- jobs_raw |>
  zap_labels()

#--------------------------#
#  LINEAR MODEL ESTIMATOR  #
#--------------------------#
# Linear model with D x M interaction

out_lin <- linmed(
  data = jobs,
  D = D,
  M = M,
  Y = Y,
  C = C,
  m = m,
  interaction_DM = TRUE,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------------------------------------#
#  SIMULATION & REGRESSION IMPUTATION ESTIMATORS  #
#-------------------------------------------------#
# M model: Additive linear model
# Y model: Logit model with D x M interaction

# Mediator model formula
predictors_M <- paste(c(D,C), collapse = " + ")
formula_M_string <- paste(M, "~", predictors_M)

# Outcome model formula
## main effects
predictors_Y <- paste(c(D,M,C), collapse = " + ")

## D x M interaction
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(D, M, sep = ":", collapse = " + ")
)

## full formula
formula_Y_string <- paste(Y, "~", predictors_Y)

# Define model specifications
sim_specs <- list(
  list(
    func = "lm",
    formula = as.formula(formula_M_string)
  ),
  list(
    func = "glm",
    formula = as.formula(formula_Y_string),
    args = list(family = "binomial")
  )
)

# Estimate ATE(1,0), NDE(1,0), and NIE(1,0) by simulation estimator
out_sim <- medsim(
  data = jobs,
  num_sim = n_sims,
  treatment = D,
  intv_med = NULL,
  model_spec = sim_specs,
  boot = TRUE,
  boot_reps = n_reps,
  seed = 3308004
)

# Estimate CDE(1,0,4) by regression imputation estimator
mod_Y <- glm(
  work1 ~ treat*job_seek + econ_hard + sex + age + nonwhite + educ + income,
  family = binomial(link = "logit"),
  data = jobs
)

out_imp_cde <- impcde(
  data = jobs,
  model_y = mod_Y,
  D = D,
  M = M,
  m = m,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-----------------#
#  IPW ESTIMATOR  #
#-----------------#
# Additive logit models

# D model 1 formula: f(D|C)
predictors1_D <- paste(C, collapse = " + ")
formula1_D_string <- paste(D, "~", predictors1_D)

# D model 2 formula: s(D|C,M)
predictors2_D <- paste(c(M,C), collapse = " + ")
formula2_D_string <- paste(D, "~", predictors2_D)

# M model formula: g(M|C,D)
#formula_M_string defined above

# Estimate ATE(1,0), NDE(1,0), NIE(1,0)
out_ipw <- ipwmed(
  data = jobs,
  D = D,
  M = M,
  Y = Y,
  formula1_string = formula1_D_string,
  formula2_string = formula2_D_string,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#---------------------#
#  IPW CDE ESTIMATOR  #
#---------------------#
# The IPW CDE estimator function ipwcde() only supports binary exposures and
# mediators. But the job_seek mediator is not binary. We will create a custom
# IPW CDE function here, treating the mediator as pseudo-continuous.

# Define inner custom IPW function
custom_ipwcde_inner <- function(
    data,
    D,
    M,
    Y,
    m = 0,
    formula_D_string,
    formula_M_string,
    stabilize = TRUE,
    censor = TRUE,
    censor_low = 0.01,
    censor_high = 0.99,
    minimal = FALSE
) {
  
  # preliminaries
  d <- 1
  dstar <- 0

  # load data
  df <- data

  # fit specified models
  d_model <- glm(
    as.formula(formula_D_string),
    data = df,
    family = binomial(link = "logit")
  )
  m_model <- lm(
    as.formula(formula_M_string),
    data = df
  )

  # additionally fit mediator model without covariates
  m_model_no_cov <- lm(
    as.formula(paste0(M,"~",D)),
    data = df
  )

  # predict exposure and mediator probabilities
  ps_D1_C <- predict(d_model, newdata = df, type = "response")
  ps_D_C <- ifelse(as.logical(df[[D]]), ps_D1_C, 1-ps_D1_C)
  marg_prob_D1 <- mean(df[[D]])
  marg_prob_D <- ifelse(
    as.logical(df[[D]]),
    marg_prob_D1,
    1 - marg_prob_D1
  )
  Ehat_M_CD <- predict(m_model, newdata = df)
  sighat_M_CD <- sqrt(mean(m_model$residuals^2))
  ps_M_CD <- dnorm(
    df[[M]],
    mean = Ehat_M_CD,
    sd = sighat_M_CD
  )
  Ehat_M_D <- predict(m_model_no_cov, newdata = df)
  sighat_M_D <- sqrt(mean(m_model_no_cov$residuals^2))
  marg_dens_M_D <- dnorm(
    df[[M]],
    mean = Ehat_M_D,
    sd = sighat_M_D
  )

  # create IPWs
  w4 <- 1 / (ps_M_CD * ps_D_C)

  # stabilize IPWs
  if (stabilize) {
    w4 <- w4 * marg_dens_M_D * marg_prob_D
  }

  # censor IPWs
  if (censor) {
    w4 <- trimQ(w4, low = censor_low, high = censor_high)
  }

  # estimate effects
  y_model <- lm(
    as.formula(paste0(Y,"~",D,"*",M)),
    data = df,
    weights = w4
  )
  CDE <-
    y_model$coefficients[[D]] +
    y_model$coefficients[[paste0(D,":",M)]] * m

  # compile and output
  if (minimal) {
    out <- CDE
  }
  else {
    out <- list(
      CDE = CDE,
      weights = w4,
      model_d = d_model,
      model_m = m_model
    )
  }
  return(out)
}

# Define outer custom IPW function (bootstrapping the inner function)
trimQ <- function(x, low = 0.01, high = 0.99) {
  min <- quantile(x, low)
  max <- quantile(x, high)
  
  x[x<min] <- min
  x[x>max] <- max
  x
}

custom_ipwcde <- function(
    data,
    D,
    M,
    Y,
    m = 0,
    formula_D_string,
    formula_M_string,
    base_weights_name = NULL,
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

  # other error/warning checks
  if (length(M)>1) {
    stop(paste(strwrap("Error: Unlike the ipwmed() function, ipwcde() requires a single mediator. Multiple mediators are not supported."), collapse = "\n"))
  }
  if (length(m)>1) {
    stop(paste(strwrap("Error: Please specify a single value for the argument m."), collapse = "\n"))
  }
  if (!is.numeric(data_outer[[D]])) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be numeric."), collapse = "\n"))
  }
  if (!is.numeric(data_outer[[M]])) {
    stop(paste(strwrap("Error: The mediator variable (identified by the string argument M in data) must be numeric."), collapse = "\n"))
  }
  if (!is.numeric(data_outer[[Y]])) {
    stop(paste(strwrap("Error: The outcome variable (identified by the string argument Y in data) must be numeric."), collapse = "\n"))
  }
  if (any(is.na(data_outer[[D]]))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the exposure variable (identified by the string argument D in data)."), collapse = "\n"))
  }
  if (any(! data_outer[[D]] %in% c(0,1))) {
    stop(paste(strwrap("Error: The exposure variable (identified by the string argument D in data) must be a numeric variable consisting only of the values 0 or 1. There is at least one observation in the data that does not meet this criteria."), collapse = "\n"))
  }
  if (any(is.na(data_outer[[M]]))) {
    stop(paste(strwrap("Error: There is at least one observation with a missing/NA value for the mediator variable (identified by the string argument M in data)."), collapse = "\n"))
  }
  if (grepl(pattern = M, x = formula_D_string, fixed = TRUE)) {
    warning(paste(strwrap("Warning: Check whether the mediator variable is among the predictors in the formula_D_string. The mediator should not be among the predictors in the formula_D_string."), collapse = "\n"))
  }
  if (!grepl(pattern = D, x = formula_M_string, fixed = TRUE)) {
    warning(paste(strwrap("Warning: Check whether the exposure variable is among the predictors in the formula_M_string. The exposure should be among the predictors in the formula_M_string."), collapse = "\n"))
  }

  # compute point estimates
  est <- custom_ipwcde_inner(
    data = data_outer,
    D = D,
    M = M,
    Y = Y,
    m = m,
    formula_D_string = formula_D_string,
    formula_M_string = formula_M_string,
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
      custom_ipwcde_inner(
        data = boot_data,
        D = D,
        M = M,
        Y = Y,
        m = m,
        formula_D_string = formula_D_string,
        formula_M_string = formula_M_string,
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
        varlist = c("custom_ipwcde_inner", "trimQ"),
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
      boot_CDE <- foreach::foreach(i = 1:boot_reps, .combine = c) %dopar% {
        boot_fnc()
      }
    }
    else {
      boot_CDE <- rep(NA_real_, boot_reps)
      for (i in seq_len(boot_reps)) {
        boot_CDE[i] <- boot_fnc()
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
    ci_CDE <- boot_ci(boot_CDE)

    # compute two-tailed bootstrap p-values
    boot_pval <- function(x) {
      2 * min(
        mean(x < 0),
        mean(x > 0)
      )
    }
    pvalue_CDE <- boot_pval(boot_CDE)

    # compile bootstrap results
    boot_out <- list(
      ci_CDE = ci_CDE,
      pvalue_CDE = pvalue_CDE,
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
out_ipw_cde <- custom_ipwcde(
  data = jobs,
  D = D,
  M = M,
  Y = Y,
  m = m,
  formula_D_string = formula1_D_string,
  formula_M_string = formula_M_string,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------#
#  COLLATE RESULTS  #
#-------------------#
master <- data.frame(
  param = c("ATE(1,0)", "NDE(1,0)", "NIE(1,0)", "CDE(1,0,4)"),

  # linear models
  lin_est = c(
    out_lin$ATE,
    out_lin$NDE,
    out_lin$NIE,
    out_lin$CDE
  ),
  lin_ci_low = c(
    out_lin$ci_ATE[1],
    out_lin$ci_NDE[1],
    out_lin$ci_NIE[1],
    out_lin$ci_CDE[1]
  ),
  lin_ci_high = c(
    out_lin$ci_ATE[2],
    out_lin$ci_NDE[2],
    out_lin$ci_NIE[2],
    out_lin$ci_CDE[2]
  ),

  # simulation and regression imputation
  sim_est = c(
    out_sim$results[1,1],
    out_sim$results[2,1],
    out_sim$results[3,1],
    out_imp_cde$CDE
  ),
  sim_ci_low = c(
    out_sim$results[1,3],
    out_sim$results[2,3],
    out_sim$results[3,3],
    out_imp_cde$ci_CDE[1]
  ),
  sim_ci_high = c(
    out_sim$results[1,4],
    out_sim$results[2,4],
    out_sim$results[3,4],
    out_imp_cde$ci_CDE[2]
  ),

  # IPW
  ipw_est = c(
    out_ipw$ATE,
    out_ipw$NDE,
    out_ipw$NIE,
    out_ipw_cde$CDE
  ),
  ipw_ci_low = c(
    out_ipw$ci_ATE[1],
    out_ipw$ci_NDE[1],
    out_ipw$ci_NIE[1],
    out_ipw_cde$ci_CDE[1]
  ),
  ipw_ci_high = c(
    out_ipw$ci_ATE[2],
    out_ipw$ci_NDE[2],
    out_ipw$ci_NIE[2],
    out_ipw_cde$ci_CDE[2]
  )
)

# Open log
sink(log_path, split = TRUE)

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
