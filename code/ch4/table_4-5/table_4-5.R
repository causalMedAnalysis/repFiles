# Preliminaries
chapter <- "ch4"
title <- "table_4-5"
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

#-------------------------------------------------------------------------------
# Causal Mediation Analysis Replication Files

# GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main

# Script:      .../code/ch4/table_4-5.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta

# Outputs:     .../code/ch4/_LOGS/table_4-5_log.txt

# Description: Replicates Chapter 4, Table 4.5: Bootstrap Inferential Statistics
#              for the Interventional Effects of College Attendance on CES-D
#              Scores Computed from the NLSY.
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
Y <- "std_cesd_age40"

# exposure
D <- "att22"

# mediator
M <- "log_faminc_adj_age3539"

# exposure-induced confounder
L <- "ever_unemp_age3539"

# baseline confounders
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

# number of simulations for simulation estimator
n_sims <- 2000

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

#------------------#
#  MODEL FORMULAE  #
#------------------#
# D model: Additive
# L model: With D x C interactions
# M model (for RWR and simulation): With D x C interactions
# M model (for IPW): With D x C, D x L interactions (and L main effect)
# Y model: With D x M, D x C, M x C, M x L interactions
# Note that M here is the log form

# D model formula
predictors_D <- paste(C, collapse = " + ")
(formula_D_string <- paste(D, "~", predictors_D))
formula_D <- as.formula(formula_D_string)

# L and M model formulae
## main effects
predictors_LM <- paste(c(D,C), collapse = " + ")
## D x C interactions
predictors_LM <- paste(
  predictors_LM,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## full formula
(formula_L_string <- paste(L, "~", predictors_LM))
(formula_M_string <- paste(M, "~", predictors_LM))
formula_L <- as.formula(formula_L_string)
formula_M <- as.formula(formula_M_string)

# M model formula for IPW
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

# Y model formula
## main effects
predictors_Y <- paste(c(D,M,L,C), collapse = " + ")
## D x M interaction
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(D, M, sep = ":", collapse = " + ")
)
## D x C interactions
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(D, C, sep = ":", collapse = " + ")
)
## M x C interactions
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(M, C, sep = ":", collapse = " + ")
)
## M x L interaction
predictors_Y <- paste(
  predictors_Y,
  "+",
  paste(M, L, sep = ":", collapse = " + ")
)
## full formula
(formula_Y_string <- paste(Y, "~", predictors_Y))
formula_Y <- as.formula(formula_Y_string)

#-----------------#
#  RWR ESTIMATOR  #
#-----------------#
out_rwr <- rwrlite(
  data = nlsy,
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
)

#------------------------#
#  SIMULATION ESTIMATOR  #
#------------------------#
# Define model specifications
out_sim_specs <- list(
  ## L model
  list(
    func = "glm",
    formula = formula_L,
    args = list(family = "binomial")
  ),
  ## M model
  list(
    func = "lm",
    formula = formula_M
  ),
  ## Y model
  list(
    func = "lm",
    formula = formula_Y
  )
)

# Estimate interventional effects
out_sim <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = M,
  model_spec = out_sim_specs,
  boot = TRUE,
  boot_reps = n_reps,
  seed = 3308004
)

# Estimate CDE(1,0,ln(50K))
out_sim_cde <- medsim(
  data = nlsy,
  num_sim = n_sims,
  treatment = D,
  intv_med = paste0(M,"=",m),
  model_spec = out_sim_specs,
  boot = TRUE,
  boot_reps = n_reps,
  seed = 3308004
)

#-----------------#
#  IPW ESTIMATOR  #
#-----------------#

# Define inner custom IPW function
trimQ <- function(x, low = 0.01, high = 0.99) {
  min <- quantile(x, low)
  max <- quantile(x, high)
  
  x[x<min] <- min
  x[x>max] <- max
  x
}

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

# function to combine lists of vectors with the same structure
# (used in the parallelized bootstraps)
comb_list_vec <- function(...) {
  mapply(c, ..., SIMPLIFY = FALSE)
}

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
        varlist = c("custom_ipwvent_inner", "trimQ"),
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

# Run custom IPW function
out_ipw <- custom_ipwvent(
  data = nlsy,
  D = D,
  M = M,
  Y = Y,
  L = L,
  m = m,
  D_formula = formula_D,
  L_formula = formula_L,
  M_formula = formula2_M,
  boot = TRUE,
  boot_reps = n_reps,
  boot_seed = 3308004,
  boot_parallel = TRUE
)

#-------------------#
#  COLLATE RESULTS  #
#-------------------#
master <- data.frame(
  param = c("OE(1,0)", "IDE(1,0)", "IIE(1,0)", "CDE(1,0,ln(50K))"),
  # RWR
  rwr_pvalue = c(
    out_rwr$pvalue_OE,
    out_rwr$pvalue_IDE,
    out_rwr$pvalue_IIE,
    out_rwr$pvalue_CDE
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
  sim_pvalue = c(
    out_sim$results[3,2],
    out_sim$results[1,2],
    out_sim$results[2,2],
    out_sim_cde$results[1,2]
  ),
  sim_ci_low = c(
    out_sim$results[3,3],
    out_sim$results[1,3],
    out_sim$results[2,3],
    out_sim_cde$results[1,3]
  ),
  sim_ci_high = c(
    out_sim$results[3,4],
    out_sim$results[1,4],
    out_sim$results[2,4],
    out_sim_cde$results[1,4]
  ),
  # IPW
  ipw_pvalue = c(
    out_ipw$pvalue_OE,
    out_ipw$pvalue_IDE,
    out_ipw$pvalue_IIE,
    out_ipw$pvalue_CDE
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
