#----------Preliminaries----------#
rm(list = ls())
chapter <- "ch6"
title <- "table_6-1"

# Specify the root directory:
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text" 

# Define subdirectories for logs and figures:
dir_log <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log, "/", title, "_log.txt")

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

# Script:      .../code/ch6/table_6-1.R

# Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.RDS

# Outputs:     .../code/ch6/_LOGS/table_6-1_log.txt

# Description: Replicates Chapter 6, Table 6.1: Estimates for the Average Total Effect 
#              of College Attendance on Depression (CES-D scores) from the NLSY.
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#
packages <-
  c(
    "survey", 
    "gbm",
    "ranger",
    "glmnet",
    "rsample",
    "caret",
    "rlang",
    "tidyverse",
    "Hmisc",
    "SuperLearner",
    "scales",
    "haven",
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

# mediators
M <- list(
  "ever_unemp_age3539",
  "log_faminc_adj_age3539"
)

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
  "cesd_age40",
  D,
  unlist(M),
  C
)

# number of bootstrap replications
n_reps <- 2000

# set seed
seed <- 02138

#-----------------------------#
#        PREPARE DATA         #
#-----------------------------#

nlsy_raw <- readRDS(
  url("https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.RDS")
)

df <- 
  nlsy_raw[complete.cases(nlsy_raw[, key_vars]),] |>
  mutate(std_cesd_age40 = (cesd_age40 - mean(cesd_age40)) / sd(cesd_age40))

#------------------------------------------------------------------------------#
#                            REPLICATE TABLE 6.1                               #
#------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
#                  Section 1: Parametric MR estimation of ATE                 #
#-----------------------------------------------------------------------------#

# Specify form of the models for the outcome and exposure
D_model_f <- as.formula(paste(D, " ~ ", paste(C, collapse= "+")))
Y_model_f <- as.formula(paste(Y, " ~ ", paste(c(C, D), collapse= "+")))

# Define function for implementing parametric MR estimator for ATE

trimQ <- function(x, low = 0.01, high = 0.99) {
  min <- quantile(x, low)
  max <- quantile(x, high)
  
  x[x<min] <- min
  x[x>max] <- max
  x
}

paraMR <- function(
    D,
    Y,
    Y_model_f,
    D_model_f,
    data,
    boot_seed,
    boot_reps = 2000
){
  
  paraMR_inner <- function(
    D,
    Y,
    Y_model_f,
    D_model_f,
    df
  ){
     
    # Design matrix for outcome model: μD(C) = E(Y|D,C)
    dm_mu <- model.matrix(Y_model_f, data = df)[, -1] %>% as_tibble()
    
    # Design matrices for outcome prediction
    dm_mu1 <- model.matrix(Y_model_f, data = mutate(df, !!sym(D) := 1)) %>% as_tibble() 
    dm_mu0 <- model.matrix(Y_model_f, data = mutate(df, !!sym(D) := 0)) %>% as_tibble() 
  
    # Estimate IP weights
    D_model <- glm(
     D_model_f, 
     family = binomial("logit"), 
     data = df
    )
  
    df <- df %>% 
      mutate(
        pi_hat = D_model$fitted.values,
        W0 = as.double(.data[[D]] == 0)/(1 - pi_hat),
        W1 = as.double(.data[[D]] == 1)/pi_hat
      )
  
    df$W0[df[[D]] == 0] <- trimQ(df$W0[df[[D]] == 0])
    df$W1[df[[D]] == 1] <- trimQ(df$W1[df[[D]] == 1])
    
  # Impute outcomes
  Y_model <- lm(Y_model_f, data = df)
  
  df <- df %>%
    mutate(
      mu0_hat = predict(Y_model, newdata = dm_mu0),
      mu1_hat = predict(Y_model, newdata = dm_mu1)
    )
  
  # Compute MR estimates
  final_results <- df %>%
     mutate(
      DR_1 = mu1_hat + W1 * (!!sym(Y) - mu1_hat),
      DR_0 = mu0_hat + W0 * (!!sym(Y) - mu0_hat),
      DR_ATE = DR_1 - DR_0,
      IPW_1 = W1 * !!sym(Y),
      IPW_0 = W0 * !!sym(Y),
      IPW_ATE = IPW_1 - IPW_0,
      RI_1 = mu1_hat,
      RI_0 = mu0_hat,
      RI_ATE = RI_1 - RI_0
    ) %>%
    summarise_at(
      vars(contains(
        c("DR_", "IPW_", "RI_"))), 
      list(est = ~ wtd.mean(.x))) %>%
    pivot_longer(everything()) %>%
    separate(name, into = c("estimator","estimand", "measure")) %>%
    pivot_wider(names_from = measure, values_from = value)
  }
  
  baseline_estimates <- 
    paraMR_inner(
      D ,
      Y,
      Y_model_f,
      D_model_f,
      df
    )
  
  # Compute bootstrap distribution #
  set.seed(boot_seed)
  
  boot_rst <-
    lapply(
      seq(1, boot_reps),
      function(i){
        if (i %% 100 == 0){
          cat("bootstrap sample ", i, "\n")
        }
        dfi <- df %>% sample_frac(replace = TRUE)
        baseline_estimates <- 
          paraMR_inner(
            D ,
            Y,
            Y_model_f,
            D_model_f,
            dfi
          )
      }
    )
  
  final_df <-
    left_join(
      bind_rows(boot_rst) %>%
      group_by(estimator, estimand) %>%
      summarise(
        lower = quantile(est,0.025),
        upper = quantile(est,0.975),
        .groups = "keep"
      ),
      baseline_estimates,
      by = c("estimand","estimator")
    ) %>%
    mutate(across(everything(), ~ round(as.numeric(.x), 3))) %>%
    mutate(
      results = paste0(est,"[", lower, ",", upper, "]")
    ) %>%
    dplyr::select(
      estimator,
      estimand,
      results
    ) %>%
    pivot_wider(
      names_from = estimator,  
      values_from = results     
    ) %>%
    rename(
      Estimand = estimand,
      `Regression Imputation` = RI,
      Weighting = IPW,
      `Parametric DR` = DR
    )
  
  return(final_df)
  
}

paraMR_results <-
  paraMR(
    D,
    Y,
    Y_model_f,
    D_model_f,
    data = df,
    boot_seed = seed,
    boot_reps = 2000
  ) %>%
  dplyr::select(
    "Estimand",
    "Regression Imputation",
    "Weighting",
    "Parametric DR"
  )

#-----------------------------------------------------------------------------#
#                  Section 2: DML estimation of ATE                           #
#-----------------------------------------------------------------------------#

DMLest <- function(
    D,
    Y,
    Y_model_f,
    D_model_f,
    df,
    K,
    seed
){
  
    # Design matrix for exposure model: πD(C)
    dm_pi <- model.matrix(D_model_f, data = df)[, -1] %>% as_tibble()
    
    # Design matrix for outcome model: μD(C)
    dm_mu <- model.matrix(Y_model_f, data = df)[, -1] %>% as_tibble()
    
    # Design matrix for outcome prediction
    dm_mu1 <- model.matrix(Y_model_f, data = mutate(df, !!sym(D) := 1))[, -1] %>% as_tibble() 
    dm_mu0 <- model.matrix(Y_model_f, data = mutate(df, !!sym(D) := 0))[, -1] %>% as_tibble() 
    
    # Initialize repeated cross-fitting
    set.seed(seed)

    cf_folds <- createFolds(df$std_cesd_age40, K)
    
    cf_rst_lst <-
      lapply(
        cf_folds,
        function(predict_fold){
          
          # Define the training and estimation data
          train_k <- df[-predict_fold,]
          pred_k <- df[predict_fold,]

          train_k_mu <- dm_mu[-predict_fold,]
          pred_k_mu <- dm_mu[predict_fold,]
          pred_k_mu0 <- dm_mu0[predict_fold,]
          pred_k_mu1 <- dm_mu1[predict_fold,]
          Y_k <- train_k[[Y]]

          train_k_pi <- dm_pi[-predict_fold,]
          pred_k_pi <- dm_pi[predict_fold,]
          D_k <- train_k[[D]]

          # Train the exposure model
          D_model <- SuperLearner(
            Y          = D_k,
            X          = train_k_pi,
            newX       = pred_k_pi,
            family     = binomial(),
            SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
            control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
            cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
          )
          
          # Estimate weights
          pred_k <-
            pred_k %>% 
            mutate(
              pi_hat = as.numeric(D_model$SL.predict),
              W0 = as.double(.data[[D]] == 0)/(1 - pi_hat),
              W1 = as.double(.data[[D]] == 1)/pi_hat
            )
          
          # Train the outcome model
          Y_model <- SuperLearner(
            Y          = Y_k,
            X          = train_k_mu,
            family     = gaussian(),
            SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
            control    = list(saveFitLibrary = TRUE),
            cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
          )
          
          # Impute outcomes
          pred_k$mu0_hat <- as.numeric(predict.SuperLearner(Y_model, newdata = pred_k_mu0)$pred)
          pred_k$mu1_hat <- as.numeric(predict.SuperLearner(Y_model, newdata = pred_k_mu1)$pred)
          
          return(pred_k)
          
        }
      )
    
    final_df <- bind_rows(cf_rst_lst)
    final_df$W0[final_df[[D]] == 0] <- trimQ(final_df$W0[final_df[[D]] == 0])
    final_df$W1[final_df[[D]] == 1] <- trimQ(final_df$W1[final_df[[D]] == 1])
    
    # Compute DML estimates
    final_results <-
      final_df %>%
      mutate(
        DR_1 = mu1_hat + W1 * (!!sym(Y) - mu1_hat),
        DR_0 = mu0_hat + W0 * (!!sym(Y) - mu0_hat),
        DR_ATE = DR_1 - DR_0
      ) %>% 
      summarise_at(
        vars(contains(c("DR_"))), 
        list(
          est = ~ wtd.mean(.x),
          se  = ~ sqrt(wtd.var(.x) / length(.x))
        )
      ) %>%
      mutate(across(everything(), ~ round(as.numeric(.x), 3))) %>%
      pivot_longer(everything()) %>%
      separate(name, into = c("estimator","estimand", "measure")) %>%
      pivot_wider(names_from = measure, values_from = value) %>%
      mutate(
        result = paste0(
          round(est,3),
          "[", 
          round(est - 1.96 * se, 3),
          ",", 
          round(est + 1.96 * se, 3), 
          "]")
      ) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      dplyr::select(
        estimator,
        estimand,
        result
      ) %>%
      pivot_wider(
        names_from = estimator,  
        values_from = result    
      ) %>%
      rename(
        Estimand = estimand,
        DML = DR
      )
    
    return(final_results)
    
}

DMLest_results <-
  DMLest(
    D,
    Y,
    Y_model_f,
    D_model_f,
    df = df,
    K=5,
    seed = 02138
  )

#-------------------------------#
#        COLLATE RESULTS        #
#-------------------------------#

master <-
  left_join(
    paraMR_results,
    DMLest_results,
    by = "Estimand"
  )

# Open log
sink(log_path, split = TRUE)

width_curr <- getOption("width")
options(width = 300)
print(master)

# Close log
sink()
