rm(list = ls())
library(survey)
library(gbm)
library(ranger)
library(glmnet)
library(rsample)
library(caret)
library(rlang)
library(tidyverse)
library(Hmisc)
library(SuperLearner)
library(scales)
source("../ch5/utils.R")
set.seed(02138)

# set my ggplot theme
mytheme <- theme_minimal(base_size = 18) + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(color = "grey30"))
theme_set(mytheme)

##office
datadir <- "../../data/" 

##input data
nlsy_raw <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

a <- "att22"
m1 <- "ever_unemp_age3539"
m2 <- "log_faminc_adj_age3539"
y <- "std_cesd_age40"
x <- c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")

nlsy <- nlsy_raw[complete.cases(nlsy_raw[, c(a, "cesd_age40", m1, m2, x)]),] %>% 
  mutate(std_cesd_age40 = as.numeric(scale(cesd_age40)))

df <- nlsy
n <- nrow(df)

##########################################################
# Formulas for treatment and outcome models
##########################################################

a0_form <- as.formula(paste(a, " ~ ", paste(x, collapse= "+")))
y0_form <- as.formula(paste(y, " ~ ", paste(c(x, a), collapse= "+")))

##########################################################
# Main analyses
##########################################################

estimands <- expand.grid(c(0, 1)) %>%
  `colnames<-`(c("a"))

S <- nrow(estimands)

K <- 5

df_p0 <- model.matrix(a0_form, data = df)[, -1] %>% as_tibble()
df_mu0 <- model.matrix(y0_form, data = df)[, -1] %>% as_tibble()

df_mu0n <- model.matrix(y0_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_mu0y <- model.matrix(y0_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()

# create cross-fitting split

cf_fold <- createFolds(df$std_cesd_age40, K)

main_list <- vector(mode = "list", K)

for(k in 1:K){
  
  cat(" cross-fitting fold ", k, "\n")
  
  #################################################
  # Design matrices for different models
  #################################################
  
  # auxiliary and main data
  aux <- df[-cf_fold[[k]], ]
  main <- df[cf_fold[[k]], ]
  
  aux_p0 <- df_p0[-cf_fold[[k]], ]
  main_p0 <- df_p0[cf_fold[[k]], ]
  
  aux_mu0 <- df_mu0[-cf_fold[[k]], ]
  main_mu0 <- df_mu0[cf_fold[[k]], ]
  
  #################################################
  # Treatment Models
  #################################################
  
  p0_sl <- SuperLearner(
    Y          = aux$att22,
    X          = aux_p0,
    newX       = df_p0,
    family     = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  df <- df %>% mutate(
    p0_fit = p0_sl$SL.predict,
    w0_n = I(att22 == 0)/(1 - p0_fit),
    w0_y = I(att22 == 1)/p0_fit,
  )
  
  #################################################
  # Outcome models
  #################################################
  
  mu0_sl <- SuperLearner(
    Y          = aux$std_cesd_age40,
    X          = aux_mu0,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu0_fit_n <- predict.SuperLearner(mu0_sl, newdata = df_mu0n)$pred
  df$mu0_fit_y <- predict.SuperLearner(mu0_sl, newdata = df_mu0y)$pred
  
  main_list[[k]] <- df[cf_fold[[k]], ]
}

main_df <- reduce(main_list, bind_rows)

for (s in 1:S){
  
  a <- estimands$a[[s]]
  
  main_df <- main_df %>%
    mutate(
      
      wt0_deno = a * p0_fit + (1 - a) * (1 - p0_fit),
      
      !!sym(paste0("w0_", a)) := as.double(att22==a)/wt0_deno,
      
    )
  
  main_df[main_df$att22 == a, paste0("w0_", a)] <- trimQ(main_df[main_df$att22 == a, paste0("w0_", a)])
  
  main_df <- main_df %>% 
    mutate(
    
      !!sym(paste0("mu0fit_", a)) := a * mu0_fit_y + (1 - a) * mu0_fit_n,
      
      !!sym(paste0("www_", a)) := !!sym(paste0("w0_", a)) * std_cesd_age40,
      !!sym(paste0("iii_", a)) := !!sym(paste0("mu0fit_", a)),
      !!sym(paste0("eif_", a)) := !!sym(paste0("w0_", a)) *
        (std_cesd_age40 - !!sym(paste0("mu0fit_", a))) +
        !!sym(paste0("mu0fit_", a))
      
    )
  
}

out_df <- main_df %>%
  mutate(www_ate = www_1 - www_0,
         iii_ate = iii_1 - iii_0,
         eif_ate = eif_1 - eif_0) %>%
  summarise_at(vars(contains("eif_")), list(est = ~ wtd.mean(.x),
                                            se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator","estimand", "measure")) %>%
  pivot_wider(names_from = measure, values_from = value)

table6_1np <-  out_df %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se) %>% 
  mutate_at(c("est", "se", "lower", "upper"), ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")")) %>% 
  mutate(out = paste(est, intv)) %>% 
  dplyr::select(-lower, -upper, -intv) 

write_csv(table6_1np, file = "table6-1np.csv")

save.image(file = "table6-1np.RData")

