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
a1_form <- as.formula(paste(a, " ~ ", paste(c(x, m1), collapse= "+")))
a2_form <- as.formula(paste(a, " ~ ", paste(c(x, m1, m2), collapse= "+")))

y0_form <- as.formula(paste(y, " ~ ", paste(c(x, a), collapse= "+")))
y1_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1), collapse= "+")))
y2_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1, m2), collapse= "+")))

##########################################################
# Main analyses
##########################################################

estimands <- expand.grid(c(0, 1), c(0, 1), c(0, 1)) %>%
  `colnames<-`(c("a1", "a2", "a3"))

S <- nrow(estimands)

K <- 5
  
df_p0 <- model.matrix(a0_form, data = df)[, -1] %>% as_tibble()
df_p1 <- model.matrix(a1_form, data = df)[, -1] %>% as_tibble()
df_p2 <- model.matrix(a2_form, data = df)[, -1] %>% as_tibble()

df_mu0 <- model.matrix(y0_form, data = df)[, -1] %>% as_tibble()
df_mu1 <- model.matrix(y1_form, data = df)[, -1] %>% as_tibble()
df_mu2 <- model.matrix(y2_form, data = df)[, -1] %>% as_tibble()

df_mu2n <- model.matrix(y2_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_mu1n <- model.matrix(y1_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_mu0n <- model.matrix(y0_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()

df_mu2y <- model.matrix(y2_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()
df_mu1y <- model.matrix(y1_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()
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
  aux_p1 <- df_p1[-cf_fold[[k]], ]
  aux_p2 <- df_p2[-cf_fold[[k]], ]
  
  main_p0 <- df_p0[cf_fold[[k]], ]
  main_p1 <- df_p1[cf_fold[[k]], ]
  main_p2 <- df_p2[cf_fold[[k]], ]
  
  aux_mu0 <- df_mu0[-cf_fold[[k]], ]
  aux_mu1 <- df_mu1[-cf_fold[[k]], ]
  aux_mu2 <- df_mu2[-cf_fold[[k]], ]
  
  main_mu0 <- df_mu0[cf_fold[[k]], ]
  main_mu1 <- df_mu1[cf_fold[[k]], ]
  main_mu2 <- df_mu2[cf_fold[[k]], ]
  
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
  
  p1_sl <- SuperLearner(
    Y          = aux$att22,
    X          = aux_p1,
    newX       = df_p1,
    family     = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  p2_sl <- SuperLearner(
    Y          = aux$att22,
    X          = aux_p2,
    newX       = df_p2,
    family     = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  df <- df %>% mutate(
    p0_fit = p0_sl$SL.predict,
    p1_fit = p1_sl$SL.predict,
    p2_fit = p2_sl$SL.predict,
    w0_n = I(att22 == 0)/(1 - p0_fit),
    w0_y = I(att22 == 1)/p0_fit,
    w1_nn = I(att22 == 0)/(1 - p0_fit) * 1,
    w1_ny = I(att22 == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit,
    w1_yn = I(att22 == 0)/p0_fit * p1_fit/(1 - p1_fit),
    w1_yy = I(att22 == 1)/p0_fit * 1,
    w2_000 = I(att22 == 0)/(1 - p0_fit) * 1 * 1,
    w2_001 = I(att22 == 1)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit,
    w2_010 = I(att22 == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit),
    w2_011 = I(att22 == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1,
    w2_100 = I(att22 == 0)/p0_fit * p1_fit/(1 - p1_fit) * 1,
    w2_101 = I(att22 == 1)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit,
    w2_110 = I(att22 == 0)/p0_fit * 1 * p2_fit/(1 - p2_fit),
    w2_111 = I(att22 == 1)/p0_fit * 1 * 1,
  )
  
  #################################################
  # Outcome models
  #################################################
  
  mu2_sl <- SuperLearner(
    Y          = aux$std_cesd_age40,
    X          = aux_mu2,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu2_fit_a3n <- predict.SuperLearner(mu2_sl, newdata = df_mu2n)$pred
  df$mu2_fit_a3y <- predict.SuperLearner(mu2_sl, newdata = df_mu2y)$pred
  
  mu1_sl_a3n <- SuperLearner(
    Y          = df$mu2_fit_a3n[-cf_fold[[k]]],
    X          = aux_mu1,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu1_sl_a3y <- SuperLearner(
    Y          = df$mu2_fit_a3y[-cf_fold[[k]]],
    X          = aux_mu1,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu1_fit_a3n_a2n <- predict.SuperLearner(mu1_sl_a3n, newdata = df_mu1n)$pred
  df$mu1_fit_a3n_a2y <- predict.SuperLearner(mu1_sl_a3n, newdata = df_mu1y)$pred
  df$mu1_fit_a3y_a2n <- predict.SuperLearner(mu1_sl_a3y, newdata = df_mu1n)$pred
  df$mu1_fit_a3y_a2y <- predict.SuperLearner(mu1_sl_a3y, newdata = df_mu1y)$pred
  
  mu0_sl_a3n_a2n <- SuperLearner(
    Y          = df$mu1_fit_a3n_a2n[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a3n_a2y <- SuperLearner(
    Y          = df$mu1_fit_a3n_a2y[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a3y_a2n <- SuperLearner(
    Y          = df$mu1_fit_a3y_a2n[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a3y_a2y <- SuperLearner(
    Y          = df$mu1_fit_a3y_a2y[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu0_fit_a3n_a2n_a1n <- predict.SuperLearner(mu0_sl_a3n_a2n, newdata = df_mu0n)$pred
  df$mu0_fit_a3n_a2n_a1y <- predict.SuperLearner(mu0_sl_a3n_a2n, newdata = df_mu0y)$pred
  
  df$mu0_fit_a3n_a2y_a1n <- predict.SuperLearner(mu0_sl_a3n_a2y, newdata = df_mu0n)$pred
  df$mu0_fit_a3n_a2y_a1y <- predict.SuperLearner(mu0_sl_a3n_a2y, newdata = df_mu0y)$pred
  
  df$mu0_fit_a3y_a2n_a1n <- predict.SuperLearner(mu0_sl_a3y_a2n, newdata = df_mu0n)$pred
  df$mu0_fit_a3y_a2n_a1y <- predict.SuperLearner(mu0_sl_a3y_a2n, newdata = df_mu0y)$pred
  
  df$mu0_fit_a3y_a2y_a1n <- predict.SuperLearner(mu0_sl_a3y_a2y, newdata = df_mu0n)$pred
  df$mu0_fit_a3y_a2y_a1y <- predict.SuperLearner(mu0_sl_a3y_a2y, newdata = df_mu0y)$pred
  
  
  main_list[[k]] <- df[cf_fold[[k]], ]
}

main_df <- reduce(main_list, bind_rows)
  
for (s in 1:S){
  
  a1 <- estimands$a1[[s]]
  a2 <- estimands$a2[[s]]
  a3 <- estimands$a3[[s]]
  
  main_df <- main_df %>%
    mutate(
      
      wt0_deno = a1 * p0_fit + (1 - a1) * (1 - p0_fit),
      wt1_nume = a1 * p1_fit + (1 - a1) * (1 - p1_fit),
      wt1_deno = a2 * p1_fit + (1 - a2) * (1 - p1_fit),
      wt2_nume = a2 * p2_fit + (1 - a2) * (1 - p2_fit),
      wt2_deno = a3 * p2_fit + (1 - a3) * (1 - p2_fit),
      
      !!sym(paste0("w0_", a1, a2, a3)) := as.double(att22==a1) / wt0_deno,
      
      !!sym(paste0("w1_", a1, a2, a3)) := as.double(att22==a2) * wt1_nume/wt1_deno/wt0_deno,
      
      !!sym(paste0("w2_", a1, a2, a3)) := as.double(att22==a3) * wt2_nume/wt2_deno * wt1_nume/wt1_deno/wt0_deno
      
    )
  
  main_df[main_df$att22 == a1, paste0("w0_", a1, a2, a3)] <- trimQ(main_df[main_df$att22 == a1, paste0("w0_", a1, a2, a3)])
  main_df[main_df$att22 == a2, paste0("w1_", a1, a2, a3)] <- trimQ(main_df[main_df$att22 == a2, paste0("w1_", a1, a2, a3)])
  main_df[main_df$att22 == a3, paste0("w2_", a1, a2, a3)] <- trimQ(main_df[main_df$att22 == a3, paste0("w2_", a1, a2, a3)])
  
  main_df <- main_df %>%
    mutate(
      
      !!sym(paste0("mu2fit_", a1, a2, a3)) := a3 * mu2_fit_a3y + (1 - a3) * mu2_fit_a3n,
      
      !!sym(paste0("mu1fit_", a1, a2, a3)) := a3 * a2 * mu1_fit_a3y_a2y + a3 * (1 - a2) * mu1_fit_a3y_a2n +
        (1 - a3) * a2 * mu1_fit_a3n_a2y + (1 - a3) * (1 - a2) * mu1_fit_a3n_a2n,
      
      !!sym(paste0("mu0fit_", a1, a2, a3)) := a3 * a2 * a1 * mu0_fit_a3y_a2y_a1y + a3 * a2 * (1 - a1) * mu0_fit_a3y_a2y_a1n +
        a3 * (1 - a2) * a1 * mu0_fit_a3y_a2n_a1y + a3 * (1 - a2) * (1 - a1) * mu0_fit_a3y_a2n_a1n +
        (1 - a3) * a2 * a1 * mu0_fit_a3n_a2y_a1y + (1 - a3) * a2 * (1 - a1) * mu0_fit_a3n_a2y_a1n +
        (1 - a3) * (1 - a2) * a1 * mu0_fit_a3n_a2n_a1y + (1 - a3) * (1 - a2) * (1 - a1) * mu0_fit_a3n_a2n_a1n,
      
      !!sym(paste0("www_", a1, a2, a3)) := !!sym(paste0("w2_", a1, a2, a3)) * std_cesd_age40,
      !!sym(paste0("iii_", a1, a2, a3)) := !!sym(paste0("mu0fit_", a1, a2, a3)),
      !!sym(paste0("eif_", a1, a2, a3)) := !!sym(paste0("w2_", a1, a2, a3)) *
        (std_cesd_age40 - !!sym(paste0("mu2fit_", a1, a2, a3))) +
        !!sym(paste0("w1_", a1, a2, a3)) * (!!sym(paste0("mu2fit_", a1, a2, a3)) -
                                                                   !!sym(paste0("mu1fit_", a1, a2, a3))) +
        !!sym(paste0("w0_", a1, a2, a3)) * (!!sym(paste0("mu1fit_", a1, a2, a3)) -
                                                                   !!sym(paste0("mu0fit_", a1, a2, a3))) +
        !!sym(paste0("mu0fit_", a1, a2, a3))
    )
}

# set my ggplot theme
mytheme <- theme_minimal(base_size = 18) + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(color = "grey30"))
theme_set(mytheme)

out_df <- main_df %>%
  mutate(eif_type1_ate = eif_111 - eif_000,
         eif_type2_ate = eif_111 - eif_000,
         eif_type3_ate = eif_111 - eif_000,
         eif_type1_pse3 = eif_001 - eif_000,
         eif_type1_pse2 = eif_011 - eif_001,
         eif_type1_pse1 = eif_111 - eif_011,
         eif_type2_pse1 = eif_100 - eif_000,
         eif_type2_pse2 = eif_110 - eif_100,
         eif_type2_pse3 = eif_111 - eif_110,
         eif_type3_pse3 = eif_001 - eif_000,
         eif_type3_pse1 = eif_101 - eif_001,
         eif_type3_pse2 = eif_111 - eif_101) %>%
  summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                            se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator", "type", "estimand", "measure")) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  filter(type == "type1")   %>%
  mutate(estimand = factor(estimand,
                           levels = rev(c("ate", "pse3", "pse2", "pse1")),
                           labels = rev(c(expression(paste("Total Effect (", psi[`111`]-psi[`000`], ")")),
                                          expression(paste("Direct Effect (", psi[`001`]-psi[`000`], ")")),
                                          expression(paste("via Household Income (", psi[`011`]-psi[`001`], ")")),
                                          expression(paste("via Unemployment (", psi[`111`]-psi[`011`], ")")))))) 

ggplot(out_df, aes(x = estimand, y = est, shape = estimator)) +
  geom_pointrange(aes(ymin = est - 1.96 * se,  ymax = est + 1.96 * se),
                  position = position_dodge(width = - 0.5), size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape("", labels = parse_format()) +
  scale_color_discrete("", labels = parse_format()) +
  scale_x_discrete("", labels = parse_format()) +
  scale_y_continuous("Effects of College Attendance on Depression") +
  coord_flip()

table6_4np <-  out_df %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se) %>% 
  mutate_at(c("est", "se", "lower", "upper"), ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")")) %>% 
  mutate(out = paste(est, intv)) %>% 
  dplyr::select(-lower, -upper, -intv) 

write_csv(table6_4np, file = "table6-4np.csv")

save.image(file = "table6-4np.RData")
