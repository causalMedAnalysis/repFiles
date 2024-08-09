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

y0_form <- as.formula(paste(y, " ~ ", paste(c(x, a), collapse= "+")))
y1_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1), collapse= "+")))

m_form <- as.formula(paste(m1, " ~ ", paste(c(x, a), collapse= "+")))

##########################################################
# Main analyses
##########################################################

estimands <- expand.grid(c(0, 1), c(0, 1)) %>%
  `colnames<-`(c("a1", "a2"))

S <- nrow(estimands)

K <- 5

df_p0 <- model.matrix(a0_form, data = df)[, -1] %>% as_tibble()
df_mu1 <- model.matrix(y1_form, data = df)[, -1] %>% as_tibble()

df_mu1_a2n <- model.matrix(y1_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_mu1_a2y <- model.matrix(y1_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()

df_mu1_a2n_m0 <- model.matrix(y1_form, data = mutate(df, att22 = 0, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
df_mu1_a2y_m0 <- model.matrix(y1_form, data = mutate(df, att22 = 1, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()

df_mu1_a2n_m1 <- model.matrix(y1_form, data = mutate(df, att22 = 0, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()
df_mu1_a2y_m1 <- model.matrix(y1_form, data = mutate(df, att22 = 1, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()

df_M <- model.matrix(m_form, data = df)[, -1] %>% as_tibble()
df_M_a1n <- model.matrix(m_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_M_a1y <- model.matrix(m_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()

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

  aux_mu1 <- df_mu1[-cf_fold[[k]], ]
  main_mu1 <- df_mu1[cf_fold[[k]], ]
  
  aux_M <- df_M[-cf_fold[[k]], ]
  main_M <- df_M[cf_fold[[k]], ]
  
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
  
  m_sl <- SuperLearner(
    Y          = aux$ever_unemp_age3539,
    X          = aux_M,
    newX       = df_M,
    family     = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  df <- df %>% mutate(
    
    p0_fit = p0_sl$SL.predict,
    
    M1_fit_a1n = predict.SuperLearner(m_sl, newdata = df_M_a1n)$pred,
    M0_fit_a1n = 1 - M1_fit_a1n,
    
    M1_fit_a1y = predict.SuperLearner(m_sl, newdata = df_M_a1y)$pred,
    M0_fit_a1y = 1 - M1_fit_a1y,
    
    fM_fit_a1n = ifelse(ever_unemp_age3539 == 0, M0_fit_a1n, M1_fit_a1n),
    fM_fit_a1y = ifelse(ever_unemp_age3539 == 0, M0_fit_a1y, M1_fit_a1y),
    
  )
  
  #################################################
  # Outcome models
  #################################################
  
  mu1_sl <- SuperLearner(
    Y          = aux$std_cesd_age40,
    X          = aux_mu1,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu1_fit_a2n <- predict.SuperLearner(mu1_sl, newdata = df_mu1_a2n)$pred
  df$mu1_fit_a2y <- predict.SuperLearner(mu1_sl, newdata = df_mu1_a2y)$pred
  
  df$mu1_fit_a2n_m0 <- predict.SuperLearner(mu1_sl, newdata = df_mu1_a2n_m0)$pred
  df$mu1_fit_a2y_m0 <- predict.SuperLearner(mu1_sl, newdata = df_mu1_a2y_m0)$pred
  df$mu1_fit_a2n_m1 <- predict.SuperLearner(mu1_sl, newdata = df_mu1_a2n_m1)$pred
  df$mu1_fit_a2y_m1 <- predict.SuperLearner(mu1_sl, newdata = df_mu1_a2y_m1)$pred
  
  df <- df %>%
    mutate(
      
      mu0_fit_a2n_a1n = mu1_fit_a2n_m0 * M0_fit_a1n + mu1_fit_a2n_m1 * M1_fit_a1n,
      mu0_fit_a2n_a1y = mu1_fit_a2n_m0 * M0_fit_a1y + mu1_fit_a2n_m1 * M1_fit_a1y,
      
      mu0_fit_a2y_a1n = mu1_fit_a2y_m0 * M0_fit_a1n + mu1_fit_a2y_m1 * M1_fit_a1n,
      mu0_fit_a2y_a1y = mu1_fit_a2y_m0 * M0_fit_a1y + mu1_fit_a2y_m1 * M1_fit_a1y
    )
  
  main_list[[k]] <- df[cf_fold[[k]], ]
}

main_df <- reduce(main_list, bind_rows)

for (s in 1:S){
  
  a1 <- estimands$a1[[s]]
  a2 <- estimands$a2[[s]]
  
  main_df <- main_df %>%
    mutate(
      
      wt0_deno = a2 * p0_fit + (1 - a2) * (1 - p0_fit),
      wt1_nume = a1 * fM_fit_a1y + (1 - a1) * fM_fit_a1n,
      wt1_deno = a2 * fM_fit_a1y + (1 - a2) * fM_fit_a1n,
      
      !!sym(paste0("w0_", a1, a2)) := as.double(att22==a1) / (a1 * p0_fit + (1 - a1) * (1 - p0_fit)),
      
      !!sym(paste0("w1_", a1, a2)) := as.double(att22==a2) * wt1_nume/wt1_deno/wt0_deno,
      
    )
  
  main_df[main_df$att22 == a1, paste0("w0_", a1, a2)] <- trimQ(main_df[main_df$att22 == a1, paste0("w0_", a1, a2)])
  main_df[main_df$att22 == a2, paste0("w1_", a1, a2)] <- trimQ(main_df[main_df$att22 == a2, paste0("w1_", a1, a2)])
  
  main_df <- main_df %>%
    mutate(
      
      !!sym(paste0("mu1fit_", a1, a2)) := a2 * mu1_fit_a2y + (1 - a2) * mu1_fit_a2n,
      
      !!sym(paste0("mu0fit_", a1, a2)) := a2 * a1 * mu0_fit_a2y_a1y + a2 * (1 - a1) * mu0_fit_a2y_a1n +
        (1 - a2) * a1 * mu0_fit_a2n_a1y + (1 - a2) * (1 - a1) * mu0_fit_a2n_a1n,
      
      !!sym(paste0("www_", a1, a2)) := !!sym(paste0("w1_", a1, a2)) * std_cesd_age40,
      !!sym(paste0("iii_", a1, a2)) := !!sym(paste0("mu0fit_", a1, a2)),
      !!sym(paste0("eif_", a1, a2)) := !!sym(paste0("w1_", a1, a2)) *
        (std_cesd_age40 - !!sym(paste0("mu1fit_", a1, a2))) +
        !!sym(paste0("w0_", a1, a2)) * (!!sym(paste0("mu1fit_", a1, a2)) -
                                          !!sym(paste0("mu0fit_", a1, a2))) +
        !!sym(paste0("mu0fit_", a1, a2))
    )
}

out_df <- main_df %>%
  mutate(eif_type1_ate = eif_11 - eif_00,
         eif_type2_ate = eif_11 - eif_00,
         eif_type1_pse2 = eif_01 - eif_00,
         eif_type1_pse1 = eif_11 - eif_01,
         eif_type2_pse1 = eif_10 - eif_00,
         eif_type2_pse2 = eif_11 - eif_10) %>%
  summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                            se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator", "type", "estimand", "measure")) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  filter(type == "type1")   %>%
  mutate(estimand = factor(estimand,
                           levels = rev(c("ate", "pse2", "pse1")),
                           labels = rev(c(expression(paste("Total Effect (", psi[`11`]-psi[`00`], ")")),
                                          expression(paste("Direct Effect (", psi[`01`]-psi[`00`], ")")),
                                          expression(paste("via Unemployment (", psi[`11`]-psi[`01`], ")")))))) 


# set my ggplot theme
mytheme <- theme_minimal(base_size = 18) + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(color = "grey30"))
theme_set(mytheme)

ggplot(out_df, aes(x = estimand, y = est, shape = estimator)) +
  geom_pointrange(aes(ymin = est - 1.96 * se,  ymax = est + 1.96 * se),
                  position = position_dodge(width = - 0.5), size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape("", labels = parse_format()) +
  scale_color_discrete("", labels = parse_format()) +
  scale_x_discrete("", labels = parse_format()) +
  scale_y_continuous("Effects of College Attendance on Depression") +
  coord_flip()

table6_2np_a <-  out_df %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se) %>% 
  mutate_at(c("est", "se", "lower", "upper"), ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")")) %>% 
  mutate(out = paste(est, intv)) %>% 
  dplyr::select(-lower, -upper, -intv) 

write_csv(table6_2np_a, file = "table6-2np_a.csv")

save.image(file = "table6-2np_a.RData")
