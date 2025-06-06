rm(list = ls())

library(survey)
library(gbm)
library(ranger)
library(glmnet)
library(rsample)
library(caret)
library(rlang)
library(tidyverse)
library(dplyr)
library(Hmisc)
library(SuperLearner)
library(scales)
library(haven)

source("https://raw.githubusercontent.com/causalMedAnalysis/causalMedR/refs/heads/main/utils.R")

set.seed(02138)

##input data
nlsy_raw <- as.data.frame(read_stata(
  file = "https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta"
  )
)

a <- "att22"
z <- "ever_unemp_age3539"
m <- "log_faminc_adj_age3539"
y <- "std_cesd_age40"
x <- c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")

nlsy <- nlsy_raw[complete.cases(nlsy_raw[, c(a, "cesd_age40", z, m, x)]),] %>% 
  mutate(std_cesd_age40 = as.numeric(scale(cesd_age40)))

df <- nlsy
n <- nrow(df)

##########################################################
# Formulas for a, z, y models
##########################################################

# b(x, a, z,m) = E[Y|x, a, z, m]
b_form <- as.formula(paste(y, " ~ ", paste(c(x, a, z, m), collapse= "+")))

# g(a|x) = Pr[A = a|x]
g_form <- as.formula(paste(a, " ~ ", paste(c(x), collapse= "+")))

# h(a|x, m) = Pr[A = a|x, m]
h_form <- as.formula(paste(a, " ~ ", paste(c(x, m), collapse= "+")))

# q(z|x, a) = Pr[Z = z|x, a]
q_form <- as.formula(paste(z, " ~ ", paste(c(x, a), collapse= "+")))

# r(z|x, a, m) = Pr[Z = z|x, a, m]
r_form <- as.formula(paste(z, " ~ ", paste(c(x, a, m), collapse= "+")))

# c(x, a, z, m) = g(a|w)/g(a1|w) * q(z|x,a)/r(z|x,a,m) * h(a1|x,m)/h(a|x,m)

# u(x, a, z) = E[b(X,A,Z,M)c(X,A,Z,M)|x, a, z]
u_form <- as.formula(paste(y, " ~ ", paste(c(x, a, z), collapse= "+")))

##########################################################
# Main analyses
##########################################################

estimands <- expand.grid(c(0, 1), c(0, 1)) %>%
  `colnames<-`(c("a1", "a2"))

S <- nrow(estimands)

#################################################
# Design matrices for different models
#################################################

df_b <- model.matrix(b_form, data = df)[, -1] %>% as_tibble()
df_g <- model.matrix(g_form, data = df)[, -1] %>% as_tibble()
df_h <- model.matrix(h_form, data = df)[, -1] %>% as_tibble()
df_q <- model.matrix(q_form, data = df)[, -1] %>% as_tibble()
df_r <- model.matrix(r_form, data = df)[, -1] %>% as_tibble()
df_u <- model.matrix(u_form, data = df)[, -1] %>% as_tibble()

df_b0 <- model.matrix(b_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_q0 <- model.matrix(q_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_r0 <- model.matrix(r_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_u0 <- model.matrix(u_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()

df_b1 <- model.matrix(b_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()
df_q1 <- model.matrix(q_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()
df_r1 <- model.matrix(r_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()
df_u1 <- model.matrix(u_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()

df_b00 <- model.matrix(b_form, data = mutate(df, att22 = 0, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
df_b01 <- model.matrix(b_form, data = mutate(df, att22 = 0, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()
df_b10 <- model.matrix(b_form, data = mutate(df, att22 = 1, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
df_b11 <- model.matrix(b_form, data = mutate(df, att22 = 1, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()

df_u00 <- model.matrix(u_form, data = mutate(df, att22 = 0, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
df_u01 <- model.matrix(u_form, data = mutate(df, att22 = 0, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()
df_u10 <- model.matrix(u_form, data = mutate(df, att22 = 1, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
df_u11 <- model.matrix(u_form, data = mutate(df, att22 = 1, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()

#################################################
# Outcome model
#################################################

b_glm <- lm(b_form, data = df)

#################################################
# Treatment Models
#################################################

g_glm <- glm(g_form, family = binomial("logit"), data = df)
h_glm <- glm(h_form, family = binomial("logit"), data = df)

#################################################
# Z Models
#################################################

q_glm <- glm(q_form, family = binomial("logit"), data = df)
r_glm <- glm(r_form, family = binomial("logit"), data = df)

df <- df %>% mutate(
  
  g0_fit = 1 - g_glm$fitted.values,
  g1_fit = g_glm$fitted.values,
  gA_fit = att22 * g1_fit + (1 - att22) * g0_fit,
  
  h0_fit = 1 - h_glm$fitted.values,
  h1_fit = h_glm$fitted.values,
  hA_fit = att22 * h1_fit + (1 - att22) * h0_fit,
  
  q01_fit = predict(q_glm, newdata = df_q0, type = "response"),
  q11_fit = predict(q_glm, newdata = df_q1, type = "response"),
  q00_fit = 1 - q01_fit,
  q10_fit = 1 - q11_fit,
  
  r01_fit = predict(r_glm, newdata = df_r0, type = "response"),
  r11_fit = predict(r_glm, newdata = df_r1, type = "response"),
  r00_fit = 1 - r01_fit,
  r10_fit = 1 - r11_fit,
  
  q0Z_fit = ifelse(ever_unemp_age3539 == 1, q01_fit, q00_fit),
  q1Z_fit = ifelse(ever_unemp_age3539 == 1, q11_fit, q10_fit),
  qAZ_fit = ifelse(att22 == 1, q1Z_fit, q0Z_fit),
  
  r0Z_fit = ifelse(ever_unemp_age3539 == 1, r01_fit, r00_fit),
  r1Z_fit = ifelse(ever_unemp_age3539 == 1, r11_fit, r10_fit),
  rAZ_fit = ifelse(att22 == 1, r1Z_fit, r0Z_fit),
  
  b0_fit = predict(b_glm, newdata = df_b0, type = "response"),
  b1_fit = predict(b_glm, newdata = df_b1, type = "response"),
  bA_fit = ifelse(att22 == 1, b1_fit, b0_fit),
  
  b00_fit = predict(b_glm, newdata = df_b00, type = "response"),
  b01_fit = predict(b_glm, newdata = df_b01, type = "response"),
  b10_fit = predict(b_glm, newdata = df_b10, type = "response"),
  b11_fit = predict(b_glm, newdata = df_b11, type = "response"),
  
  c_a1n_fit = gA_fit/g0_fit * qAZ_fit/rAZ_fit * h0_fit/hA_fit,
  c_a1y_fit = gA_fit/g1_fit * qAZ_fit/rAZ_fit * h1_fit/hA_fit,
  
  c00_fit = q0Z_fit/r0Z_fit,
  c01_fit = g1_fit/g0_fit * q1Z_fit/r1Z_fit * h0_fit/h1_fit,
  c10_fit = g0_fit/g1_fit * q0Z_fit/r0Z_fit * h1_fit/h0_fit,
  c11_fit = q1Z_fit/r1Z_fit,
  
  dep_umod_a1n = bA_fit * c_a1n_fit,
  dep_umod_a1y = bA_fit * c_a1y_fit,
  
  dep_vmod_a2n = b00_fit * q00_fit + b01_fit * q01_fit,
  dep_vmod_a2y = b10_fit * q10_fit + b11_fit * q11_fit,
)

#################################################
# U and V models
#################################################

u_a1n_form <- as.formula(paste("dep_umod_a1n", " ~ ", paste(c(x, a, z), collapse= "+")))
u_a1y_form <- as.formula(paste("dep_umod_a1y", " ~ ", paste(c(x, a, z), collapse= "+")))

v_a2n_form <- as.formula(paste("dep_vmod_a2n", " ~ ", paste(c(x, a), collapse= "+")))
v_a2y_form <- as.formula(paste("dep_vmod_a2y", " ~ ", paste(c(x, a), collapse= "+")))

u_a1n_glm <- lm(u_a1n_form, data = df)
u_a1y_glm <- lm(u_a1y_form, data = df)
v_a2n_glm <- lm(v_a2n_form, data = df)
v_a2y_glm <- lm(v_a2y_form, data = df)

df <- df %>% mutate(
  
  u0_a1n_fit =  predict(u_a1n_glm, newdata = df_u0),
  u1_a1n_fit =  predict(u_a1n_glm, newdata = df_u1),
  
  u0_a1y_fit =  predict(u_a1y_glm, newdata = df_u0),
  u1_a1y_fit =  predict(u_a1y_glm, newdata = df_u1),
  
  v0_a2n_fit =  predict(v_a2n_glm, newdata = df_q0),
  v1_a2n_fit =  predict(v_a2n_glm, newdata = df_q1),
  
  v0_a2y_fit =  predict(v_a2y_glm, newdata = df_q0),
  v1_a2y_fit =  predict(v_a2y_glm, newdata = df_q1),
  
  u00_a1n_fit =  predict(u_a1n_glm, newdata = df_u00),
  u01_a1n_fit =  predict(u_a1n_glm, newdata = df_u01),
  u10_a1n_fit =  predict(u_a1n_glm, newdata = df_u10),
  u11_a1n_fit =  predict(u_a1n_glm, newdata = df_u11),
  
  u00_a1y_fit =  predict(u_a1y_glm, newdata = df_u00),
  u01_a1y_fit =  predict(u_a1y_glm, newdata = df_u01),
  u10_a1y_fit =  predict(u_a1y_glm, newdata = df_u10),
  u11_a1y_fit =  predict(u_a1y_glm, newdata = df_u11),
  
  uqOz_nn_fit = u00_a1n_fit * q00_fit +  u01_a1n_fit * q01_fit,
  uqOz_ny_fit = u10_a1n_fit * q10_fit +  u11_a1n_fit * q11_fit,
  uqOz_yn_fit = u00_a1y_fit * q00_fit +  u01_a1y_fit * q01_fit,
  uqOz_yy_fit = u10_a1y_fit * q10_fit +  u11_a1y_fit * q11_fit,
  
  bqOz_a2n_fit = b00_fit * q00_fit +  b01_fit * q01_fit,
  bqOz_a2y_fit = b10_fit * q10_fit +  b11_fit * q11_fit,
)

for (s in 1:S){
  
  a1 <- estimands$a1[[s]]
  a2 <- estimands$a2[[s]]
  
  df <- df %>%
    mutate(
      
      g_a2_fit = a2 * g1_fit + (1 - a2) * g0_fit,
      g_a1_fit = a1 * g1_fit + (1 - a1) * g0_fit,
      
      c_fit = a1 * a2 * c11_fit +
        a1 * (1 - a2) * c10_fit +
        (1 - a1) * a2 * c01_fit +
        (1 - a1) * (1 - a2) * c00_fit,
      
      b_fit = a2 * b1_fit + (1 - a2) * b0_fit,
      
      u_fit = a1 * a2 * u1_a1y_fit +
        a1 * (1 - a2) * u0_a1y_fit +
        (1 - a1) * a2 * u1_a1n_fit +
        (1 - a1) * (1 - a2) * u0_a1n_fit,
      
      v_fit = a1 * a2 * v1_a2y_fit +
        a1 * (1 - a2) * v1_a2n_fit +
        (1 - a1) * a2 * v0_a2y_fit +
        (1 - a1) * (1 - a2) * v0_a2n_fit,
      
      uqOz_fit = a1 * a2 * uqOz_yy_fit +
        a1 * (1 - a2) * uqOz_yn_fit +
        (1 - a1) * a2 * uqOz_ny_fit +
        (1 - a1) * (1 - a2) * uqOz_nn_fit,
      
      bqOz_fit = a2 * bqOz_a2y_fit + (1 - a2) * bqOz_a2n_fit,
      
      !!sym(paste0("w2_", a1, a2)) := as.double(att22==a2) * c_fit/g_a2_fit,
      
      !!sym(paste0("w1_", a1, a2)) := as.double(att22==a2) /g_a2_fit,
      
      !!sym(paste0("w0_", a1, a2)) := as.double(att22==a1) /g_a1_fit     
    )
  
  df[df$att22 == a2, paste0("w2_", a1, a2)] <- trimQ(df[df$att22 == a2, paste0("w2_", a1, a2)])
  df[df$att22 == a2, paste0("w1_", a1, a2)] <- trimQ(df[df$att22 == a2, paste0("w1_", a1, a2)])
  df[df$att22 == a1, paste0("w0_", a1, a2)] <- trimQ(df[df$att22 == a1, paste0("w0_", a1, a2)])
  
  df <- df %>%
    mutate(
      
      !!sym(paste0("eif_", a1, a2)) := !!sym(paste0("w2_", a1, a2)) * (std_cesd_age40 - b_fit) +
        !!sym(paste0("w1_", a1, a2)) * (u_fit - uqOz_fit) +
        !!sym(paste0("w0_", a1, a2)) * (bqOz_fit - v_fit) +
        v_fit
    )
}

out_df <- df %>%
  mutate(eif_type1_ate = eif_11 - eif_00,
         eif_type2_ate = eif_11 - eif_00,
         eif_type1_iie = eif_10 - eif_00,
         eif_type1_ide = eif_11 - eif_10,
         eif_type2_iie = eif_11 - eif_01,
         eif_type2_ide = eif_01 - eif_00) %>%
  summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                            se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator", "type", "estimand", "measure")) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  filter(type == "type1")   %>%
  mutate(estimand = factor(estimand,
                           levels = rev(c("ate", "ide", "iie")),
                           labels = rev(c(expression(paste("Total Effect (", psi[`11`]-psi[`00`], ")")),
                                          expression(paste("Direct Effect (", psi[`01`]-psi[`00`], ")")),
                                          expression(paste("via Household Income (", psi[`11`]-psi[`01`], ")")))))) 

#################################################
# Nonparametric Bootstrap
#################################################

B <- 2000

boots <- matrix(NA, nrow = B, ncol = 3)

for (b in 1:B){
  
  if (b %% 10 ==0){
    cat("bootstrap sample ", b, "\n")
  }
  
  dfi <- nlsy %>% sample_frac(replace = TRUE)
  
  #################################################
  # Design matrices for different models
  #################################################
  
  dfi_b <- model.matrix(b_form, data = dfi)[, -1] %>% as_tibble()
  dfi_g <- model.matrix(g_form, data = dfi)[, -1] %>% as_tibble()
  dfi_h <- model.matrix(h_form, data = dfi)[, -1] %>% as_tibble()
  dfi_q <- model.matrix(q_form, data = dfi)[, -1] %>% as_tibble()
  dfi_r <- model.matrix(r_form, data = dfi)[, -1] %>% as_tibble()
  dfi_u <- model.matrix(u_form, data = dfi)[, -1] %>% as_tibble()
  
  dfi_b0 <- model.matrix(b_form, data = mutate(dfi, att22 = 0))[, -1] %>% as_tibble()
  dfi_q0 <- model.matrix(q_form, data = mutate(dfi, att22 = 0))[, -1] %>% as_tibble()
  dfi_r0 <- model.matrix(r_form, data = mutate(dfi, att22 = 0))[, -1] %>% as_tibble()
  dfi_u0 <- model.matrix(u_form, data = mutate(dfi, att22 = 0))[, -1] %>% as_tibble()
  
  dfi_b1 <- model.matrix(b_form, data = mutate(dfi, att22 = 1))[, -1] %>% as_tibble()
  dfi_q1 <- model.matrix(q_form, data = mutate(dfi, att22 = 1))[, -1] %>% as_tibble()
  dfi_r1 <- model.matrix(r_form, data = mutate(dfi, att22 = 1))[, -1] %>% as_tibble()
  dfi_u1 <- model.matrix(u_form, data = mutate(dfi, att22 = 1))[, -1] %>% as_tibble()
  
  dfi_b00 <- model.matrix(b_form, data = mutate(dfi, att22 = 0, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
  dfi_b01 <- model.matrix(b_form, data = mutate(dfi, att22 = 0, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()
  dfi_b10 <- model.matrix(b_form, data = mutate(dfi, att22 = 1, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
  dfi_b11 <- model.matrix(b_form, data = mutate(dfi, att22 = 1, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()
  
  dfi_u00 <- model.matrix(u_form, data = mutate(dfi, att22 = 0, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
  dfi_u01 <- model.matrix(u_form, data = mutate(dfi, att22 = 0, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()
  dfi_u10 <- model.matrix(u_form, data = mutate(dfi, att22 = 1, ever_unemp_age3539 = 0))[, -1] %>% as_tibble()
  dfi_u11 <- model.matrix(u_form, data = mutate(dfi, att22 = 1, ever_unemp_age3539 = 1))[, -1] %>% as_tibble()
  
    #################################################
  # Outcome model
  #################################################
  
  b_glm <- lm(b_form, data = dfi)
  
  #################################################
  # Treatment Models
  #################################################
  
  g_glm <- glm(g_form, family = binomial("logit"), data = dfi)
  h_glm <- glm(h_form, family = binomial("logit"), data = dfi)
  
  #################################################
  # Z Models
  #################################################
  
  q_glm <- glm(q_form, family = binomial("logit"), data = dfi)
  r_glm <- glm(r_form, family = binomial("logit"), data = dfi)
  
  dfi <- dfi %>% mutate(
    
    g0_fit = 1 - g_glm$fitted.values,
    g1_fit = g_glm$fitted.values,
    gA_fit = att22 * g1_fit + (1 - att22) * g0_fit,
    
    h0_fit = 1 - h_glm$fitted.values,
    h1_fit = h_glm$fitted.values,
    hA_fit = att22 * h1_fit + (1 - att22) * h0_fit,
    
    q01_fit = predict(q_glm, newdata = dfi_q0, type = "response"),
    q11_fit = predict(q_glm, newdata = dfi_q1, type = "response"),
    q00_fit = 1 - q01_fit,
    q10_fit = 1 - q11_fit,
    
    r01_fit = predict(r_glm, newdata = dfi_r0, type = "response"),
    r11_fit = predict(r_glm, newdata = dfi_r1, type = "response"),
    r00_fit = 1 - r01_fit,
    r10_fit = 1 - r11_fit,
    
    q0Z_fit = ifelse(ever_unemp_age3539 == 1, q01_fit, q00_fit),
    q1Z_fit = ifelse(ever_unemp_age3539 == 1, q11_fit, q10_fit),
    qAZ_fit = ifelse(att22 == 1, q1Z_fit, q0Z_fit),
    
    r0Z_fit = ifelse(ever_unemp_age3539 == 1, r01_fit, r00_fit),
    r1Z_fit = ifelse(ever_unemp_age3539 == 1, r11_fit, r10_fit),
    rAZ_fit = ifelse(att22 == 1, r1Z_fit, r0Z_fit),
    
    b0_fit = predict(b_glm, newdata = dfi_b0, type = "response"),
    b1_fit = predict(b_glm, newdata = dfi_b1, type = "response"),
    bA_fit = ifelse(att22 == 1, b1_fit, b0_fit),
    
    b00_fit = predict(b_glm, newdata = dfi_b00, type = "response"),
    b01_fit = predict(b_glm, newdata = dfi_b01, type = "response"),
    b10_fit = predict(b_glm, newdata = dfi_b10, type = "response"),
    b11_fit = predict(b_glm, newdata = dfi_b11, type = "response"),
    
    c_a1n_fit = gA_fit/g0_fit * qAZ_fit/rAZ_fit * h0_fit/hA_fit,
    c_a1y_fit = gA_fit/g1_fit * qAZ_fit/rAZ_fit * h1_fit/hA_fit,
    
    c00_fit = q0Z_fit/r0Z_fit,
    c01_fit = g1_fit/g0_fit * q1Z_fit/r1Z_fit * h0_fit/h1_fit,
    c10_fit = g0_fit/g1_fit * q0Z_fit/r0Z_fit * h1_fit/h0_fit,
    c11_fit = q1Z_fit/r1Z_fit,
    
    dep_umod_a1n = bA_fit * c_a1n_fit,
    dep_umod_a1y = bA_fit * c_a1y_fit,
    
    dep_vmod_a2n = b00_fit * q00_fit + b01_fit * q01_fit,
    dep_vmod_a2y = b10_fit * q10_fit + b11_fit * q11_fit,
  )
  
  #################################################
  # U and V models
  #################################################
  
  u_a1n_form <- as.formula(paste("dep_umod_a1n", " ~ ", paste(c(x, a, z), collapse= "+")))
  u_a1y_form <- as.formula(paste("dep_umod_a1y", " ~ ", paste(c(x, a, z), collapse= "+")))
  
  v_a2n_form <- as.formula(paste("dep_vmod_a2n", " ~ ", paste(c(x, a), collapse= "+")))
  v_a2y_form <- as.formula(paste("dep_vmod_a2y", " ~ ", paste(c(x, a), collapse= "+")))
  
  u_a1n_glm <- lm(u_a1n_form, data = dfi)
  u_a1y_glm <- lm(u_a1y_form, data = dfi)
  v_a2n_glm <- lm(v_a2n_form, data = dfi)
  v_a2y_glm <- lm(v_a2y_form, data = dfi)
  
  dfi <- dfi %>% mutate(
    
    u0_a1n_fit =  predict(u_a1n_glm, newdata = dfi_u0),
    u1_a1n_fit =  predict(u_a1n_glm, newdata = dfi_u1),
    
    u0_a1y_fit =  predict(u_a1y_glm, newdata = dfi_u0),
    u1_a1y_fit =  predict(u_a1y_glm, newdata = dfi_u1),
    
    v0_a2n_fit =  predict(v_a2n_glm, newdata = dfi_q0),
    v1_a2n_fit =  predict(v_a2n_glm, newdata = dfi_q1),
    
    v0_a2y_fit =  predict(v_a2y_glm, newdata = dfi_q0),
    v1_a2y_fit =  predict(v_a2y_glm, newdata = dfi_q1),
    
    u00_a1n_fit =  predict(u_a1n_glm, newdata = dfi_u00),
    u01_a1n_fit =  predict(u_a1n_glm, newdata = dfi_u01),
    u10_a1n_fit =  predict(u_a1n_glm, newdata = dfi_u10),
    u11_a1n_fit =  predict(u_a1n_glm, newdata = dfi_u11),
    
    u00_a1y_fit =  predict(u_a1y_glm, newdata = dfi_u00),
    u01_a1y_fit =  predict(u_a1y_glm, newdata = dfi_u01),
    u10_a1y_fit =  predict(u_a1y_glm, newdata = dfi_u10),
    u11_a1y_fit =  predict(u_a1y_glm, newdata = dfi_u11),
    
    uqOz_nn_fit = u00_a1n_fit * q00_fit +  u01_a1n_fit * q01_fit,
    uqOz_ny_fit = u10_a1n_fit * q10_fit +  u11_a1n_fit * q11_fit,
    uqOz_yn_fit = u00_a1y_fit * q00_fit +  u01_a1y_fit * q01_fit,
    uqOz_yy_fit = u10_a1y_fit * q10_fit +  u11_a1y_fit * q11_fit,
    
    bqOz_a2n_fit = b00_fit * q00_fit +  b01_fit * q01_fit,
    bqOz_a2y_fit = b10_fit * q10_fit +  b11_fit * q11_fit,
  )
  
  for (s in 1:S){
    
    a1 <- estimands$a1[[s]]
    a2 <- estimands$a2[[s]]
    
    dfi <- dfi %>%
      mutate(
        
        g_a2_fit = a2 * g1_fit + (1 - a2) * g0_fit,
        g_a1_fit = a1 * g1_fit + (1 - a1) * g0_fit,
        
        c_fit = a1 * a2 * c11_fit +
          a1 * (1 - a2) * c10_fit +
          (1 - a1) * a2 * c01_fit +
          (1 - a1) * (1 - a2) * c00_fit,
        
        b_fit = a2 * b1_fit + (1 - a2) * b0_fit,
        
        u_fit = a1 * a2 * u1_a1y_fit +
          a1 * (1 - a2) * u0_a1y_fit +
          (1 - a1) * a2 * u1_a1n_fit +
          (1 - a1) * (1 - a2) * u0_a1n_fit,
        
        v_fit = a1 * a2 * v1_a2y_fit +
          a1 * (1 - a2) * v1_a2n_fit +
          (1 - a1) * a2 * v0_a2y_fit +
          (1 - a1) * (1 - a2) * v0_a2n_fit,
        
        uqOz_fit = a1 * a2 * uqOz_yy_fit +
          a1 * (1 - a2) * uqOz_yn_fit +
          (1 - a1) * a2 * uqOz_ny_fit +
          (1 - a1) * (1 - a2) * uqOz_nn_fit,
        
        bqOz_fit = a2 * bqOz_a2y_fit + (1 - a2) * bqOz_a2n_fit,
        
        !!sym(paste0("w2_", a1, a2)) := as.double(att22==a2) * c_fit/g_a2_fit,
        
        !!sym(paste0("w1_", a1, a2)) := as.double(att22==a2) /g_a2_fit,
        
        !!sym(paste0("w0_", a1, a2)) := as.double(att22==a1) /g_a1_fit     
      )
    
    dfi[dfi$att22 == a2, paste0("w2_", a1, a2)] <- trimQ(dfi[dfi$att22 == a2, paste0("w2_", a1, a2)])
    dfi[dfi$att22 == a2, paste0("w1_", a1, a2)] <- trimQ(dfi[dfi$att22 == a2, paste0("w1_", a1, a2)])
    dfi[dfi$att22 == a1, paste0("w0_", a1, a2)] <- trimQ(dfi[dfi$att22 == a1, paste0("w0_", a1, a2)])
    
    dfi <- dfi %>%
      mutate(
        
        !!sym(paste0("eif_", a1, a2)) := !!sym(paste0("w2_", a1, a2)) * (std_cesd_age40 - b_fit) +
          !!sym(paste0("w1_", a1, a2)) * (u_fit - uqOz_fit) +
          !!sym(paste0("w0_", a1, a2)) * (bqOz_fit - v_fit) +
          v_fit
      )
  }
  
  out_dfi <- dfi %>%
    mutate(eif_type1_ate = eif_11 - eif_00,
           eif_type2_ate = eif_11 - eif_00,
           eif_type1_iie = eif_10 - eif_00,
           eif_type1_ide = eif_11 - eif_10,
           eif_type2_iie = eif_11 - eif_01,
           eif_type2_ide = eif_01 - eif_00) %>%
    summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                              se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
    pivot_longer(everything()) %>%
    separate(name, into = c("estimator", "type", "estimand", "measure")) %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    filter(type == "type1")   %>%
    mutate(estimand = factor(estimand,
                             levels = rev(c("ate", "ide", "iie")),
                             labels = rev(c(expression(paste("Total Effect (", psi[`11`]-psi[`00`], ")")),
                                            expression(paste("Direct Effect (", psi[`01`]-psi[`00`], ")")),
                                            expression(paste("via Household Income (", psi[`11`]-psi[`01`], ")")))))) 
  
  boots[b, ] <- out_dfi$est
  
}

out_df <- out_df %>% 
  mutate(
    se_naive = se,
    se = apply(boots, 2, sd),
    lower = apply(boots, 2, quantile, 0.025),
    upper = apply(boots, 2, quantile, 0.975)
  )

table6_3par <-  out_df %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se) %>% 
  mutate_at(c("est", "se", "lower", "upper"), ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")")) %>% 
  mutate(out = paste(est, intv)) %>% 
  dplyr::select(-lower, -upper, -intv) %>% 
  arrange(desc(estimand))

print(table6_3par)
