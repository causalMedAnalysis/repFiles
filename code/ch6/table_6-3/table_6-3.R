#----------Preliminaries----------#
rm(list = ls())
chapter <- "ch6"
title   <- "table_6-3"

# Specify the root directory:
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"

# Define subdirectories for logs and figures:
dir_log  <- file.path(dir_root, "code", chapter, "_LOGS")
log_path <- file.path(dir_log,  paste0(title, "_log.txt"))

# Ensure all necessary directories exist under your root folder
# If not, the following function will create folders for you.

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
#  Causal Mediation Analysis Replication Files
#
#  GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main
#
#  Script:      .../code/ch6/table_6-3.R
#
#  Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta
#
#  Outputs:     .../code/ch6/_LOGS/table_6-3_log.txt
#
#  Description: Replicates Chapter 6, Table 6.3: Multiply Robust Estimates for the
#               Interventional Effects of College Attendance on CES-D scores, as
#               Mediated by Household Income, from the NLSY
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#

packages <- c(
  "tidyverse", 
  "rlang", 
  "Hmisc",
  "survey", 
  "gbm", 
  "ranger", 
  "glmnet",
  "rsample", 
  "caret", 
  "SuperLearner",
  "scales", 
  "haven", 
  "devtools", 
  "dplyr"
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

#----------------------#
#    SPECIFICATIONS    #
#----------------------#

# outcome
y <- "std_cesd_age40"

# exposure
a <- "att22"

# exposure-induced confounder
z <- "ever_unemp_age3539"

# mediator 
m <- "log_faminc_adj_age3539"

# baseline confounders
x <- c("female",
       "black",
       "hispan",
       "paredu",
       "parprof",
       "parinc_prank",
       "famsize",
       "afqt3"
)

#set seed
set.seed(02138)

#-----------------------------#
#        PREPARE DATA         #
#-----------------------------#

nlsy_raw <- as.data.frame(
  read_stata("https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/NLSY79/nlsy79BK_ed2.dta")
)

nlsy <- nlsy_raw[complete.cases(nlsy_raw[, c(a, "cesd_age40", z, m, x)]), ] %>%
  mutate(std_cesd_age40 = as.numeric(scale(cesd_age40)))

df <- nlsy

n  <- nrow(df)

#------------------------#
#     MODEL FORMULAE     #
#------------------------#

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

#------------------------#
#     MAIN ANALYSES      #
#------------------------#

estimands <- expand.grid(c(0, 1), c(0, 1)) %>%
  `colnames<-`(c("a1", "a2"))

S <- nrow(estimands)

#pre-compute design matrices used in both DML and MR
df_b  <- model.matrix(b_form, data = df)[, -1] %>% as_tibble()
df_g  <- model.matrix(g_form, data = df)[, -1] %>% as_tibble()
df_h  <- model.matrix(h_form, data = df)[, -1] %>% as_tibble()
df_q  <- model.matrix(q_form, data = df)[, -1] %>% as_tibble()
df_r  <- model.matrix(r_form, data = df)[, -1] %>% as_tibble()
df_u  <- model.matrix(u_form, data = df)[, -1] %>% as_tibble()

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

#define trimQ locally
trimQ <- function(x, low = 0.01, high = 0.99) {
  q <- quantile(x, c(low, high), na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

#====================================#
#      Parametric MR Estimation      #
#====================================#

#----------------------------------#
#       Parametric MR model        #
#----------------------------------#
b_glm <- lm (b_form, data = df)
g_glm <- glm(g_form, family = binomial("logit"), data = df)
h_glm <- glm(h_form, family = binomial("logit"), data = df)
q_glm <- glm(q_form, family = binomial("logit"), data = df)
r_glm <- glm(r_form, family = binomial("logit"), data = df)

#--------------------------------------------#
#    Parametric MR: EIF & point estimates    #
#--------------------------------------------#
df <- df %>% mutate(
  g0_fit = 1 - g_glm$fitted.values,
  g1_fit = g_glm$fitted.values,
  gA_fit = att22 * g1_fit + (1 - att22) * g0_fit,
  
  h0_fit = 1 - h_glm$fitted.values,
  h1_fit = h_glm$fitted.values,
  hA_fit = att22 * h1_fit + (1 - att22) * h0_fit,
  
  q01_fit = predict(q_glm, newdata = df_q0, type = "response"),
  q11_fit = predict(q_glm, newdata = df_q1, type = "response"),
  q00_fit = 1 - q01_fit, q10_fit = 1 - q11_fit,
  
  r01_fit = predict(r_glm, newdata = df_r0, type = "response"),
  r11_fit = predict(r_glm, newdata = df_r1, type = "response"),
  r00_fit = 1 - r01_fit, r10_fit = 1 - r11_fit,
  
  q0Z_fit = ifelse(ever_unemp_age3539 == 1, q01_fit, q00_fit),
  q1Z_fit = ifelse(ever_unemp_age3539 == 1, q11_fit, q10_fit),
  qAZ_fit = ifelse(att22 == 1, q1Z_fit, q0Z_fit),
  
  r0Z_fit = ifelse(ever_unemp_age3539 == 1, r01_fit, r00_fit),
  r1Z_fit = ifelse(ever_unemp_age3539 == 1, r11_fit, r10_fit),
  rAZ_fit = ifelse(att22 == 1, r1Z_fit, r0Z_fit),
  
  b0_fit  = predict(b_glm, newdata = df_b0, type = "response"),
  b1_fit  = predict(b_glm, newdata = df_b1, type = "response"),
  bA_fit  = ifelse(att22 == 1, b1_fit, b0_fit),
  
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
  dep_vmod_a2y = b10_fit * q10_fit + b11_fit * q11_fit
)

u_a1n_glm <- lm(as.formula(paste("dep_umod_a1n", " ~ ", paste(c(x, a, z), collapse = "+"))), data = df)
u_a1y_glm <- lm(as.formula(paste("dep_umod_a1y", " ~ ", paste(c(x, a, z), collapse = "+"))), data = df)
v_a2n_glm <- lm(as.formula(paste("dep_vmod_a2n", " ~ ", paste(c(x, a),    collapse = "+"))), data = df)
v_a2y_glm <- lm(as.formula(paste("dep_vmod_a2y", " ~ ", paste(c(x, a),    collapse = "+"))), data = df)

df <- df %>% mutate(
  
  u0_a1n_fit = predict(u_a1n_glm, newdata = df_u0),
  u1_a1n_fit = predict(u_a1n_glm, newdata = df_u1),
  
  u0_a1y_fit = predict(u_a1y_glm, newdata = df_u0),
  u1_a1y_fit = predict(u_a1y_glm, newdata = df_u1),
  
  v0_a2n_fit = predict(v_a2n_glm, newdata = df_q0),
  v1_a2n_fit = predict(v_a2n_glm, newdata = df_q1),
  
  v0_a2y_fit = predict(v_a2y_glm, newdata = df_q0),
  v1_a2y_fit = predict(v_a2y_glm, newdata = df_q1),
  
  u00_a1n_fit = predict(u_a1n_glm, newdata = df_u00),
  u01_a1n_fit = predict(u_a1n_glm, newdata = df_u01),
  u10_a1n_fit = predict(u_a1n_glm, newdata = df_u10),
  u11_a1n_fit = predict(u_a1n_glm, newdata = df_u11),
  
  u00_a1y_fit = predict(u_a1y_glm, newdata = df_u00),
  u01_a1y_fit = predict(u_a1y_glm, newdata = df_u01),
  u10_a1y_fit = predict(u_a1y_glm, newdata = df_u10),
  u11_a1y_fit = predict(u_a1y_glm, newdata = df_u11),
  
  uqOz_nn_fit = u00_a1n_fit * q00_fit +  u01_a1n_fit * q01_fit,
  uqOz_ny_fit = u10_a1n_fit * q10_fit +  u11_a1n_fit * q11_fit,
  uqOz_yn_fit = u00_a1y_fit * q00_fit +  u01_a1y_fit * q01_fit,
  uqOz_yy_fit = u10_a1y_fit * q10_fit +  u11_a1y_fit * q11_fit,
  
  bqOz_a2n_fit = b00_fit * q00_fit +  b01_fit * q01_fit,
  bqOz_a2y_fit = b10_fit * q10_fit +  b11_fit * q11_fit
  
)

df_mr <- df

for (s in 1:S) {
  
  a1 <- estimands$a1[[s]]; a2 <- estimands$a2[[s]]
  
  df_mr <- df_mr %>% mutate(
    g_a2_fit = a2 * g1_fit + (1 - a2) * g0_fit,
    g_a1_fit = a1 * g1_fit + (1 - a1) * g0_fit,
    
    c_fit = a1 * a2 * c11_fit + 
      a1 * (1-a2) * c10_fit + 
      (1-a1) * a2 * c01_fit + 
      (1-a1) * (1-a2) * c00_fit,
    
    b_fit = a2 * b1_fit + (1-a2) * b0_fit,
    
    u_fit = a1 * a2 * u1_a1y_fit + 
      a1 * (1-a2) * u0_a1y_fit + 
      (1-a1) * a2 * u1_a1n_fit + 
      (1-a1) * (1-a2) * u0_a1n_fit,
    
    v_fit = a1 * a2 * v1_a2y_fit + 
      a1 * (1-a2) * v1_a2n_fit + 
      (1-a1) * a2 * v0_a2y_fit + 
      (1-a1) * (1-a2) * v0_a2n_fit,
    
    uqOz_fit = a1 * a2 * uqOz_yy_fit + 
      a1 * (1-a2) * uqOz_yn_fit + 
      (1-a1) * a2 * uqOz_ny_fit + 
      (1-a1) * (1-a2) * uqOz_nn_fit,
    
    bqOz_fit = a2 * bqOz_a2y_fit + (1-a2) * bqOz_a2n_fit,
    
    !!sym(paste0("w2_", a1, a2)) := as.double(att22==a2)  *  c_fit/g_a2_fit,
    !!sym(paste0("w1_", a1, a2)) := as.double(att22==a2) /g_a2_fit,
    !!sym(paste0("w0_", a1, a2)) := as.double(att22==a1) /g_a1_fit
    
  )
  
  df_mr[df_mr$att22==a2, paste0("w2_", a1, a2)] <- trimQ(df_mr[df_mr$att22==a2, paste0("w2_", a1, a2)])
  df_mr[df_mr$att22==a2, paste0("w1_", a1, a2)] <- trimQ(df_mr[df_mr$att22==a2, paste0("w1_", a1, a2)])
  df_mr[df_mr$att22==a1, paste0("w0_", a1, a2)] <- trimQ(df_mr[df_mr$att22==a1, paste0("w0_", a1, a2)])
  
  df_mr <- df_mr %>% mutate(
    
    !!sym(paste0("eif_", a1, a2)) := !!sym(paste0("w2_", a1, a2)) * (std_cesd_age40 - b_fit) +
      !!sym(paste0("w1_", a1, a2)) * (u_fit - uqOz_fit) +
      !!sym(paste0("w0_", a1, a2)) * (bqOz_fit - v_fit) +
      v_fit
  )
}

out_df_par_base <- df_mr %>%
  mutate(eif_type1_ate = eif_11 - eif_00,
         eif_type2_ate = eif_11 - eif_00,
         eif_type1_iie = eif_10 - eif_00,
         eif_type1_ide = eif_11 - eif_10,
         eif_type2_iie = eif_11 - eif_01,
         eif_type2_ide = eif_01 - eif_00) %>%
  
  summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                            se  = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator","type","estimand","measure")) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  filter(type == "type1")

#-----------------------------#
#   Non-parametric bootstrap  #
#-----------------------------#
B <- 2000

boots <- matrix(NA_real_, nrow = B, ncol = 3)

for (b in 1:B) {
  
  if (b %% 100 == 0) cat(" bootstrap sample ", b, "\n")
  dfi <- nlsy %>% sample_frac(replace = TRUE)
  
  #design matrices for the different model
  dfi_b  <- model.matrix(b_form, data = dfi)[, -1] %>% as_tibble()
  dfi_g  <- model.matrix(g_form, data = dfi)[, -1] %>% as_tibble()
  dfi_h  <- model.matrix(h_form, data = dfi)[, -1] %>% as_tibble()
  dfi_q  <- model.matrix(q_form, data = dfi)[, -1] %>% as_tibble()
  dfi_r  <- model.matrix(r_form, data = dfi)[, -1] %>% as_tibble()
  dfi_u  <- model.matrix(u_form, data = dfi)[, -1] %>% as_tibble()
  
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
  
  #fit GLMs in the bootstrap draw
  b_glm <- lm (b_form, data = dfi)
  g_glm <- glm(g_form, family = binomial("logit"), data = dfi)
  h_glm <- glm(h_form, family = binomial("logit"), data = dfi)
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
    q00_fit = 1 - q01_fit, q10_fit = 1 - q11_fit,
    
    r01_fit = predict(r_glm, newdata = dfi_r0, type = "response"),
    r11_fit = predict(r_glm, newdata = dfi_r1, type = "response"),
    r00_fit = 1 - r01_fit, r10_fit = 1 - r11_fit,
    
    q0Z_fit = ifelse(ever_unemp_age3539 == 1, q01_fit, q00_fit),
    q1Z_fit = ifelse(ever_unemp_age3539 == 1, q11_fit, q10_fit),
    qAZ_fit = ifelse(att22 == 1, q1Z_fit, q0Z_fit),
    
    r0Z_fit = ifelse(ever_unemp_age3539 == 1, r01_fit, r00_fit),
    r1Z_fit = ifelse(ever_unemp_age3539 == 1, r11_fit, r10_fit),
    rAZ_fit = ifelse(att22 == 1, r1Z_fit, r0Z_fit),
    
    b0_fit  = predict(b_glm, newdata = dfi_b0, type = "response"),
    b1_fit  = predict(b_glm, newdata = dfi_b1, type = "response"),
    bA_fit  = ifelse(att22 == 1, b1_fit, b0_fit),
    
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
    dep_vmod_a2y = b10_fit * q10_fit + b11_fit * q11_fit
  )
  
  u_a1n_glm <- lm(as.formula(paste("dep_umod_a1n", " ~ ", paste(c(x, a, z), collapse = "+"))), data = dfi)
  u_a1y_glm <- lm(as.formula(paste("dep_umod_a1y", " ~ ", paste(c(x, a, z), collapse = "+"))), data = dfi)
  v_a2n_glm <- lm(as.formula(paste("dep_vmod_a2n", " ~ ", paste(c(x, a),    collapse = "+"))), data = dfi)
  v_a2y_glm <- lm(as.formula(paste("dep_vmod_a2y", " ~ ", paste(c(x, a),    collapse = "+"))), data = dfi)
  
  dfi <- dfi %>% mutate(
    u0_a1n_fit = predict(u_a1n_glm, newdata = dfi_u0),
    u1_a1n_fit = predict(u_a1n_glm, newdata = dfi_u1),
    
    u0_a1y_fit = predict(u_a1y_glm, newdata = dfi_u0),
    u1_a1y_fit = predict(u_a1y_glm, newdata = dfi_u1),
    
    v0_a2n_fit = predict(v_a2n_glm, newdata = dfi_q0),
    v1_a2n_fit = predict(v_a2n_glm, newdata = dfi_q1),
    
    v0_a2y_fit = predict(v_a2y_glm, newdata = dfi_q0),
    v1_a2y_fit = predict(v_a2y_glm, newdata = dfi_q1),
    
    u00_a1n_fit = predict(u_a1n_glm, newdata = dfi_u00),
    u01_a1n_fit = predict(u_a1n_glm, newdata = dfi_u01),
    
    u10_a1n_fit = predict(u_a1n_glm, newdata = dfi_u10),
    u11_a1n_fit = predict(u_a1n_glm, newdata = dfi_u11),
    
    u00_a1y_fit = predict(u_a1y_glm, newdata = dfi_u00),
    u01_a1y_fit = predict(u_a1y_glm, newdata = dfi_u01),
    
    u10_a1y_fit = predict(u_a1y_glm, newdata = dfi_u10),
    u11_a1y_fit = predict(u_a1y_glm, newdata = dfi_u11),
    
    uqOz_nn_fit = u00_a1n_fit * q00_fit + u01_a1n_fit * q01_fit,
    uqOz_ny_fit = u10_a1n_fit * q10_fit + u11_a1n_fit * q11_fit,
    
    uqOz_yn_fit = u00_a1y_fit * q00_fit + u01_a1y_fit * q01_fit,
    uqOz_yy_fit = u10_a1y_fit * q10_fit + u11_a1y_fit * q11_fit,
    
    bqOz_a2n_fit = b00_fit * q00_fit + b01_fit * q01_fit,
    bqOz_a2y_fit = b10_fit * q10_fit + b11_fit * q11_fit
  )
  
  for (s in 1:S) {
    a1 <- estimands$a1[[s]]; a2 <- estimands$a2[[s]]
    dfi <- dfi %>% mutate(
      g_a2_fit = a2 * g1_fit + (1 - a2) * g0_fit,
      g_a1_fit = a1 * g1_fit + (1 - a1) * g0_fit,
      
      c_fit = a1 * a2 * c11_fit + 
        a1 * (1-a2) * c10_fit + 
        (1-a1) * a2 * c01_fit + 
        (1-a1) * (1-a2) * c00_fit,
      
      b_fit = a2 * b1_fit + (1-a2) * b0_fit,
      
      u_fit = a1 * a2 * u1_a1y_fit + 
        a1 * (1-a2) * u0_a1y_fit + 
        (1-a1) * a2 * u1_a1n_fit + 
        (1-a1) * (1-a2) * u0_a1n_fit,
      
      v_fit = a1 * a2 * v1_a2y_fit + 
        a1 * (1-a2) * v1_a2n_fit + 
        (1-a1) * a2 * v0_a2y_fit + 
        (1-a1) * (1-a2) * v0_a2n_fit,
      
      uqOz_fit = a1 * a2 * uqOz_yy_fit + 
        a1 * (1-a2) * uqOz_yn_fit + 
        (1-a1) * a2 * uqOz_ny_fit + 
        (1-a1) * (1-a2) * uqOz_nn_fit,
      
      bqOz_fit = a2 * bqOz_a2y_fit 
      + (1-a2) * bqOz_a2n_fit,
      
      !!sym(paste0("w2_", a1, a2)) := as.double(att22==a2)  *  c_fit/g_a2_fit,
      !!sym(paste0("w1_", a1, a2)) := as.double(att22==a2) /g_a2_fit,
      !!sym(paste0("w0_", a1, a2)) := as.double(att22==a1) /g_a1_fit
    )
    
    dfi[dfi$att22==a2, paste0("w2_", a1, a2)] <- trimQ(dfi[dfi$att22==a2, paste0("w2_", a1, a2)])
    dfi[dfi$att22==a2, paste0("w1_", a1, a2)] <- trimQ(dfi[dfi$att22==a2, paste0("w1_", a1, a2)])
    dfi[dfi$att22==a1, paste0("w0_", a1, a2)] <- trimQ(dfi[dfi$att22==a1, paste0("w0_", a1, a2)])
    
    dfi <- dfi %>% mutate(
      
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
    
    summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x))) %>%
    pivot_longer(everything()) %>%
    separate(name, into = c("estimator","type","estimand","measure")) %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    filter(type == "type1")
  
  boots[b, ] <- out_dfi$est
}

out_df_par <- out_df_par_base %>%
  mutate(se_naive = se,
         se = apply(boots, 2, sd),
         lower = apply(boots, 2, quantile, 0.025),
         upper = apply(boots, 2, quantile, 0.975))

label_map <- c(ate = "OE(1,0)", ide = "IDE(1,0)", iie = "IIE(1,0)")

table6_3par <- out_df_par %>%
  mutate(Estimand = label_map[as.character(estimand)],
         lower = est - 1.96 * se, upper = est + 1.96 * se) %>%
  mutate(across(c(est, se, lower, upper), ~ round(.x, 3))) %>%
  mutate(intv = paste0("(", lower, ", ", upper, ")"),
         out  = paste(est, intv)) %>%
  dplyr::select(estimator, type, estimand, est, se, se_naive, out, Estimand)

#==============================#
#        DML Estimation        #
#==============================#

#----------------------------#
#  DML cross-fitting setup   #
#----------------------------#
K <- 5
cf_fold <- caret::createFolds(df$std_cesd_age40, K)
main_list <- vector(mode = "list", K)

for (k in 1:K) {
  
  cat(" cross-fitting fold ", k, "\n")
  
  # auxiliary and main data
  aux  <- df[-cf_fold[[k]], ]
  main <- df[ cf_fold[[k]], ]
  
  # design matrices by role
  aux_b <- df_b[-cf_fold[[k]], ]; main_b <- df_b[cf_fold[[k]], ]
  aux_g <- df_g[-cf_fold[[k]], ]; main_g <- df_g[cf_fold[[k]], ]
  aux_h <- df_h[-cf_fold[[k]], ]; main_h <- df_h[cf_fold[[k]], ]
  aux_q <- df_q[-cf_fold[[k]], ]; main_q <- df_q[cf_fold[[k]], ]
  aux_r <- df_r[-cf_fold[[k]], ]; main_r <- df_r[cf_fold[[k]], ]
  aux_u <- df_u[-cf_fold[[k]], ]; main_u <- df_u[cf_fold[[k]], ]
  
  #-----------------------------#
  #        DML estimation       #
  #-----------------------------#
  
  #outcome model
  b_sl <- SuperLearner(
    Y = aux$std_cesd_age40, 
    X = aux_b, 
    family = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  #treatment models
  g_sl <- SuperLearner(
    Y = aux$att22, 
    X = aux_g, 
    newX = df_g, 
    family = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  h_sl <- SuperLearner(
    Y = aux$att22, 
    X = aux_h, 
    newX = df_h, 
    family = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  #z models
  q_sl <- SuperLearner(
    Y = aux$ever_unemp_age3539, 
    X = aux_q, 
    family = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  r_sl <- SuperLearner(
    Y = aux$ever_unemp_age3539, 
    X = aux_r, 
    family = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  #---------------------------------#
  #           DML plug-in           #
  #---------------------------------#
  df <- df %>% mutate(
    g0_fit = 1 - g_sl$SL.predict,
    g1_fit = g_sl$SL.predict,
    gA_fit = ifelse(att22 == 1, g1_fit, g0_fit),
    
    h0_fit = 1 - h_sl$SL.predict,
    h1_fit = h_sl$SL.predict,
    hA_fit = ifelse(att22 == 1, h1_fit, h0_fit),
    
    q01_fit = predict.SuperLearner(q_sl, newdata = df_q0)$pred,
    q11_fit = predict.SuperLearner(q_sl, newdata = df_q1)$pred,
    q00_fit = 1 - q01_fit, q10_fit = 1 - q11_fit,
    
    r01_fit = predict.SuperLearner(r_sl, newdata = df_r0)$pred,
    r11_fit = predict.SuperLearner(r_sl, newdata = df_r1)$pred,
    r00_fit = 1 - r01_fit, r10_fit = 1 - r11_fit,
    
    q0Z_fit = ifelse(ever_unemp_age3539 == 1, q01_fit, q00_fit),
    q1Z_fit = ifelse(ever_unemp_age3539 == 1, q11_fit, q10_fit),
    qAZ_fit = ifelse(att22 == 1, q1Z_fit, q0Z_fit),
    
    r0Z_fit = ifelse(ever_unemp_age3539 == 1, r01_fit, r00_fit),
    r1Z_fit = ifelse(ever_unemp_age3539 == 1, r11_fit, r10_fit),
    rAZ_fit = ifelse(att22 == 1, r1Z_fit, r0Z_fit),
    
    b0_fit  = as.numeric(predict.SuperLearner(b_sl, newdata = df_b0)$pred),
    b1_fit  = as.numeric(predict.SuperLearner(b_sl, newdata = df_b1)$pred),
    bA_fit  = ifelse(att22 == 1, b1_fit, b0_fit),
    
    b00_fit = as.numeric(predict.SuperLearner(b_sl, newdata = df_b00)$pred),
    b01_fit = as.numeric(predict.SuperLearner(b_sl, newdata = df_b01)$pred),
    b10_fit = as.numeric(predict.SuperLearner(b_sl, newdata = df_b10)$pred),
    b11_fit = as.numeric(predict.SuperLearner(b_sl, newdata = df_b11)$pred),
    
    c_a1n_fit = gA_fit/g0_fit * qAZ_fit/rAZ_fit * h0_fit/hA_fit,
    c_a1y_fit = gA_fit/g1_fit * qAZ_fit/rAZ_fit * h1_fit/hA_fit,
    
    c00_fit = q0Z_fit/r0Z_fit,
    c01_fit = g1_fit/g0_fit * q1Z_fit/r1Z_fit * h0_fit/h1_fit,
    c10_fit = g0_fit/g1_fit * q0Z_fit/r0Z_fit * h1_fit/h0_fit,
    c11_fit = q1Z_fit/r1Z_fit,
    
    dep_umod_a1n = bA_fit * c_a1n_fit,
    dep_umod_a1y = bA_fit * c_a1y_fit,
    dep_vmod_a2n = b00_fit * q00_fit + b01_fit * q01_fit,
    dep_vmod_a2y = b10_fit * q10_fit + b11_fit * q11_fit
  )
  
  # u- and v-models for DML targets
  u_a1n_sl <- SuperLearner(
    Y = df$dep_umod_a1n[-cf_fold[[k]]], 
    X = aux_u, 
    family = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  u_a1y_sl <- SuperLearner(
    Y = df$dep_umod_a1y[-cf_fold[[k]]], 
    X = aux_u, 
    family = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  v_a2n_sl <- SuperLearner(
    Y = df$dep_vmod_a2n[-cf_fold[[k]]], 
    X = aux_q, 
    family = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  v_a2y_sl <- SuperLearner(
    Y = df$dep_vmod_a2y[-cf_fold[[k]]], 
    X = aux_q, 
    family = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df <- df %>% mutate(
    u0_a1n_fit = predict.SuperLearner(u_a1n_sl, newdata = df_u0)$pred,
    u1_a1n_fit = predict.SuperLearner(u_a1n_sl, newdata = df_u1)$pred,
    
    u0_a1y_fit = predict.SuperLearner(u_a1y_sl, newdata = df_u0)$pred,
    u1_a1y_fit = predict.SuperLearner(u_a1y_sl, newdata = df_u1)$pred,
    
    v0_a2n_fit = predict.SuperLearner(v_a2n_sl, newdata = df_q0)$pred,
    v1_a2n_fit = predict.SuperLearner(v_a2n_sl, newdata = df_q1)$pred,
    
    v0_a2y_fit = predict.SuperLearner(v_a2y_sl, newdata = df_q0)$pred,
    v1_a2y_fit = predict.SuperLearner(v_a2y_sl, newdata = df_q1)$pred,
    
    u00_a1n_fit = predict.SuperLearner(u_a1n_sl, newdata = df_u00)$pred,
    u01_a1n_fit = predict.SuperLearner(u_a1n_sl, newdata = df_u01)$pred,
    u10_a1n_fit = predict.SuperLearner(u_a1n_sl, newdata = df_u10)$pred,
    u11_a1n_fit = predict.SuperLearner(u_a1n_sl, newdata = df_u11)$pred,
    
    u00_a1y_fit = predict.SuperLearner(u_a1y_sl, newdata = df_u00)$pred,
    u01_a1y_fit = predict.SuperLearner(u_a1y_sl, newdata = df_u01)$pred,
    u10_a1y_fit = predict.SuperLearner(u_a1y_sl, newdata = df_u10)$pred,
    u11_a1y_fit = predict.SuperLearner(u_a1y_sl, newdata = df_u11)$pred,
    
    uqOz_nn_fit = u00_a1n_fit * q00_fit + u01_a1n_fit * q01_fit,
    uqOz_ny_fit = u10_a1n_fit * q10_fit + u11_a1n_fit * q11_fit,
    uqOz_yn_fit = u00_a1y_fit * q00_fit + u01_a1y_fit * q01_fit,
    uqOz_yy_fit = u10_a1y_fit * q10_fit + u11_a1y_fit * q11_fit,
    
    bqOz_a2n_fit = b00_fit * q00_fit + b01_fit * q01_fit,
    bqOz_a2y_fit = b10_fit * q10_fit + b11_fit * q11_fit
  )
  
  main_list[[k]] <- df[cf_fold[[k]], ]
}

main_df <- purrr::reduce(main_list, bind_rows)
for (s in 1:S) {
  a1 <- estimands$a1[[s]]; a2 <- estimands$a2[[s]]
  main_df <- main_df %>% mutate(
    g_a2_fit = a2 * g1_fit + (1 - a2) * g0_fit,
    g_a1_fit = a1 * g1_fit + (1 - a1) * g0_fit,
    
    c_fit = a1 * a2 * c11_fit +
      a1 * (1-a2) * c10_fit + 
      (1-a1) * a2 * c01_fit + 
      (1-a1) * (1-a2) * c00_fit,
    
    b_fit = a2 * b1_fit + (1-a2) * b0_fit,
    
    u_fit = a1 * a2 * u1_a1y_fit +
      a1 * (1-a2) * u0_a1y_fit + 
      (1-a1) * a2 * u1_a1n_fit + 
      (1-a1) * (1-a2) * u0_a1n_fit,
    
    v_fit = a1 * a2 * v1_a2y_fit + 
      a1 * (1-a2) * v1_a2n_fit + 
      (1-a1) * a2 * v0_a2y_fit + 
      (1-a1) * (1-a2) * v0_a2n_fit,
    
    uqOz_fit = a1 * a2 * uqOz_yy_fit + 
      a1 * (1-a2) * uqOz_yn_fit + 
      (1-a1) * a2 * uqOz_ny_fit + 
      (1-a1) * (1-a2) * uqOz_nn_fit,
    
    bqOz_fit = a2 * bqOz_a2y_fit + (1-a2) * bqOz_a2n_fit,
    
    !!sym(paste0("w2_", a1, a2)) := as.double(att22==a2) * c_fit/g_a2_fit,
    !!sym(paste0("w1_", a1, a2)) := as.double(att22==a2) /g_a2_fit,
    !!sym(paste0("w0_", a1, a2)) := as.double(att22==a1) /g_a1_fit
  )
  
  main_df[main_df$att22==a2, paste0("w2_", a1, a2)] <- trimQ(main_df[main_df$att22==a2, paste0("w2_", a1, a2)])
  main_df[main_df$att22==a2, paste0("w1_", a1, a2)] <- trimQ(main_df[main_df$att22==a2, paste0("w1_", a1, a2)])
  main_df[main_df$att22==a1, paste0("w0_", a1, a2)] <- trimQ(main_df[main_df$att22==a1, paste0("w0_", a1, a2)])
  
  main_df <- main_df %>% mutate(
    
    !!sym(paste0("eif_", a1, a2)) := !!sym(paste0("w2_", a1, a2)) * (std_cesd_age40 - b_fit) +
      !!sym(paste0("w1_", a1, a2)) * (u_fit - uqOz_fit) +
      !!sym(paste0("w0_", a1, a2)) * (bqOz_fit - v_fit) +
      v_fit
  )
}

out_df_dml <- main_df %>%
  mutate(eif_type1_ate = eif_11 - eif_00,
         eif_type2_ate = eif_11 - eif_00,
         eif_type1_iie = eif_10 - eif_00,
         eif_type1_ide = eif_11 - eif_10,
         eif_type2_iie = eif_11 - eif_01,
         eif_type2_ide = eif_01 - eif_00) %>%
  
  summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                            se  = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator","type","estimand","measure")) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  filter(type == "type1")

#table for dml
table6_3np <- out_df_dml %>%
  mutate(Estimand = label_map[as.character(estimand)],
         lower = est - 1.96 * se, upper = est + 1.96 * se) %>%
  mutate(across(c(est, se, lower, upper), ~ round(.x, 3))) %>%
  mutate(intv = paste0("(", lower, ", ", upper, ")"),
         out  = paste(est, intv)) %>%
  dplyr::select(estimator, type, estimand, est, se, out, Estimand)

#-------------------------------#
#        COMBINE RESULTS        #
#-------------------------------#

par_tbl <- table6_3par %>%
  dplyr::filter(type == "type1") %>%
  dplyr::transmute(Estimand, `Parametric MR` = out)

dml_tbl <- table6_3np %>%
  dplyr::filter(type == "type1") %>%
  dplyr::transmute(Estimand, DML = out)

final_rst <- tibble::tibble(Estimand = c("OE(1,0)", "IDE(1,0)", "IIE(1,0)")) %>%
  dplyr::left_join(par_tbl, by = "Estimand") %>%
  dplyr::left_join(dml_tbl, by = "Estimand")

# Open log
sink(log_path, split = TRUE)

width_curr <- getOption("width"); options(width = 300)
print(final_rst, right = FALSE, row.names = FALSE)
options(width = width_curr)

# Close log
sink()
