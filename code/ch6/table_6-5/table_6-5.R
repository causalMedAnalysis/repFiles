#----------Preliminaries----------#
rm(list = ls())
chapter <- "ch6"
title   <- "table_6-5"

# Specify the root directory:
dir_root <- "C:/Users/Geoffrey Wodtke/Dropbox/D/projects/causal_mediation_text"

# Define subdirectories for logs and figures:
dir_log  <- paste0(dir_root, "/code/", chapter, "/_LOGS")
log_path <- paste0(dir_log,  "/", title, "_log.txt")

# Ensure required folders exist
create_dir_if_missing <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

create_dir_if_missing(dir_root)
create_dir_if_missing(dir_log)

#-------------------------------------------------------------------------------
#  Causal Mediation Analysis Replication Files
#
#  GitHub Repo: https://github.com/causalMedAnalysis/repFiles/tree/main
#
#  Script:      .../code/ch6/table_6-5.R
#
#  Inputs:      https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/Tatar/tatar.rds
#
#  Outputs:     .../code/ch6/_LOGS/table_6-5_log.txt
#
#  Description: Replicates Chapter 6, Table 6.5: Multiply Robust Estimates for the Total
#               and Path-Specific Effects of Ancestor Victimization on Regime Support as
#               Mediated by Political Identities
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
  "devtools"
)

install_and_load <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      install.packages(p, dependencies = TRUE)
    }
    library(p, character.only = TRUE)
  }
}
install_and_load(packages)

install_github("causalMedAnalysis/cmedR")

library(cmedR)

#----------------------#
#    SPECIFICATIONS    #
#----------------------#

# exposure
a <- "violence"

# outcome
y <- "annex"

# mediator1
m1 <- c("trust_g1", "victim_g1", "fear_g1")
# mediator2
m2 <- c("trust_g2", "victim_g2", "fear_g2")
# mediator3
m3 <- c("trust_g3", "victim_g3", "fear_g3")
m <- list(m1, m2, m3)

# baseline confounders
x <- c("kulak", 
       "prosoviet_pre", 
       "religiosity_pre", 
       "land_pre",
       "orchard_pre", 
       "animals_pre", 
       "carriage_pre", 
       "otherprop_pre")

#set seed
set.seed(02138)

#-----------------------------#
#        PREPARE DATA         #
#-----------------------------#

tatar<-
  as.data.frame(
    readRDS(
      url("https://raw.githubusercontent.com/causalMedAnalysis/repFiles/refs/heads/main/data/Tatar/tatar.rds")
    )
  )

df <- tatar
n <- nrow(df)

#-------------------------------------------------#
#   FORMULAS FOR TREATMENT AND OUTCOME MODELS     #
#-------------------------------------------------#

a0_form <- as.formula(paste(a, " ~ ", paste(x, collapse= "+")))
a1_form <- as.formula(paste(a, " ~ ", paste(c(x, m1), collapse= "+")))
a2_form <- as.formula(paste(a, " ~ ", paste(c(x, m1, m2), collapse= "+")))
a3_form <- as.formula(paste(a, " ~ ", paste(c(x, m1, m2, m3), collapse= "+")))

y0_form <- as.formula(paste(y, " ~ ", paste(c(x, a), collapse= "+")))
y1_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1), collapse= "+")))
y2_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1, m2), collapse= "+")))
y3_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1, m2, m3), collapse= "+")))

#-----------------------------#
#       MAIN ANALYSES         #
#-----------------------------#

estimands <- expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1)) %>%
  `colnames<-`(c("a1", "a2", "a3", "a4"))

S <- nrow(estimands)

#pre-compute design matrices used in both parametric MR and DML
df_p0 <- model.matrix(a0_form, data = df)[, -1] %>% as_tibble()
df_p1 <- model.matrix(a1_form, data = df)[, -1] %>% as_tibble()
df_p2 <- model.matrix(a2_form, data = df)[, -1] %>% as_tibble()
df_p3 <- model.matrix(a3_form, data = df)[, -1] %>% as_tibble()

df_mu0 <- model.matrix(y0_form, data = df)[, -1] %>% as_tibble()
df_mu1 <- model.matrix(y1_form, data = df)[, -1] %>% as_tibble()
df_mu2 <- model.matrix(y2_form, data = df)[, -1] %>% as_tibble()
df_mu3 <- model.matrix(y3_form, data = df)[, -1] %>% as_tibble()

df_mu3n <- model.matrix(y3_form, data = mutate(df, violence = 0))[, -1] %>% as_tibble()
df_mu2n <- model.matrix(y2_form, data = mutate(df, violence = 0))[, -1] %>% as_tibble()
df_mu1n <- model.matrix(y1_form, data = mutate(df, violence = 0))[, -1] %>% as_tibble()
df_mu0n <- model.matrix(y0_form, data = mutate(df, violence = 0))[, -1] %>% as_tibble()

df_mu3y <- model.matrix(y3_form, data = mutate(df, violence = 1))[, -1] %>% as_tibble()
df_mu2y <- model.matrix(y2_form, data = mutate(df, violence = 1))[, -1] %>% as_tibble()
df_mu1y <- model.matrix(y1_form, data = mutate(df, violence = 1))[, -1] %>% as_tibble()
df_mu0y <- model.matrix(y0_form, data = mutate(df, violence = 1))[, -1] %>% as_tibble()

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

#--------------------------#
#      Treatment Model     #
#--------------------------#

p0_glm <- glm(a0_form, family = binomial("logit"), data = df)
p1_glm <- glm(a1_form, family = binomial("logit"), data = df)
p2_glm <- glm(a2_form, family = binomial("logit"), data = df)
p3_glm <- glm(a3_form, family = binomial("logit"), data = df)

df <- df %>% mutate(
  p0_fit = p0_glm$fitted.values,
  p1_fit = p1_glm$fitted.values,
  p2_fit = p2_glm$fitted.values,
  p3_fit = p3_glm$fitted.values,
  
  w0_n = I(violence == 0)/(1 - p0_fit),
  w0_y = I(violence == 1)/p0_fit,
  w1_nn = I(violence == 0)/(1 - p0_fit) * 1,
  w1_ny = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit,
  w1_yn = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit),
  w1_yy = I(violence == 1)/p0_fit * 1,
  w2_000 = I(violence == 0)/(1 - p0_fit) * 1 * 1,
  w2_001 = I(violence == 1)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit,
  w2_010 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit),
  w2_011 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1,
  w2_100 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * 1,
  w2_101 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit,
  w2_110 = I(violence == 0)/p0_fit * 1 * p2_fit/(1 - p2_fit),
  w2_111 = I(violence == 1)/p0_fit * 1 * 1,
  
  w3_0000 = I(violence == 0)/(1 - p0_fit) * 1 * 1 * 1,
  w3_0001 = I(violence == 1)/(1 - p0_fit) * 1 * 1 * (1 - p3_fit)/p3_fit,
  w3_0010 = I(violence == 0)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit * p3_fit/(1 - p3_fit),
  w3_0011 = I(violence == 1)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit * 1,
  
  w3_0100 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit) * 1,
  w3_0101 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit) * (1 - p3_fit)/p3_fit,
  w3_0110 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1 * p3_fit/(1 - p3_fit),
  w3_0111 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1 * 1,
  
  w3_1000 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * 1 * 1,
  w3_1001 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * 1 * (1 - p3_fit)/p3_fit,
  w3_1010 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit * p3_fit/(1 - p3_fit),
  w3_1011 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit * 1,
  
  w3_1100 = I(violence == 0)/p0_fit * 1 * p2_fit/(1 - p2_fit) * 1,
  w3_1101 = I(violence == 1)/p0_fit * 1 * p2_fit/(1 - p2_fit) * (1 - p3_fit)/p3_fit,
  w3_1110 = I(violence == 0)/p0_fit * 1 * 1 * p3_fit/(1 - p3_fit),
  w3_1111 = I(violence == 1)/p0_fit * 1 * 1 * 1,
)

#-----------------------#
#     Outcome Model     #
#-----------------------#

mu3_glm <- glm(y3_form, family = binomial("logit"), data = df)

df$mu3_fit_a4n <- predict(mu3_glm, type = "response", newdata = df_mu3n)
df$mu3_fit_a4y <- predict(mu3_glm, type = "response", newdata = df_mu3y)

mu3_fit_a4y_form <- mu3_fit_a4n_form <- y2_form
mu3_fit_a4n_form[[2]] <- expr(mu3_fit_a4n)
mu3_fit_a4y_form[[2]] <- expr(mu3_fit_a4y)

mu2_glm_a4n <- glm(mu3_fit_a4n_form, family = quasibinomial("logit"), data = df)
mu2_glm_a4y <- glm(mu3_fit_a4y_form, family = quasibinomial("logit"), data = df)

df$mu2_fit_a4n_a3n <- predict(mu2_glm_a4n, type = "response", newdata = df_mu2n)
df$mu2_fit_a4n_a3y <- predict(mu2_glm_a4n, type = "response", newdata = df_mu2y)
df$mu2_fit_a4y_a3n <- predict(mu2_glm_a4y, type = "response", newdata = df_mu2n)
df$mu2_fit_a4y_a3y <- predict(mu2_glm_a4y, type = "response", newdata = df_mu2y)

mu2_fit_a4y_a3y_form <- mu2_fit_a4y_a3n_form <- y1_form
mu2_fit_a4n_a3y_form <- mu2_fit_a4n_a3n_form <- y1_form

mu2_fit_a4y_a3y_form[[2]] <- expr(mu2_fit_a4y_a3y)
mu2_fit_a4y_a3n_form[[2]] <- expr(mu2_fit_a4y_a3n)
mu2_fit_a4n_a3y_form[[2]] <- expr(mu2_fit_a4n_a3y)
mu2_fit_a4n_a3n_form[[2]] <- expr(mu2_fit_a4n_a3n)

mu1_glm_a4n_a3n <- glm(mu2_fit_a4n_a3n_form, family = quasibinomial("logit"), data = df)
mu1_glm_a4n_a3y <- glm(mu2_fit_a4n_a3y_form, family = quasibinomial("logit"), data = df)
mu1_glm_a4y_a3n <- glm(mu2_fit_a4y_a3n_form, family = quasibinomial("logit"), data = df)
mu1_glm_a4y_a3y <- glm(mu2_fit_a4y_a3y_form, family = quasibinomial("logit"), data = df)

df$mu1_fit_a4n_a3n_a2n <- predict(mu1_glm_a4n_a3n, type = "response", newdata = df_mu1n)
df$mu1_fit_a4n_a3n_a2y <- predict(mu1_glm_a4n_a3n, type = "response", newdata = df_mu1y)

df$mu1_fit_a4n_a3y_a2n <- predict(mu1_glm_a4n_a3y, type = "response", newdata = df_mu1n)
df$mu1_fit_a4n_a3y_a2y <- predict(mu1_glm_a4n_a3y, type = "response", newdata = df_mu1y)

df$mu1_fit_a4y_a3n_a2n <- predict(mu1_glm_a4y_a3n, type = "response", newdata = df_mu1n)
df$mu1_fit_a4y_a3n_a2y <- predict(mu1_glm_a4y_a3n, type = "response", newdata = df_mu1y)

df$mu1_fit_a4y_a3y_a2n <- predict(mu1_glm_a4y_a3y, type = "response", newdata = df_mu1n)
df$mu1_fit_a4y_a3y_a2y <- predict(mu1_glm_a4y_a3y, type = "response", newdata = df_mu1y)

mu1_fit_a4y_a3y_a2y_form <- mu1_fit_a4y_a3n_a2y_form <- y0_form
mu1_fit_a4n_a3y_a2y_form <- mu1_fit_a4n_a3n_a2y_form <- y0_form
mu1_fit_a4y_a3y_a2n_form <- mu1_fit_a4y_a3n_a2n_form <- y0_form
mu1_fit_a4n_a3y_a2n_form <- mu1_fit_a4n_a3n_a2n_form <- y0_form

mu1_fit_a4y_a3y_a2y_form[[2]] <- expr(mu1_fit_a4y_a3y_a2y)
mu1_fit_a4y_a3y_a2n_form[[2]] <- expr(mu1_fit_a4y_a3y_a2n)
mu1_fit_a4y_a3n_a2y_form[[2]] <- expr(mu1_fit_a4y_a3n_a2y)
mu1_fit_a4y_a3n_a2n_form[[2]] <- expr(mu1_fit_a4y_a3n_a2n)

mu1_fit_a4n_a3y_a2y_form[[2]] <- expr(mu1_fit_a4n_a3y_a2y)
mu1_fit_a4n_a3y_a2n_form[[2]] <- expr(mu1_fit_a4n_a3y_a2n)
mu1_fit_a4n_a3n_a2y_form[[2]] <- expr(mu1_fit_a4n_a3n_a2y)
mu1_fit_a4n_a3n_a2n_form[[2]] <- expr(mu1_fit_a4n_a3n_a2n)

mu0_glm_a4n_a3n_a2n <- glm(mu1_fit_a4n_a3n_a2n_form, family = quasibinomial("logit"), data = df)
mu0_glm_a4n_a3n_a2y <- glm(mu1_fit_a4n_a3n_a2y_form, family = quasibinomial("logit"), data = df)
mu0_glm_a4n_a3y_a2n <- glm(mu1_fit_a4n_a3y_a2n_form, family = quasibinomial("logit"), data = df)
mu0_glm_a4n_a3y_a2y <- glm(mu1_fit_a4n_a3y_a2y_form, family = quasibinomial("logit"), data = df)

mu0_glm_a4y_a3n_a2n <- glm(mu1_fit_a4y_a3n_a2n_form, family = quasibinomial("logit"), data = df)
mu0_glm_a4y_a3n_a2y <- glm(mu1_fit_a4y_a3n_a2y_form, family = quasibinomial("logit"), data = df)
mu0_glm_a4y_a3y_a2n <- glm(mu1_fit_a4y_a3y_a2n_form, family = quasibinomial("logit"), data = df)
mu0_glm_a4y_a3y_a2y <- glm(mu1_fit_a4y_a3y_a2y_form, family = quasibinomial("logit"), data = df)

df$mu0_fit_a4n_a3n_a2n_a1n <- predict(mu0_glm_a4n_a3n_a2n, type = "response", newdata = df_mu0n)
df$mu0_fit_a4n_a3n_a2n_a1y <- predict(mu0_glm_a4n_a3n_a2n, type = "response", newdata = df_mu0y)

df$mu0_fit_a4n_a3n_a2y_a1n <- predict(mu0_glm_a4n_a3n_a2y, type = "response", newdata = df_mu0n)
df$mu0_fit_a4n_a3n_a2y_a1y <- predict(mu0_glm_a4n_a3n_a2y, type = "response", newdata = df_mu0y)

df$mu0_fit_a4n_a3y_a2n_a1n <- predict(mu0_glm_a4n_a3y_a2n, type = "response", newdata = df_mu0n)
df$mu0_fit_a4n_a3y_a2n_a1y <- predict(mu0_glm_a4n_a3y_a2n, type = "response", newdata = df_mu0y)

df$mu0_fit_a4n_a3y_a2y_a1n <- predict(mu0_glm_a4n_a3y_a2y, type = "response", newdata = df_mu0n)
df$mu0_fit_a4n_a3y_a2y_a1y <- predict(mu0_glm_a4n_a3y_a2y, type = "response", newdata = df_mu0y)

df$mu0_fit_a4y_a3n_a2n_a1n <- predict(mu0_glm_a4y_a3n_a2n, type = "response", newdata = df_mu0n)
df$mu0_fit_a4y_a3n_a2n_a1y <- predict(mu0_glm_a4y_a3n_a2n, type = "response", newdata = df_mu0y)

df$mu0_fit_a4y_a3n_a2y_a1n <- predict(mu0_glm_a4y_a3n_a2y, type = "response", newdata = df_mu0n)
df$mu0_fit_a4y_a3n_a2y_a1y <- predict(mu0_glm_a4y_a3n_a2y, type = "response", newdata = df_mu0y)

df$mu0_fit_a4y_a3y_a2n_a1n <- predict(mu0_glm_a4y_a3y_a2n, type = "response", newdata = df_mu0n)
df$mu0_fit_a4y_a3y_a2n_a1y <- predict(mu0_glm_a4y_a3y_a2n, type = "response", newdata = df_mu0y)

df$mu0_fit_a4y_a3y_a2y_a1n <- predict(mu0_glm_a4y_a3y_a2y, type = "response", newdata = df_mu0n)
df$mu0_fit_a4y_a3y_a2y_a1y <- predict(mu0_glm_a4y_a3y_a2y, type = "response", newdata = df_mu0y)

for (s in 1:S){
  
  a1 <- estimands$a1[[s]]
  a2 <- estimands$a2[[s]]
  a3 <- estimands$a3[[s]]
  a4 <- estimands$a4[[s]]
  
  df <- df %>%
    mutate(
      
      wt0_deno = a1 * p0_fit + (1 - a1) * (1 - p0_fit),
      wt1_nume = a1 * p1_fit + (1 - a1) * (1 - p1_fit),
      wt1_deno = a2 * p1_fit + (1 - a2) * (1 - p1_fit),
      wt2_nume = a2 * p2_fit + (1 - a2) * (1 - p2_fit),
      wt2_deno = a3 * p2_fit + (1 - a3) * (1 - p2_fit),
      wt3_nume = a3 * p3_fit + (1 - a3) * (1 - p3_fit),
      wt3_deno = a4 * p3_fit + (1 - a4) * (1 - p3_fit),
      
      !!sym(paste0("w0_", a1, a2, a3, a4)) := as.double(violence==a1) / wt0_deno,
      !!sym(paste0("w1_", a1, a2, a3, a4)) := as.double(violence==a2) * wt1_nume/wt1_deno/wt0_deno,
      !!sym(paste0("w2_", a1, a2, a3, a4)) := as.double(violence==a3) * wt2_nume/wt2_deno * wt1_nume/wt1_deno/wt0_deno,
      !!sym(paste0("w3_", a1, a2, a3, a4)) := as.double(violence==a4) * wt3_nume/wt3_deno * wt2_nume/wt2_deno * wt1_nume/wt1_deno/wt0_deno
    )
  
  df[df$violence == a1, paste0("w0_", a1, a2, a3, a4)] <- trimQ(df[df$violence == a1, paste0("w0_", a1, a2, a3, a4)])
  df[df$violence == a2, paste0("w1_", a1, a2, a3, a4)] <- trimQ(df[df$violence == a2, paste0("w1_", a1, a2, a3, a4)])
  df[df$violence == a3, paste0("w2_", a1, a2, a3, a4)] <- trimQ(df[df$violence == a3, paste0("w2_", a1, a2, a3, a4)])
  df[df$violence == a4, paste0("w3_", a1, a2, a3, a4)] <- trimQ(df[df$violence == a4, paste0("w3_", a1, a2, a3, a4)])
  
  df <- df %>%
    mutate(
      
      !!sym(paste0("mu3fit_", a1, a2, a3, a4)) := a4 * mu3_fit_a4y + (1 - a4) * mu3_fit_a4n,
      
      !!sym(paste0("mu2fit_", a1, a2, a3, a4)) := a4 * a3 * mu2_fit_a4y_a3y + a4 * (1 - a3) * mu2_fit_a4y_a3n +
        (1 - a4) * a3 * mu2_fit_a4n_a3y + (1 - a4) * (1 - a3) * mu2_fit_a4n_a3n,
      
      !!sym(paste0("mu1fit_", a1, a2, a3, a4)) := a4 * a3 * a2 * mu1_fit_a4y_a3y_a2y + a4 * a3 * (1 - a2) * mu1_fit_a4y_a3y_a2n +
        a4 * (1 - a3) * a2 * mu1_fit_a4y_a3n_a2y + a4 * (1 - a3) * (1 - a2) * mu1_fit_a4y_a3n_a2n +
        (1 - a4) * a3 * a2 * mu1_fit_a4n_a3y_a2y + (1 - a4) * a3 * (1 - a2) * mu1_fit_a4n_a3y_a2n +
        (1 - a4) * (1 - a3) * a2 * mu1_fit_a4n_a3n_a2y + (1 - a4) * (1 - a3) * (1 - a2) * mu1_fit_a4n_a3n_a2n,
      
      !!sym(paste0("mu0fit_", a1, a2, a3, a4)) := a4 * a3 * a2 * a1 * mu0_fit_a4y_a3y_a2y_a1y + a4 * a3 * (1 - a2) * a1 * mu0_fit_a4y_a3y_a2n_a1y +
        a4 * (1 - a3) * a2 * a1 * mu0_fit_a4y_a3n_a2y_a1y + a4 * (1 - a3) * (1 - a2) * a1 * mu0_fit_a4y_a3n_a2n_a1y +
        (1 - a4) * a3 * a2 * a1 * mu0_fit_a4n_a3y_a2y_a1y + (1 - a4) * a3 * (1 - a2) * a1 * mu0_fit_a4n_a3y_a2n_a1y +
        (1 - a4) * (1 - a3) * a2 * a1 * mu0_fit_a4n_a3n_a2y_a1y + (1 - a4) * (1 - a3) * (1 - a2) * a1 * mu0_fit_a4n_a3n_a2n_a1y +
        a4 * a3 * a2 * (1 - a1) * mu0_fit_a4y_a3y_a2y_a1n + a4 * a3 * (1 - a2) * (1 - a1) * mu0_fit_a4y_a3y_a2n_a1n +
        a4 * (1 - a3) * a2 * (1 - a1) * mu0_fit_a4y_a3n_a2y_a1n + a4 * (1 - a3) * (1 - a2) * (1 - a1) * mu0_fit_a4y_a3n_a2n_a1n +
        (1 - a4) * a3 * a2 * (1 - a1) * mu0_fit_a4n_a3y_a2y_a1n + (1 - a4) * a3 * (1 - a2) * (1 - a1) * mu0_fit_a4n_a3y_a2n_a1n +
        (1 - a4) * (1 - a3) * a2 * (1 - a1) * mu0_fit_a4n_a3n_a2y_a1n + (1 - a4) * (1 - a3) * (1 - a2) * (1 - a1) * mu0_fit_a4n_a3n_a2n_a1n,
      
      !!sym(paste0("www_", a1, a2, a3, a4)) := !!sym(paste0("w3_", a1, a2, a3, a4)) * annex,
      !!sym(paste0("iii_", a1, a2, a3, a4)) := !!sym(paste0("mu0fit_", a1, a2, a3, a4)),
      !!sym(paste0("eif_", a1, a2, a3, a4)) := !!sym(paste0("w3_", a1, a2, a3, a4)) *
        (annex - !!sym(paste0("mu3fit_", a1, a2, a3, a4))) +
        !!sym(paste0("w2_", a1, a2, a3, a4)) * (!!sym(paste0("mu3fit_", a1, a2, a3, a4)) -
                                                  !!sym(paste0("mu2fit_", a1, a2, a3, a4))) +
        !!sym(paste0("w1_", a1, a2, a3, a4)) * (!!sym(paste0("mu2fit_", a1, a2, a3, a4)) -
                                                  !!sym(paste0("mu1fit_", a1, a2, a3, a4))) +
        !!sym(paste0("w0_", a1, a2, a3, a4)) * (!!sym(paste0("mu1fit_", a1, a2, a3, a4)) -
                                                  !!sym(paste0("mu0fit_", a1, a2, a3, a4))) +
        !!sym(paste0("mu0fit_", a1, a2, a3, a4))
    )
}

out_df <- df %>%
  mutate(
    eif_type1_ate = eif_1111 - eif_0000,
    eif_type1_pse4 = eif_0001 - eif_0000,
    eif_type1_pse3 = eif_0011 - eif_0001,
    eif_type1_pse2 = eif_0111 - eif_0011,
    eif_type1_pse1 = eif_1111 - eif_0111) %>%
  summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                            se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator", "type", "estimand", "measure")) %>%
  dplyr::select(-type) %>% 
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(estimand = factor(estimand,
                           levels = rev(c("ate", "pse4", "pse3", "pse2", "pse1")),
                           labels = rev(c(expression(paste("Total Effect (", psi[`1111`]-psi[`0000`], ")")),
                                          expression(paste("Direct Effect (", psi[`0001`]-psi[`0000`], ")")),
                                          expression(paste("via G3 Identity (", psi[`0011`]-psi[`0001`], ")")),
                                          expression(paste("via G2 Identity (", psi[`0111`]-psi[`0011`], ")")),
                                          expression(paste("via G1 Identity (", psi[`1111`]-psi[`0111`], ")")))))) 

#------------------------------#
#   Non-parametric bootstrap   #
#------------------------------#

B <- 2000

boots <- matrix(NA, nrow = B, ncol = 5)

for (b in 1:B){
  
  if (b %% 100 == 0){
    cat("bootstrap sample ", b, "\n")
  }
  
  dfi <- tatar %>% sample_frac(replace = TRUE)
  
  dfi_p0 <- model.matrix(a0_form, data = dfi)[, -1] %>% as_tibble()
  dfi_p1 <- model.matrix(a1_form, data = dfi)[, -1] %>% as_tibble()
  dfi_p2 <- model.matrix(a2_form, data = dfi)[, -1] %>% as_tibble()
  dfi_p3 <- model.matrix(a3_form, data = dfi)[, -1] %>% as_tibble()
  
  dfi_mu0 <- model.matrix(y0_form, data = dfi)[, -1] %>% as_tibble()
  dfi_mu1 <- model.matrix(y1_form, data = dfi)[, -1] %>% as_tibble()
  dfi_mu2 <- model.matrix(y2_form, data = dfi)[, -1] %>% as_tibble()
  dfi_mu3 <- model.matrix(y3_form, data = dfi)[, -1] %>% as_tibble()
  
  dfi_mu3n <- model.matrix(y3_form, data = mutate(dfi, violence = 0))[, -1] %>% as_tibble()
  dfi_mu2n <- model.matrix(y2_form, data = mutate(dfi, violence = 0))[, -1] %>% as_tibble()
  dfi_mu1n <- model.matrix(y1_form, data = mutate(dfi, violence = 0))[, -1] %>% as_tibble()
  dfi_mu0n <- model.matrix(y0_form, data = mutate(dfi, violence = 0))[, -1] %>% as_tibble()
  
  dfi_mu3y <- model.matrix(y3_form, data = mutate(dfi, violence = 1))[, -1] %>% as_tibble()
  dfi_mu2y <- model.matrix(y2_form, data = mutate(dfi, violence = 1))[, -1] %>% as_tibble()
  dfi_mu1y <- model.matrix(y1_form, data = mutate(dfi, violence = 1))[, -1] %>% as_tibble()
  dfi_mu0y <- model.matrix(y0_form, data = mutate(dfi, violence = 1))[, -1] %>% as_tibble()


  p0_glm <- glm(a0_form, family = binomial("logit"), data = dfi)
  p1_glm <- glm(a1_form, family = binomial("logit"), data = dfi)
  p2_glm <- glm(a2_form, family = binomial("logit"), data = dfi)
  p3_glm <- glm(a3_form, family = binomial("logit"), data = dfi)
  
  #--------------------------#
  #      Treatment Model     #
  #--------------------------#
  
  dfi <- dfi %>% mutate(
    p0_fit = p0_glm$fitted.values,
    p1_fit = p1_glm$fitted.values,
    p2_fit = p2_glm$fitted.values,
    p3_fit = p3_glm$fitted.values,
    
    w0_n = I(violence == 0)/(1 - p0_fit),
    w0_y = I(violence == 1)/p0_fit,
    w1_nn = I(violence == 0)/(1 - p0_fit) * 1,
    w1_ny = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit,
    w1_yn = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit),
    w1_yy = I(violence == 1)/p0_fit * 1,
    w2_000 = I(violence == 0)/(1 - p0_fit) * 1 * 1,
    w2_001 = I(violence == 1)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit,
    w2_010 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit),
    w2_011 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1,
    w2_100 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * 1,
    w2_101 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit,
    w2_110 = I(violence == 0)/p0_fit * 1 * p2_fit/(1 - p2_fit),
    w2_111 = I(violence == 1)/p0_fit * 1 * 1,
    
    w3_0000 = I(violence == 0)/(1 - p0_fit) * 1 * 1 * 1,
    w3_0001 = I(violence == 1)/(1 - p0_fit) * 1 * 1 * (1 - p3_fit)/p3_fit,
    w3_0010 = I(violence == 0)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit * p3_fit/(1 - p3_fit),
    w3_0011 = I(violence == 1)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit * 1,
    
    w3_0100 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit) * 1,
    w3_0101 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit) * (1 - p3_fit)/p3_fit,
    w3_0110 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1 * p3_fit/(1 - p3_fit),
    w3_0111 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1 * 1,
    
    w3_1000 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * 1 * 1,
    w3_1001 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * 1 * (1 - p3_fit)/p3_fit,
    w3_1010 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit * p3_fit/(1 - p3_fit),
    w3_1011 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit * 1,
    
    w3_1100 = I(violence == 0)/p0_fit * 1 * p2_fit/(1 - p2_fit) * 1,
    w3_1101 = I(violence == 1)/p0_fit * 1 * p2_fit/(1 - p2_fit) * (1 - p3_fit)/p3_fit,
    w3_1110 = I(violence == 0)/p0_fit * 1 * 1 * p3_fit/(1 - p3_fit),
    w3_1111 = I(violence == 1)/p0_fit * 1 * 1 * 1,
  )

  #-----------------------#
  #     Outcome Model     #
  #-----------------------#
  
  mu3_glm <- glm(y3_form, family = binomial("logit"), data = dfi)
  
  dfi$mu3_fit_a4n <- predict(mu3_glm, type = "response", newdata = dfi_mu3n)
  dfi$mu3_fit_a4y <- predict(mu3_glm, type = "response", newdata = dfi_mu3y)
  
  mu3_fit_a4y_form <- mu3_fit_a4n_form <- y2_form
  mu3_fit_a4n_form[[2]] <- expr(mu3_fit_a4n)
  mu3_fit_a4y_form[[2]] <- expr(mu3_fit_a4y)
  
  mu2_glm_a4n <- glm(mu3_fit_a4n_form, family = quasibinomial("logit"), data = dfi)
  mu2_glm_a4y <- glm(mu3_fit_a4y_form, family = quasibinomial("logit"), data = dfi)
  
  dfi$mu2_fit_a4n_a3n <- predict(mu2_glm_a4n, type = "response", newdata = dfi_mu2n)
  dfi$mu2_fit_a4n_a3y <- predict(mu2_glm_a4n, type = "response", newdata = dfi_mu2y)
  dfi$mu2_fit_a4y_a3n <- predict(mu2_glm_a4y, type = "response", newdata = dfi_mu2n)
  dfi$mu2_fit_a4y_a3y <- predict(mu2_glm_a4y, type = "response", newdata = dfi_mu2y)
  
  mu2_fit_a4y_a3y_form <- mu2_fit_a4y_a3n_form <- y1_form
  mu2_fit_a4n_a3y_form <- mu2_fit_a4n_a3n_form <- y1_form
  
  mu2_fit_a4y_a3y_form[[2]] <- expr(mu2_fit_a4y_a3y)
  mu2_fit_a4y_a3n_form[[2]] <- expr(mu2_fit_a4y_a3n)
  mu2_fit_a4n_a3y_form[[2]] <- expr(mu2_fit_a4n_a3y)
  mu2_fit_a4n_a3n_form[[2]] <- expr(mu2_fit_a4n_a3n)
  
  mu1_glm_a4n_a3n <- glm(mu2_fit_a4n_a3n_form, family = quasibinomial("logit"), data = dfi)
  mu1_glm_a4n_a3y <- glm(mu2_fit_a4n_a3y_form, family = quasibinomial("logit"), data = dfi)
  mu1_glm_a4y_a3n <- glm(mu2_fit_a4y_a3n_form, family = quasibinomial("logit"), data = dfi)
  mu1_glm_a4y_a3y <- glm(mu2_fit_a4y_a3y_form, family = quasibinomial("logit"), data = dfi)
  
  dfi$mu1_fit_a4n_a3n_a2n <- predict(mu1_glm_a4n_a3n, type = "response", newdata = dfi_mu1n)
  dfi$mu1_fit_a4n_a3n_a2y <- predict(mu1_glm_a4n_a3n, type = "response", newdata = dfi_mu1y)
  
  dfi$mu1_fit_a4n_a3y_a2n <- predict(mu1_glm_a4n_a3y, type = "response", newdata = dfi_mu1n)
  dfi$mu1_fit_a4n_a3y_a2y <- predict(mu1_glm_a4n_a3y, type = "response", newdata = dfi_mu1y)
  
  dfi$mu1_fit_a4y_a3n_a2n <- predict(mu1_glm_a4y_a3n, type = "response", newdata = dfi_mu1n)
  dfi$mu1_fit_a4y_a3n_a2y <- predict(mu1_glm_a4y_a3n, type = "response", newdata = dfi_mu1y)
  
  dfi$mu1_fit_a4y_a3y_a2n <- predict(mu1_glm_a4y_a3y, type = "response", newdata = dfi_mu1n)
  dfi$mu1_fit_a4y_a3y_a2y <- predict(mu1_glm_a4y_a3y, type = "response", newdata = dfi_mu1y)
  
  mu1_fit_a4y_a3y_a2y_form <- mu1_fit_a4y_a3n_a2y_form <- y0_form
  mu1_fit_a4n_a3y_a2y_form <- mu1_fit_a4n_a3n_a2y_form <- y0_form
  mu1_fit_a4y_a3y_a2n_form <- mu1_fit_a4y_a3n_a2n_form <- y0_form
  mu1_fit_a4n_a3y_a2n_form <- mu1_fit_a4n_a3n_a2n_form <- y0_form
  
  mu1_fit_a4y_a3y_a2y_form[[2]] <- expr(mu1_fit_a4y_a3y_a2y)
  mu1_fit_a4y_a3y_a2n_form[[2]] <- expr(mu1_fit_a4y_a3y_a2n)
  mu1_fit_a4y_a3n_a2y_form[[2]] <- expr(mu1_fit_a4y_a3n_a2y)
  mu1_fit_a4y_a3n_a2n_form[[2]] <- expr(mu1_fit_a4y_a3n_a2n)
  
  mu1_fit_a4n_a3y_a2y_form[[2]] <- expr(mu1_fit_a4n_a3y_a2y)
  mu1_fit_a4n_a3y_a2n_form[[2]] <- expr(mu1_fit_a4n_a3y_a2n)
  mu1_fit_a4n_a3n_a2y_form[[2]] <- expr(mu1_fit_a4n_a3n_a2y)
  mu1_fit_a4n_a3n_a2n_form[[2]] <- expr(mu1_fit_a4n_a3n_a2n)
  
  mu0_glm_a4n_a3n_a2n <- glm(mu1_fit_a4n_a3n_a2n_form, family = quasibinomial("logit"), data = dfi)
  mu0_glm_a4n_a3n_a2y <- glm(mu1_fit_a4n_a3n_a2y_form, family = quasibinomial("logit"), data = dfi)
  mu0_glm_a4n_a3y_a2n <- glm(mu1_fit_a4n_a3y_a2n_form, family = quasibinomial("logit"), data = dfi)
  mu0_glm_a4n_a3y_a2y <- glm(mu1_fit_a4n_a3y_a2y_form, family = quasibinomial("logit"), data = dfi)
  
  mu0_glm_a4y_a3n_a2n <- glm(mu1_fit_a4y_a3n_a2n_form, family = quasibinomial("logit"), data = dfi)
  mu0_glm_a4y_a3n_a2y <- glm(mu1_fit_a4y_a3n_a2y_form, family = quasibinomial("logit"), data = dfi)
  mu0_glm_a4y_a3y_a2n <- glm(mu1_fit_a4y_a3y_a2n_form, family = quasibinomial("logit"), data = dfi)
  mu0_glm_a4y_a3y_a2y <- glm(mu1_fit_a4y_a3y_a2y_form, family = quasibinomial("logit"), data = dfi)
  
  dfi$mu0_fit_a4n_a3n_a2n_a1n <- predict(mu0_glm_a4n_a3n_a2n, type = "response", newdata = dfi_mu0n)
  dfi$mu0_fit_a4n_a3n_a2n_a1y <- predict(mu0_glm_a4n_a3n_a2n, type = "response", newdata = dfi_mu0y)
  
  dfi$mu0_fit_a4n_a3n_a2y_a1n <- predict(mu0_glm_a4n_a3n_a2y, type = "response", newdata = dfi_mu0n)
  dfi$mu0_fit_a4n_a3n_a2y_a1y <- predict(mu0_glm_a4n_a3n_a2y, type = "response", newdata = dfi_mu0y)
  
  dfi$mu0_fit_a4n_a3y_a2n_a1n <- predict(mu0_glm_a4n_a3y_a2n, type = "response", newdata = dfi_mu0n)
  dfi$mu0_fit_a4n_a3y_a2n_a1y <- predict(mu0_glm_a4n_a3y_a2n, type = "response", newdata = dfi_mu0y)
  
  dfi$mu0_fit_a4n_a3y_a2y_a1n <- predict(mu0_glm_a4n_a3y_a2y, type = "response", newdata = dfi_mu0n)
  dfi$mu0_fit_a4n_a3y_a2y_a1y <- predict(mu0_glm_a4n_a3y_a2y, type = "response", newdata = dfi_mu0y)
  
  dfi$mu0_fit_a4y_a3n_a2n_a1n <- predict(mu0_glm_a4y_a3n_a2n, type = "response", newdata = dfi_mu0n)
  dfi$mu0_fit_a4y_a3n_a2n_a1y <- predict(mu0_glm_a4y_a3n_a2n, type = "response", newdata = dfi_mu0y)
  
  dfi$mu0_fit_a4y_a3n_a2y_a1n <- predict(mu0_glm_a4y_a3n_a2y, type = "response", newdata = dfi_mu0n)
  dfi$mu0_fit_a4y_a3n_a2y_a1y <- predict(mu0_glm_a4y_a3n_a2y, type = "response", newdata = dfi_mu0y)
  
  dfi$mu0_fit_a4y_a3y_a2n_a1n <- predict(mu0_glm_a4y_a3y_a2n, type = "response", newdata = dfi_mu0n)
  dfi$mu0_fit_a4y_a3y_a2n_a1y <- predict(mu0_glm_a4y_a3y_a2n, type = "response", newdata = dfi_mu0y)
  
  dfi$mu0_fit_a4y_a3y_a2y_a1n <- predict(mu0_glm_a4y_a3y_a2y, type = "response", newdata = dfi_mu0n)
  dfi$mu0_fit_a4y_a3y_a2y_a1y <- predict(mu0_glm_a4y_a3y_a2y, type = "response", newdata = dfi_mu0y)
  
  for (s in 1:S){
    
    a1 <- estimands$a1[[s]]
    a2 <- estimands$a2[[s]]
    a3 <- estimands$a3[[s]]
    a4 <- estimands$a4[[s]]
    
    dfi <- dfi %>%
      mutate(
        
        wt0_deno = a1 * p0_fit + (1 - a1) * (1 - p0_fit),
        wt1_nume = a1 * p1_fit + (1 - a1) * (1 - p1_fit),
        wt1_deno = a2 * p1_fit + (1 - a2) * (1 - p1_fit),
        wt2_nume = a2 * p2_fit + (1 - a2) * (1 - p2_fit),
        wt2_deno = a3 * p2_fit + (1 - a3) * (1 - p2_fit),
        wt3_nume = a3 * p3_fit + (1 - a3) * (1 - p3_fit),
        wt3_deno = a4 * p3_fit + (1 - a4) * (1 - p3_fit),
        
        !!sym(paste0("w0_", a1, a2, a3, a4)) := as.double(violence==a1) / wt0_deno,
        !!sym(paste0("w1_", a1, a2, a3, a4)) := as.double(violence==a2) * wt1_nume/wt1_deno/wt0_deno,
        !!sym(paste0("w2_", a1, a2, a3, a4)) := as.double(violence==a3) * wt2_nume/wt2_deno * wt1_nume/wt1_deno/wt0_deno,
        !!sym(paste0("w3_", a1, a2, a3, a4)) := as.double(violence==a4) * wt3_nume/wt3_deno * wt2_nume/wt2_deno * wt1_nume/wt1_deno/wt0_deno
      )
    
    dfi[dfi$violence == a1, paste0("w0_", a1, a2, a3, a4)] <- trimQ(dfi[dfi$violence == a1, paste0("w0_", a1, a2, a3, a4)])
    dfi[dfi$violence == a2, paste0("w1_", a1, a2, a3, a4)] <- trimQ(dfi[dfi$violence == a2, paste0("w1_", a1, a2, a3, a4)])
    dfi[dfi$violence == a3, paste0("w2_", a1, a2, a3, a4)] <- trimQ(dfi[dfi$violence == a3, paste0("w2_", a1, a2, a3, a4)])
    dfi[dfi$violence == a4, paste0("w3_", a1, a2, a3, a4)] <- trimQ(dfi[dfi$violence == a4, paste0("w3_", a1, a2, a3, a4)])
    
    dfi <- dfi %>%
      mutate(
        
        !!sym(paste0("mu3fit_", a1, a2, a3, a4)) := a4 * mu3_fit_a4y + (1 - a4) * mu3_fit_a4n,
        
        !!sym(paste0("mu2fit_", a1, a2, a3, a4)) := a4 * a3 * mu2_fit_a4y_a3y + a4 * (1 - a3) * mu2_fit_a4y_a3n +
          (1 - a4) * a3 * mu2_fit_a4n_a3y + (1 - a4) * (1 - a3) * mu2_fit_a4n_a3n,
        
        !!sym(paste0("mu1fit_", a1, a2, a3, a4)) := a4 * a3 * a2 * mu1_fit_a4y_a3y_a2y + a4 * a3 * (1 - a2) * mu1_fit_a4y_a3y_a2n +
          a4 * (1 - a3) * a2 * mu1_fit_a4y_a3n_a2y + a4 * (1 - a3) * (1 - a2) * mu1_fit_a4y_a3n_a2n +
          (1 - a4) * a3 * a2 * mu1_fit_a4n_a3y_a2y + (1 - a4) * a3 * (1 - a2) * mu1_fit_a4n_a3y_a2n +
          (1 - a4) * (1 - a3) * a2 * mu1_fit_a4n_a3n_a2y + (1 - a4) * (1 - a3) * (1 - a2) * mu1_fit_a4n_a3n_a2n,
        
        !!sym(paste0("mu0fit_", a1, a2, a3, a4)) := a4 * a3 * a2 * a1 * mu0_fit_a4y_a3y_a2y_a1y + a4 * a3 * (1 - a2) * a1 * mu0_fit_a4y_a3y_a2n_a1y +
          a4 * (1 - a3) * a2 * a1 * mu0_fit_a4y_a3n_a2y_a1y + a4 * (1 - a3) * (1 - a2) * a1 * mu0_fit_a4y_a3n_a2n_a1y +
          (1 - a4) * a3 * a2 * a1 * mu0_fit_a4n_a3y_a2y_a1y + (1 - a4) * a3 * (1 - a2) * a1 * mu0_fit_a4n_a3y_a2n_a1y +
          (1 - a4) * (1 - a3) * a2 * a1 * mu0_fit_a4n_a3n_a2y_a1y + (1 - a4) * (1 - a3) * (1 - a2) * a1 * mu0_fit_a4n_a3n_a2n_a1y +
          a4 * a3 * a2 * (1 - a1) * mu0_fit_a4y_a3y_a2y_a1n + a4 * a3 * (1 - a2) * (1 - a1) * mu0_fit_a4y_a3y_a2n_a1n +
          a4 * (1 - a3) * a2 * (1 - a1) * mu0_fit_a4y_a3n_a2y_a1n + a4 * (1 - a3) * (1 - a2) * (1 - a1) * mu0_fit_a4y_a3n_a2n_a1n +
          (1 - a4) * a3 * a2 * (1 - a1) * mu0_fit_a4n_a3y_a2y_a1n + (1 - a4) * a3 * (1 - a2) * (1 - a1) * mu0_fit_a4n_a3y_a2n_a1n +
          (1 - a4) * (1 - a3) * a2 * (1 - a1) * mu0_fit_a4n_a3n_a2y_a1n + (1 - a4) * (1 - a3) * (1 - a2) * (1 - a1) * mu0_fit_a4n_a3n_a2n_a1n,
        
        !!sym(paste0("www_", a1, a2, a3, a4)) := !!sym(paste0("w3_", a1, a2, a3, a4)) * annex,
        !!sym(paste0("iii_", a1, a2, a3, a4)) := !!sym(paste0("mu0fit_", a1, a2, a3, a4)),
        !!sym(paste0("eif_", a1, a2, a3, a4)) := !!sym(paste0("w3_", a1, a2, a3, a4)) *
          (annex - !!sym(paste0("mu3fit_", a1, a2, a3, a4))) +
          !!sym(paste0("w2_", a1, a2, a3, a4)) * (!!sym(paste0("mu3fit_", a1, a2, a3, a4)) -
                                                    !!sym(paste0("mu2fit_", a1, a2, a3, a4))) +
          !!sym(paste0("w1_", a1, a2, a3, a4)) * (!!sym(paste0("mu2fit_", a1, a2, a3, a4)) -
                                                    !!sym(paste0("mu1fit_", a1, a2, a3, a4))) +
          !!sym(paste0("w0_", a1, a2, a3, a4)) * (!!sym(paste0("mu1fit_", a1, a2, a3, a4)) -
                                                    !!sym(paste0("mu0fit_", a1, a2, a3, a4))) +
          !!sym(paste0("mu0fit_", a1, a2, a3, a4))
      )
  }
  
  out_dfi <- dfi %>%
    mutate(
      eif_type1_ate = eif_1111 - eif_0000,
      eif_type1_pse4 = eif_0001 - eif_0000,
      eif_type1_pse3 = eif_0011 - eif_0001,
      eif_type1_pse2 = eif_0111 - eif_0011,
      eif_type1_pse1 = eif_1111 - eif_0111) %>%
    summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                              se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
    pivot_longer(everything()) %>%
    separate(name, into = c("estimator", "type", "estimand", "measure")) %>%
    dplyr::select(-type) %>% 
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(estimand = factor(estimand,
                             levels = rev(c("ate", "pse4", "pse3", "pse2", "pse1")),
                             labels = rev(c(expression(paste("Total Effect (", psi[`1111`]-psi[`0000`], ")")),
                                            expression(paste("Direct Effect (", psi[`0001`]-psi[`0000`], ")")),
                                            expression(paste("via G3 Identity (", psi[`0011`]-psi[`0001`], ")")),
                                            expression(paste("via G2 Identity (", psi[`0111`]-psi[`0011`], ")")),
                                            expression(paste("via G1 Identity (", psi[`1111`]-psi[`0111`], ")")))))) 
  
  boots[b, ] <- out_dfi$est
}

out_df <- out_df %>% 
  mutate(
    se_naive = se,
    se = apply(boots, 2, sd),
    lower = apply(boots, 2, quantile, 0.025),
    upper = apply(boots, 2, quantile, 0.975)
  )

table6_5par <-  out_df %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se) %>% 
  mutate_at(c("est", "se", "lower", "upper"), ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")")) %>% 
  mutate(out = paste(est, intv)) %>% 
  dplyr::select(-lower, -upper, -intv) 

print(table6_5par)

#==============================#
#        DML Estimation        #
#==============================#

#----------------------------#
#  DML cross-fitting setup   #
#----------------------------#

K <- 5

# create cross-fitting split
cf_fold <- createFolds(df[[y]], K)

main_list <- vector(mode = "list", K)

k <- 1

for(k in 1:K){
  
  cat(" cross-fitting fold ", k, "\n")

  #design matrices for the different model
  
  # auxiliary and main data
  aux <- df[-cf_fold[[k]], ]
  main <- df[cf_fold[[k]], ]
  
  aux_p0 <- df_p0[-cf_fold[[k]], ]
  aux_p1 <- df_p1[-cf_fold[[k]], ]
  aux_p2 <- df_p2[-cf_fold[[k]], ]
  aux_p3 <- df_p3[-cf_fold[[k]], ]
  
  main_p0 <- df_p0[cf_fold[[k]], ]
  main_p1 <- df_p1[cf_fold[[k]], ]
  main_p2 <- df_p2[cf_fold[[k]], ]
  main_p3 <- df_p3[cf_fold[[k]], ]
  
  aux_mu0 <- df_mu0[-cf_fold[[k]], ]
  aux_mu1 <- df_mu1[-cf_fold[[k]], ]
  aux_mu2 <- df_mu2[-cf_fold[[k]], ]
  aux_mu3 <- df_mu3[-cf_fold[[k]], ]
  
  main_mu0 <- df_mu0[cf_fold[[k]], ]
  main_mu1 <- df_mu1[cf_fold[[k]], ]
  main_mu2 <- df_mu2[cf_fold[[k]], ]
  main_mu3 <- df_mu3[cf_fold[[k]], ]

  #--------------------------#
  #      Treatment Model     #
  #--------------------------#
  
  p0_sl <- SuperLearner(
    Y          = aux$violence,
    X          = aux_p0,
    newX       = df_p0,
    family     = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  p1_sl <- SuperLearner(
    Y          = aux$violence,
    X          = aux_p1,
    newX       = df_p1,
    family     = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  p2_sl <- SuperLearner(
    Y          = aux$violence,
    X          = aux_p2,
    newX       = df_p2,
    family     = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE, trimLogit = 0.001),
    cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
  )
  
  p3_sl <- SuperLearner(
    Y          = aux$violence,
    X          = aux_p3,
    newX       = df_p3,
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
    p3_fit = p3_sl$SL.predict,
    w0_n = I(violence == 0)/(1 - p0_fit),
    w0_y = I(violence == 1)/p0_fit,
    w1_nn = I(violence == 0)/(1 - p0_fit) * 1,
    w1_ny = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit,
    w1_yn = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit),
    w1_yy = I(violence == 1)/p0_fit * 1,
    w2_000 = I(violence == 0)/(1 - p0_fit) * 1 * 1,
    w2_001 = I(violence == 1)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit,
    w2_010 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit),
    w2_011 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1,
    w2_100 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * 1,
    w2_101 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit,
    w2_110 = I(violence == 0)/p0_fit * 1 * p2_fit/(1 - p2_fit),
    w2_111 = I(violence == 1)/p0_fit * 1 * 1,
    
    w3_0000 = I(violence == 0)/(1 - p0_fit) * 1 * 1 * 1,
    w3_0001 = I(violence == 1)/(1 - p0_fit) * 1 * 1 * (1 - p3_fit)/p3_fit,
    w3_0010 = I(violence == 0)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit * p3_fit/(1 - p3_fit),
    w3_0011 = I(violence == 1)/(1 - p0_fit) * 1 * (1 - p2_fit)/p2_fit * 1,
    
    w3_0100 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit) * 1,
    w3_0101 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * p2_fit/(1 - p2_fit) * (1 - p3_fit)/p3_fit,
    w3_0110 = I(violence == 0)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1 * p3_fit/(1 - p3_fit),
    w3_0111 = I(violence == 1)/(1 - p0_fit) * (1 - p1_fit)/p1_fit * 1 * 1,
    
    w3_1000 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * 1 * 1,
    w3_1001 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * 1 * (1 - p3_fit)/p3_fit,
    w3_1010 = I(violence == 0)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit * p3_fit/(1 - p3_fit),
    w3_1011 = I(violence == 1)/p0_fit * p1_fit/(1 - p1_fit) * (1 - p2_fit)/p2_fit * 1,
    
    w3_1100 = I(violence == 0)/p0_fit * 1 * p2_fit/(1 - p2_fit) * 1,
    w3_1101 = I(violence == 1)/p0_fit * 1 * p2_fit/(1 - p2_fit) * (1 - p3_fit)/p3_fit,
    w3_1110 = I(violence == 0)/p0_fit * 1 * 1 * p3_fit/(1 - p3_fit),
    w3_1111 = I(violence == 1)/p0_fit * 1 * 1 * 1,
  )
  
  #-----------------------#
  #     Outcome Model     #
  #-----------------------#

  mu3_sl <- SuperLearner(
    Y          = aux$annex,
    X          = aux_mu3,
    family     = binomial(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu3_fit_a4n <- predict.SuperLearner(mu3_sl, newdata = df_mu3n)$pred
  df$mu3_fit_a4y <- predict.SuperLearner(mu3_sl, newdata = df_mu3y)$pred
  
  mu2_sl_a4n <- SuperLearner(
    Y          = df$mu3_fit_a4n[-cf_fold[[k]]],
    X          = aux_mu2,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu2_sl_a4y <- SuperLearner(
    Y          = df$mu3_fit_a4y[-cf_fold[[k]]],
    X          = aux_mu2,
    family     = gaussian(),
    # obsWeights = aux$weight,
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu2_fit_a4n_a3n <- predict.SuperLearner(mu2_sl_a4n, newdata = df_mu2n)$pred
  df$mu2_fit_a4n_a3y <- predict.SuperLearner(mu2_sl_a4n, newdata = df_mu2y)$pred
  df$mu2_fit_a4y_a3n <- predict.SuperLearner(mu2_sl_a4y, newdata = df_mu2n)$pred
  df$mu2_fit_a4y_a3y <- predict.SuperLearner(mu2_sl_a4y, newdata = df_mu2y)$pred
  
  mu1_sl_a4n_a3n <- SuperLearner(
    Y          = df$mu2_fit_a4n_a3n[-cf_fold[[k]]],
    X          = aux_mu1,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu1_sl_a4n_a3y <- SuperLearner(
    Y          = df$mu2_fit_a4n_a3y[-cf_fold[[k]]],
    X          = aux_mu1,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu1_sl_a4y_a3n <- SuperLearner(
    Y          = df$mu2_fit_a4y_a3n[-cf_fold[[k]]],
    X          = aux_mu1,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu1_sl_a4y_a3y <- SuperLearner(
    Y          = df$mu2_fit_a4y_a3y[-cf_fold[[k]]],
    X          = aux_mu1,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu1_fit_a4n_a3n_a2n <- predict.SuperLearner(mu1_sl_a4n_a3n, newdata = df_mu1n)$pred
  df$mu1_fit_a4n_a3n_a2y <- predict.SuperLearner(mu1_sl_a4n_a3n, newdata = df_mu1y)$pred
  
  df$mu1_fit_a4n_a3y_a2n <- predict.SuperLearner(mu1_sl_a4n_a3y, newdata = df_mu1n)$pred
  df$mu1_fit_a4n_a3y_a2y <- predict.SuperLearner(mu1_sl_a4n_a3y, newdata = df_mu1y)$pred
  
  df$mu1_fit_a4y_a3n_a2n <- predict.SuperLearner(mu1_sl_a4y_a3n, newdata = df_mu1n)$pred
  df$mu1_fit_a4y_a3n_a2y <- predict.SuperLearner(mu1_sl_a4y_a3n, newdata = df_mu1y)$pred
  
  df$mu1_fit_a4y_a3y_a2n <- predict.SuperLearner(mu1_sl_a4y_a3y, newdata = df_mu1n)$pred
  df$mu1_fit_a4y_a3y_a2y <- predict.SuperLearner(mu1_sl_a4y_a3y, newdata = df_mu1y)$pred
  
  mu0_sl_a4n_a3n_a2n <- SuperLearner(
    Y          = df$mu1_fit_a4n_a3n_a2n[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a4n_a3n_a2y <- SuperLearner(
    Y          = df$mu1_fit_a4n_a3n_a2y[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a4n_a3y_a2n <- SuperLearner(
    Y          = df$mu1_fit_a4n_a3y_a2n[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a4n_a3y_a2y <- SuperLearner(
    Y          = df$mu1_fit_a4n_a3y_a2y[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a4y_a3n_a2n <- SuperLearner(
    Y          = df$mu1_fit_a4y_a3n_a2n[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a4y_a3n_a2y <- SuperLearner(
    Y          = df$mu1_fit_a4y_a3n_a2y[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a4y_a3y_a2n <- SuperLearner(
    Y          = df$mu1_fit_a4y_a3y_a2n[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  mu0_sl_a4y_a3y_a2y <- SuperLearner(
    Y          = df$mu1_fit_a4y_a3y_a2y[-cf_fold[[k]]],
    X          = aux_mu0,
    family     = gaussian(),
    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    control    = list(saveFitLibrary = TRUE),
    cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
  )
  
  df$mu0_fit_a4n_a3n_a2n_a1n <- predict.SuperLearner(mu0_sl_a4n_a3n_a2n, newdata = df_mu0n)$pred
  df$mu0_fit_a4n_a3n_a2n_a1y <- predict.SuperLearner(mu0_sl_a4n_a3n_a2n, newdata = df_mu0y)$pred
  
  df$mu0_fit_a4n_a3n_a2y_a1n <- predict.SuperLearner(mu0_sl_a4n_a3n_a2y, newdata = df_mu0n)$pred
  df$mu0_fit_a4n_a3n_a2y_a1y <- predict.SuperLearner(mu0_sl_a4n_a3n_a2y, newdata = df_mu0y)$pred
  
  df$mu0_fit_a4n_a3y_a2n_a1n <- predict.SuperLearner(mu0_sl_a4n_a3y_a2n, newdata = df_mu0n)$pred
  df$mu0_fit_a4n_a3y_a2n_a1y <- predict.SuperLearner(mu0_sl_a4n_a3y_a2n, newdata = df_mu0y)$pred
  
  df$mu0_fit_a4n_a3y_a2y_a1n <- predict.SuperLearner(mu0_sl_a4n_a3y_a2y, newdata = df_mu0n)$pred
  df$mu0_fit_a4n_a3y_a2y_a1y <- predict.SuperLearner(mu0_sl_a4n_a3y_a2y, newdata = df_mu0y)$pred
  
  df$mu0_fit_a4y_a3n_a2n_a1n <- predict.SuperLearner(mu0_sl_a4y_a3n_a2n, newdata = df_mu0n)$pred
  df$mu0_fit_a4y_a3n_a2n_a1y <- predict.SuperLearner(mu0_sl_a4y_a3n_a2n, newdata = df_mu0y)$pred
  
  df$mu0_fit_a4y_a3n_a2y_a1n <- predict.SuperLearner(mu0_sl_a4y_a3n_a2y, newdata = df_mu0n)$pred
  df$mu0_fit_a4y_a3n_a2y_a1y <- predict.SuperLearner(mu0_sl_a4y_a3n_a2y, newdata = df_mu0y)$pred
  
  df$mu0_fit_a4y_a3y_a2n_a1n <- predict.SuperLearner(mu0_sl_a4y_a3y_a2n, newdata = df_mu0n)$pred
  df$mu0_fit_a4y_a3y_a2n_a1y <- predict.SuperLearner(mu0_sl_a4y_a3y_a2n, newdata = df_mu0y)$pred
  
  df$mu0_fit_a4y_a3y_a2y_a1n <- predict.SuperLearner(mu0_sl_a4y_a3y_a2y, newdata = df_mu0n)$pred
  df$mu0_fit_a4y_a3y_a2y_a1y <- predict.SuperLearner(mu0_sl_a4y_a3y_a2y, newdata = df_mu0y)$pred
  
  main_list[[k]] <- df[cf_fold[[k]], ]
}

main_df <- reduce(main_list, bind_rows)

for (s in 1:S){
  
  a1 <- estimands$a1[[s]]
  a2 <- estimands$a2[[s]]
  a3 <- estimands$a3[[s]]
  a4 <- estimands$a4[[s]]
  
  main_df <- main_df %>%
    mutate(
      
      wt0_deno = a1 * p0_fit + (1 - a1) * (1 - p0_fit),
      wt1_nume = a1 * p1_fit + (1 - a1) * (1 - p1_fit),
      wt1_deno = a2 * p1_fit + (1 - a2) * (1 - p1_fit),
      wt2_nume = a2 * p2_fit + (1 - a2) * (1 - p2_fit),
      wt2_deno = a3 * p2_fit + (1 - a3) * (1 - p2_fit),
      wt3_nume = a3 * p3_fit + (1 - a3) * (1 - p3_fit),
      wt3_deno = a4 * p3_fit + (1 - a4) * (1 - p3_fit),
      
      !!sym(paste0("w0_", a1, a2, a3, a4)) := as.double(violence==a1) / wt0_deno,
      !!sym(paste0("w1_", a1, a2, a3, a4)) := as.double(violence==a2) * wt1_nume/wt1_deno/wt0_deno,
      !!sym(paste0("w2_", a1, a2, a3, a4)) := as.double(violence==a3) * wt2_nume/wt2_deno * wt1_nume/wt1_deno/wt0_deno,
      !!sym(paste0("w3_", a1, a2, a3, a4)) := as.double(violence==a4) * wt3_nume/wt3_deno * wt2_nume/wt2_deno * wt1_nume/wt1_deno/wt0_deno
    )
  
  main_df[main_df$violence == a1, paste0("w0_", a1, a2, a3, a4)] <- trimQ(main_df[main_df$violence == a1, paste0("w0_", a1, a2, a3, a4)])
  main_df[main_df$violence == a2, paste0("w1_", a1, a2, a3, a4)] <- trimQ(main_df[main_df$violence == a2, paste0("w1_", a1, a2, a3, a4)])
  main_df[main_df$violence == a3, paste0("w2_", a1, a2, a3, a4)] <- trimQ(main_df[main_df$violence == a3, paste0("w2_", a1, a2, a3, a4)])
  main_df[main_df$violence == a4, paste0("w3_", a1, a2, a3, a4)] <- trimQ(main_df[main_df$violence == a4, paste0("w3_", a1, a2, a3, a4)])
  
  main_df <- main_df %>%
    mutate(
      
      !!sym(paste0("mu3fit_", a1, a2, a3, a4)) := a4 * mu3_fit_a4y + (1 - a4) * mu3_fit_a4n,
      
      !!sym(paste0("mu2fit_", a1, a2, a3, a4)) := a4 * a3 * mu2_fit_a4y_a3y + a4 * (1 - a3) * mu2_fit_a4y_a3n +
        (1 - a4) * a3 * mu2_fit_a4n_a3y + (1 - a4) * (1 - a3) * mu2_fit_a4n_a3n,
      
      !!sym(paste0("mu1fit_", a1, a2, a3, a4)) := a4 * a3 * a2 * mu1_fit_a4y_a3y_a2y + a4 * a3 * (1 - a2) * mu1_fit_a4y_a3y_a2n +
        a4 * (1 - a3) * a2 * mu1_fit_a4y_a3n_a2y + a4 * (1 - a3) * (1 - a2) * mu1_fit_a4y_a3n_a2n +
        (1 - a4) * a3 * a2 * mu1_fit_a4n_a3y_a2y + (1 - a4) * a3 * (1 - a2) * mu1_fit_a4n_a3y_a2n +
        (1 - a4) * (1 - a3) * a2 * mu1_fit_a4n_a3n_a2y + (1 - a4) * (1 - a3) * (1 - a2) * mu1_fit_a4n_a3n_a2n,
      
      !!sym(paste0("mu0fit_", a1, a2, a3, a4)) := a4 * a3 * a2 * a1 * mu0_fit_a4y_a3y_a2y_a1y + a4 * a3 * (1 - a2) * a1 * mu0_fit_a4y_a3y_a2n_a1y +
        a4 * (1 - a3) * a2 * a1 * mu0_fit_a4y_a3n_a2y_a1y + a4 * (1 - a3) * (1 - a2) * a1 * mu0_fit_a4y_a3n_a2n_a1y +
        (1 - a4) * a3 * a2 * a1 * mu0_fit_a4n_a3y_a2y_a1y + (1 - a4) * a3 * (1 - a2) * a1 * mu0_fit_a4n_a3y_a2n_a1y +
        (1 - a4) * (1 - a3) * a2 * a1 * mu0_fit_a4n_a3n_a2y_a1y + (1 - a4) * (1 - a3) * (1 - a2) * a1 * mu0_fit_a4n_a3n_a2n_a1y +
        a4 * a3 * a2 * (1 - a1) * mu0_fit_a4y_a3y_a2y_a1n + a4 * a3 * (1 - a2) * (1 - a1) * mu0_fit_a4y_a3y_a2n_a1n +
        a4 * (1 - a3) * a2 * (1 - a1) * mu0_fit_a4y_a3n_a2y_a1n + a4 * (1 - a3) * (1 - a2) * (1 - a1) * mu0_fit_a4y_a3n_a2n_a1n +
        (1 - a4) * a3 * a2 * (1 - a1) * mu0_fit_a4n_a3y_a2y_a1n + (1 - a4) * a3 * (1 - a2) * (1 - a1) * mu0_fit_a4n_a3y_a2n_a1n +
        (1 - a4) * (1 - a3) * a2 * (1 - a1) * mu0_fit_a4n_a3n_a2y_a1n + (1 - a4) * (1 - a3) * (1 - a2) * (1 - a1) * mu0_fit_a4n_a3n_a2n_a1n,
      
      !!sym(paste0("www_", a1, a2, a3, a4)) := !!sym(paste0("w3_", a1, a2, a3, a4)) * annex,
      !!sym(paste0("iii_", a1, a2, a3, a4)) := !!sym(paste0("mu0fit_", a1, a2, a3, a4)),
      !!sym(paste0("eif_", a1, a2, a3, a4)) := !!sym(paste0("w3_", a1, a2, a3, a4)) *
        (annex - !!sym(paste0("mu3fit_", a1, a2, a3, a4))) +
        !!sym(paste0("w2_", a1, a2, a3, a4)) * (!!sym(paste0("mu3fit_", a1, a2, a3, a4)) -
                                                  !!sym(paste0("mu2fit_", a1, a2, a3, a4))) +
        !!sym(paste0("w1_", a1, a2, a3, a4)) * (!!sym(paste0("mu2fit_", a1, a2, a3, a4)) -
                                                  !!sym(paste0("mu1fit_", a1, a2, a3, a4))) +
        !!sym(paste0("w0_", a1, a2, a3, a4)) * (!!sym(paste0("mu1fit_", a1, a2, a3, a4)) -
                                                  !!sym(paste0("mu0fit_", a1, a2, a3, a4))) +
        !!sym(paste0("mu0fit_", a1, a2, a3, a4))
    )
}

out_df <- main_df %>%
  mutate(eif_type1_ate = eif_1111 - eif_0000,
         eif_type1_pse4 = eif_0001 - eif_0000,
         eif_type1_pse3 = eif_0011 - eif_0001,
         eif_type1_pse2 = eif_0111 - eif_0011,
         eif_type1_pse1 = eif_1111 - eif_0111) %>%
  summarise_at(vars(contains("type")), list(est = ~ wtd.mean(.x),
                                            se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator", "type", "estimand", "measure")) %>%
  dplyr::select(-type) %>% 
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(estimand = factor(estimand,
                           levels = rev(c("ate", "pse4", "pse3", "pse2", "pse1")),
                           labels = rev(c(expression(paste("Total Effect (", psi[`1111`]-psi[`0000`], ")")),
                                          expression(paste("Direct Effect (", psi[`0001`]-psi[`0000`], ")")),
                                          expression(paste("via G3 Identity (", psi[`0011`]-psi[`0001`], ")")),
                                          expression(paste("via G2 Identity (", psi[`0111`]-psi[`0011`], ")")),
                                          expression(paste("via G1 Identity (", psi[`1111`]-psi[`0111`], ")")))))) 

table6_5np <-  out_df %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se) %>% 
  mutate_at(c("est", "se", "lower", "upper"), ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")")) %>% 
  mutate(out = paste(est, intv)) %>% 
  dplyr::select(-lower, -upper, -intv) 

print(table6_5np)

#-------------------------------#
#        COMBINE RESULTS        #
#-------------------------------#

.norm65 <- function(s) {
  s <- as.character(s)
  dplyr::case_when(
    grepl("Total|ATE",           s, ignore.case = TRUE) ~ "ATE(1,0)",
    grepl("Direct",              s, ignore.case = TRUE) ~ "PSE_D→Y (1,0)",
    grepl("G3|M3",               s, ignore.case = TRUE) ~ "PSE_D→M3→Y (1,0)",
    grepl("G2|M2",               s, ignore.case = TRUE) ~ "PSE_D→M2↝Y (1,0)",
    grepl("G1|M1",               s, ignore.case = TRUE) ~ "PSE_D→M1↝Y (1,0)",
    TRUE ~ s
  )
}

par_tbl <- table6_5par %>%
  dplyr::transmute(
    Estimand = .norm65(estimand),
    `Parametric MR` = out
  )

dml_tbl <- table6_5np %>%
  dplyr::transmute(
    Estimand = .norm65(estimand),
    DML = out
  )

# enforce order
order_65 <- c("ATE(1,0)",
              "PSE_D→Y (1,0)",
              "PSE_D→M3→Y (1,0)",
              "PSE_D→M2↝Y (1,0)",
              "PSE_D→M1↝Y (1,0)")

final_rst <- tibble::tibble(Estimand = order_65) %>%
  dplyr::left_join(par_tbl, by = "Estimand") %>%
  dplyr::left_join(dml_tbl, by = "Estimand")

# Open log
sink(log_path, split = TRUE)

width_curr <- getOption("width"); options(width = 300)
print(final_rst, right = FALSE, row.names = FALSE)
options(width = width_curr)

# Close log
sink()
