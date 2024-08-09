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
tatar <- readRDS(paste(datadir, "Tatar/tatar.rds", sep=""))

# variable names
x <- c("kulak", "prosoviet_pre", "religiosity_pre", "land_pre",
       "orchard_pre", "animals_pre", "carriage_pre", "otherprop_pre")
a <- "violence"
y <- "annex"
m1 <- c("trust_g1", "victim_g1", "fear_g1")
m2 <- c("trust_g2", "victim_g2", "fear_g2")
m3 <- c("trust_g3", "victim_g3", "fear_g3")
m <- list(m1, m2, m3)
df <- tatar
n <- nrow(df)

##########################################################
# Formulas for treatment and outcome models
##########################################################

a0_form <- as.formula(paste(a, " ~ ", paste(x, collapse= "+")))
a1_form <- as.formula(paste(a, " ~ ", paste(c(x, m1), collapse= "+")))
a2_form <- as.formula(paste(a, " ~ ", paste(c(x, m1, m2), collapse= "+")))
a3_form <- as.formula(paste(a, " ~ ", paste(c(x, m1, m2, m3), collapse= "+")))

y0_form <- as.formula(paste(y, " ~ ", paste(c(x, a), collapse= "+")))
y1_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1), collapse= "+")))
y2_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1, m2), collapse= "+")))
y3_form <- as.formula(paste(y, " ~ ", paste(c(x, a, m1, m2, m3), collapse= "+")))

##########################################################
# Main analyses
##########################################################

estimands <- expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1)) %>%
  `colnames<-`(c("a1", "a2", "a3", "a4"))

S <- nrow(estimands)

#################################################
# Design matrices for different models
#################################################

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

#################################################
# Treatment Models
#################################################

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


#################################################
# Outcome models
#################################################

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


#################################################
# Nonparametric Bootstrap
#################################################

B <- 2000

boots <- matrix(NA, nrow = B, ncol = 5)

for (b in 1:B){
  
  if (b %% 10 ==0){
    cat("bootstrap sample ", b, "\n")
  }
  
  dfi <- tatar %>% sample_frac(replace = TRUE)
  
  #################################################
  # Design matrices for different models
  #################################################
  
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
  
  #################################################
  # Treatment Models
  #################################################
  
  p0_glm <- glm(a0_form, family = binomial("logit"), data = dfi)
  p1_glm <- glm(a1_form, family = binomial("logit"), data = dfi)
  p2_glm <- glm(a2_form, family = binomial("logit"), data = dfi)
  p3_glm <- glm(a3_form, family = binomial("logit"), data = dfi)
  
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
  
  
  #################################################
  # Outcome models
  #################################################
  
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
  
#################################################
# Plot
#################################################

# set my ggplot theme
mytheme <- theme_minimal(base_size = 18) + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(color = "grey30"))
theme_set(mytheme)

ggplot(out_df, aes(x = estimand, y = est, shape = estimator)) +
  geom_pointrange(aes(ymin = lower,  ymax = upper),
                  position = position_dodge(width = - 0.5), size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_shape("", labels = parse_format()) +
  scale_color_discrete("", labels = parse_format()) +
  scale_x_discrete("", labels = parse_format()) +
  scale_y_continuous("Effects of Ancestor Victimization on Regime Support") +
  coord_flip()

table6_5par <-  out_df %>%
  mutate(lower = est - 1.96 * se, upper = est + 1.96 * se) %>% 
  mutate_at(c("est", "se", "lower", "upper"), ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")")) %>% 
  mutate(out = paste(est, intv)) %>% 
  dplyr::select(-lower, -upper, -intv) 

write_csv(table6_5par, file = "table6-5par.csv")

save.image(file = "table6-5par.RData")
