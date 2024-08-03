###Table 5.3###

rm(list=ls())

packages<-c("tidyverse", "margins", "mediation", "foreach", "doParallel", "doRNG")

for (package.i in packages) {
  suppressPackageStartupMessages(library(package.i, character.only=TRUE))
}
source("utils.R")

# devtools::install_github("xiangzhou09/paths")
library(paths)

##office
datadir <- "../../data/" 
logdir <- "../../code/ch5/_LOGS/"

#sink(paste(logdir, "table_5-1_log.txt", sep=""))

##input data
load("Brader_et_al2008.RData")

# function for demeaning
demean <- function(x) x - mean(x, na.rm = TRUE)

# data preprocessing
Brader2 <- Brader %>%
  dplyr::select(immigr, emo, p_harm, tone_eth, ppage, ppeducat, ppgender, ppincimp) %>% na.omit() %>%
  mutate(immigr = as.numeric(scale(4 - immigr)),
         hs = (ppeducat == "high school"),
         sc = (ppeducat == "some college"),
         ba = (ppeducat == "bachelor's degree or higher"),
         female = (ppgender == "female")) %>%
  mutate_at(vars(emo, p_harm, ppage, female, hs, sc, ba, ppincimp), demean)

summary(Brader2)

d <- "tone_eth"
m1 <- "p_harm"
m2 <- "emo"
y <- "immigr"
x <- c("ppage", "female", "hs", "sc", "ba", "ppincimp")

df <- Brader2
m12 <- c(m1, m2)

set.seed(02138)

##########################################
## RI without interactions
##########################################

mediators <- list(m1, m2)

formula_m0 <- as.formula(paste(y, "~", paste(c(x, d), collapse = " + ")))
formula_m1 <- as.formula(paste(y, "~", paste(c(x, d, m1), collapse = " + ")))
formula_m2 <- as.formula(paste(y, "~", paste(c(x, d, m1, m2), collapse = " + ")))

# outcome models
glm_m0 <- glm(formula_m0, data = df)
glm_m1 <- glm(formula_m1, data = df)
glm_m2 <- glm(formula_m2, data = df)
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

paths_glm <- paths(a = d, y = y, m = mediators,
                   glm_ymodels, data = df, nboot = 2000)

pure <- paths_glm$pure %>% 
  filter(decomposition == "Type I") %>% 
  mutate(estimand = factor(estimand, levels = c("total", "direct", "via M2", "via M1"),
                           labels = c("ATE", "AY", "AM2Y", "AM1Y"))) %>% 
  rename(est = estimate) %>% 
  dplyr::select(-decomposition, -se, -p) %>% 
  mutate_at(c("est", "lower", "upper"),  ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")"))
  
  
hybrid <- paths_glm$hybrid %>% 
  filter(decomposition == "Type I") %>% 
  mutate(estimand = factor(estimand, levels = c("total", "direct", "via M2", "via M1"),
                           labels = c("ATE", "AY", "AM2Y", "AM1Y"))) %>% 
  rename(est = estimate) %>% 
  dplyr::select(-decomposition, -se, -p) %>% 
  mutate_at(c("est", "lower", "upper"),  ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")"))

m_out2a <- bind_rows(pure, hybrid) %>% 
  mutate(estimator = factor(estimator, levels = c("pure", "hybrid"))) %>% 
  arrange(estimator, estimand) %>% 
  mutate(out = paste(est, intv)) %>% 
  add_column(interactions = "n")


##########################################
## RI with D-M interactions
##########################################

formula_m0 <- as.formula(paste(y, "~", paste(c(x, d), collapse = " + ")))
formula_m1x <- as.formula(paste(y, "~", paste(c(x, d, m1, paste(d, m1, sep = "*")), collapse = " + ")))
formula_m2x <- as.formula(paste(y, "~", paste(c(x, d, m12, paste(d, m12, sep = "*")), collapse = " + ")))

# outcome models
glm_m0 <- glm(formula_m0, data = df)
glm_m1x <- glm(formula_m1x, data = df)
glm_m2x <- glm(formula_m2x, data = df)
glm_ymodels_x <- list(glm_m0, glm_m1x, glm_m2x)

paths_glm_x <- paths(a = d, y = y, m = mediators,
                     glm_ymodels_x, data = df, nboot = 2000)

pure <- paths_glm_x$pure %>% 
  filter(decomposition == "Type I") %>% 
  mutate(estimand = factor(estimand, levels = c("total", "direct", "via M2", "via M1"),
                           labels = c("ATE", "AY", "AM2Y", "AM1Y"))) %>% 
  rename(est = estimate) %>% 
  dplyr::select(-decomposition, -se, -p) %>% 
  mutate_at(c("est", "lower", "upper"),  ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")"))


hybrid <- paths_glm_x$hybrid %>% 
  filter(decomposition == "Type I") %>% 
  mutate(estimand = factor(estimand, levels = c("total", "direct", "via M2", "via M1"),
                           labels = c("ATE", "AY", "AM2Y", "AM1Y"))) %>% 
  rename(est = estimate) %>% 
  dplyr::select(-decomposition, -se, -p) %>% 
  mutate_at(c("est", "lower", "upper"),  ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")"))

m_out2b <- bind_rows(pure, hybrid) %>% 
  mutate(estimator = factor(estimator, levels = c("pure", "hybrid"))) %>% 
  arrange(estimator, estimand) %>% 
  mutate(out = paste(est, intv)) %>% 
  add_column(interactions = "x")

##########################################
## RI with D-M and D-C interactions
##########################################

formula_m0xx <- as.formula(paste(y, "~", paste(c(x, d, paste(d, x, sep = "*")), collapse = " + ")))
formula_m1xx <- as.formula(paste(y, "~", paste(c(x, d, m1, paste(d, x, sep = "*"),
                                                 paste(d, m1, sep = "*")), collapse = " + ")))
formula_m2xx <- as.formula(paste(y, "~", paste(c(x, d, m12, paste(d, x, sep = "*"),
                                                 paste(d, m12, sep = "*")), collapse = " + ")))

# outcome models
glm_m0xx <- glm(formula_m0xx, data = df)
glm_m1xx <- glm(formula_m1xx, data = df)
glm_m2xx <- glm(formula_m2xx, data = df)
glm_ymodels_xx <- list(glm_m0xx, glm_m1xx, glm_m2xx)

paths_glm_xx <- paths(a = d, y = y, m = mediators,
                     glm_ymodels_xx,  data = df, nboot = 2000)

pure <- paths_glm_xx$pure %>% 
  filter(decomposition == "Type I") %>% 
  mutate(estimand = factor(estimand, levels = c("total", "direct", "via M2", "via M1"),
                           labels = c("ATE", "AY", "AM2Y", "AM1Y"))) %>% 
  rename(est = estimate) %>% 
  dplyr::select(-decomposition, -se, -p) %>% 
  mutate_at(c("est", "lower", "upper"),  ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")"))


hybrid <- paths_glm_xx$hybrid %>% 
  filter(decomposition == "Type I") %>% 
  mutate(estimand = factor(estimand, levels = c("total", "direct", "via M2", "via M1"),
                           labels = c("ATE", "AY", "AM2Y", "AM1Y"))) %>% 
  rename(est = estimate) %>% 
  dplyr::select(-decomposition, -se, -p) %>% 
  mutate_at(c("est", "lower", "upper"),  ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")"))

m_out2c <- bind_rows(pure, hybrid) %>% 
  mutate(estimator = factor(estimator, levels = c("pure", "hybrid"))) %>% 
  arrange(estimator, estimand) %>% 
  mutate(out = paste(est, intv)) %>% 
  add_column(interactions = "xx")

m_out2 <- bind_rows(m_out2a, m_out2b, m_out2c) %>% 
  filter(estimator == "pure")

save.image(file = "table5-8.RData")

write_csv(m_out2, file = "table5-8.csv")

