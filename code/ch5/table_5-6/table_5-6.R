###Table 5.5-5.6###

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
nlsy_raw <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

d <- "att22"
m1 <- "ever_unemp_age3539"
m2 <- "log_faminc_adj_age3539"
y <- "std_cesd_age40"
x <- c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")

nlsy <- nlsy_raw[complete.cases(nlsy_raw[, c(d, "cesd_age40", m1, m2, x)]),] %>% 
  mutate(std_cesd_age40 = as.numeric(scale(cesd_age40)))

summary(nlsy)

df <- nlsy
m12 <- c(m1, m2)

set.seed(02138)

##########################################
## RI without interactions
##########################################

mediators <- list(m1, m2)

formula_m0 <- as.formula(paste(y, "~", paste(c(x, d), collapse = " + ")))
formula_m1 <- as.formula(paste(y, "~", paste(c(x, d, m1), collapse = " + ")))
formula_m2 <- as.formula(paste(y, "~", paste(c(x, d, m1, m2), collapse = " + ")))
formula_ps <- as.formula(paste(d, "~", paste(x, collapse = " + ")))

# outcome models
glm_m0 <- glm(formula_m0, data = df)
glm_m1 <- glm(formula_m1, data = df)
glm_m2 <- glm(formula_m2, data = df)
glm_ymodels <- list(glm_m0, glm_m1, glm_m2)

# propensity score model
glm_ps <- glm(formula_ps, family = binomial("logit"), data = df)

paths_glm <- paths(a = d, y = y, m = mediators,
                   glm_ymodels, ps_model = glm_ps, 
                   data = df, nboot = 2000)

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
## RI with interactions
##########################################

formula_m0 <- as.formula(paste(y, "~", paste(c(x, d), collapse = " + ")))
formula_m1x <- as.formula(paste(y, "~", paste(c(x, d, m1, paste(d, m1, sep = "*")), collapse = " + ")))
formula_m2x <- as.formula(paste(y, "~", paste(c(x, d, m12, paste(d, m12, sep = "*")), collapse = " + ")))
formula_ps <- as.formula(paste(d, "~", paste(x, collapse = " + ")))

# outcome models
glm_m0 <- glm(formula_m0, data = df)
glm_m1x <- glm(formula_m1x, data = df)
glm_m2x <- glm(formula_m2x, data = df)
glm_ymodels_x <- list(glm_m0, glm_m1x, glm_m2x)

# propensity score model
glm_ps <- glm(formula_ps, family = binomial("logit"), data = df)

paths_glm_x <- paths(a = d, y = y, m = mediators,
                     glm_ymodels_x, ps_model = glm_ps, 
                     data = df, nboot = 2000)


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
  add_column(interactions = "y")

m_out2 <- bind_rows(m_out2a, m_out2b)

write_csv(m_out2, file = "table5-6.csv")
