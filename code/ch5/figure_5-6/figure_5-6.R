###Figure 5.6###

rm(list=ls())

packages<-c("tidyverse", "margins", "mediation", "foreach", "doParallel", "doRNG")

for (package.i in packages) {
  suppressPackageStartupMessages(library(package.i, character.only=TRUE))
}
source("utils.R")

# devtools::install_github("xiangzhou09/paths")
library(paths)
library(latex2exp)

##office
datadir <- "../../data/" 
logdir <- "../../code/ch5/_LOGS/"

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
## Sensitivity Analysis
##########################################

paths_sens <- sens(paths_glm, confounded = "M1", estimand = "direct",
                   gamma_values = seq(0, 2, 0.005),
                   eta_values = seq(-0.5, 0.5, 0.005))

# gamma2 and eta2 for selected covariates

summary(Brader)

Brader2_ed <- Brader2 %>%
  mutate(female2 = as.double(ppgender == "female"),
         ba2 = as.double(ppeducat == "bachelor's degree or higher"))

# formulas
d <- "tone_eth"
m1 <- "p_harm"
m2 <- "emo"
y <- "immigr"
x <- c("ppage", "female2", "hs", "sc", "ba2", "ppincimp")

form_y <- as.formula(paste0(y, "~", paste0(c(x, d, m1, m2), collapse = "+")))
form_female2 <- as.formula(paste0("female2 ~", paste0(c(setdiff(x, "female"), d, m1, m2),
                                                    collapse = "+")))
form_ba2 <- as.formula(paste0("ba2 ~", paste0(c(setdiff(x, "ba"), d, m1, m2),
                                                  collapse = "+"))) 

mod_y <- lm(form_y, data = Brader2_ed)
mod_female2 <- lm(form_female2, Brader2_ed)
mod_ba2 <- lm(form_ba2, Brader2_ed)

gamma2_female2 <- mod_y$coefficients["female2"]
gamma2_ba2 <- mod_y$coefficients["ba2"]
eta2_female2 <- mod_female2$coefficients[d]
eta2_ba2 <- mod_ba2$coefficients[d]

xlab_new <- bquote("Conditional Difference in the Prevalence of U across levels of D"~(delta[DU]))
ylab_new <- bquote("Conditional Difference in the Mean of Y across levels of U"~(delta[UY]))


# sensitivity analysis plot
plot(paths_sens, outcome_name = "Support for Immigration") +
  xlab(xlab_new) +
  ylab(ylab_new) +
  annotate("point", x = eta2_female2, y = gamma2_female2) +
  annotate("text", x = eta2_female2 + 0.05, y = gamma2_female2, label = "female", size = 5) +
  annotate("point", x = eta2_ba2, y = gamma2_ba2) +
  annotate("text", x = eta2_ba2 - 0.05, y = gamma2_ba2, label = "college\n graduate", size = 5) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c(NA, "grey70"),  na.value = NA)

ggsave("figure_5-6.png", width = 10, height = 7)
