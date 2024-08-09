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
y0_form <- as.formula(paste(y, " ~ ", paste(c(x, a), collapse= "+")))

##########################################################
# Main analyses
##########################################################

estimands <- expand.grid(c(0, 1)) %>%
  `colnames<-`(c("a"))

S <- nrow(estimands)

#################################################
# Design matrices for different models
#################################################

df_p0 <- model.matrix(a0_form, data = df)[, -1] %>% as_tibble()
df_mu0 <- model.matrix(y0_form, data = df)[, -1] %>% as_tibble()

df_mu0n <- model.matrix(y0_form, data = mutate(df, att22 = 0))[, -1] %>% as_tibble()
df_mu0y <- model.matrix(y0_form, data = mutate(df, att22 = 1))[, -1] %>% as_tibble()

#################################################
# Treatment Models
#################################################

p0_glm <- glm(a0_form, family = binomial("logit"), data = df)

df <- df %>% mutate(
  p0_fit = p0_glm$fitted.values,
  w0_n = I(att22 == 0)/(1 - p0_fit),
  w0_y = I(att22 == 1)/p0_fit,
)


#################################################
# Outcome models
#################################################

mu0_glm <- lm(y0_form, data = df)

df$mu0_fit_n <- predict(mu0_glm, newdata = df_mu0n)
df$mu0_fit_y <- predict(mu0_glm, newdata = df_mu0y)

for (s in 1:S){
  
  a <- estimands$a[[s]]
  
  df <- df %>%
    mutate(
      
      wt0_deno = a * p0_fit + (1 - a) * (1 - p0_fit),
      
      !!sym(paste0("w0_", a)) := as.double(att22==a)/wt0_deno,
      
    )
  
  df[df$att22 == a, paste0("w0_", a)] <- trimQ(df[df$att22 == a, paste0("w0_", a)])
  
  df <- df %>% 
    mutate(
      
      !!sym(paste0("mu0fit_", a)) := a * mu0_fit_y + (1 - a) * mu0_fit_n,
      
      !!sym(paste0("www_", a)) := !!sym(paste0("w0_", a)) * std_cesd_age40,
      !!sym(paste0("iii_", a)) := !!sym(paste0("mu0fit_", a)),
      !!sym(paste0("eif_", a)) := !!sym(paste0("w0_", a)) *
        (std_cesd_age40 - !!sym(paste0("mu0fit_", a))) +
        !!sym(paste0("mu0fit_", a))
      
    )
}

out_df <- df %>%
  mutate(www_ate = www_1 - www_0,
         iii_ate = iii_1 - iii_0,
         eif_ate = eif_1 - eif_0) %>%
  summarise_at(vars(contains(c("iii_", "www_", "eif_"))), list(est = ~ wtd.mean(.x),
                                            se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("estimator","estimand", "measure")) %>%
  pivot_wider(names_from = measure, values_from = value)

#################################################
# Nonparametric Bootstrap
#################################################

B <- 2000

boots <- matrix(NA, nrow = B, ncol = 9)

for (b in 1:B){
  
  if (b %% 10 ==0){
    cat("bootstrap sample ", b, "\n")
  }
  
  dfi <- nlsy %>% sample_frac(replace = TRUE)
  
  #################################################
  # Design matrices for different models
  #################################################
  
  df_p0 <- model.matrix(a0_form, data = dfi)[, -1] %>% as_tibble()
  df_mu0 <- model.matrix(y0_form, data = dfi)[, -1] %>% as_tibble()
  
  df_mu0n <- model.matrix(y0_form, data = mutate(dfi, att22 = 0))[, -1] %>% as_tibble()
  df_mu0y <- model.matrix(y0_form, data = mutate(dfi, att22 = 1))[, -1] %>% as_tibble()
  
  #################################################
  # Treatment Models
  #################################################
  
  p0_glm <- glm(a0_form, family = binomial("logit"), data = dfi)
  
  dfi <- dfi %>% mutate(
    p0_fit = p0_glm$fitted.values,
    w0_n = I(att22 == 0)/(1 - p0_fit),
    w0_y = I(att22 == 1)/p0_fit,
  )
  
  
  #################################################
  # Outcome models
  #################################################
  
  mu0_glm <- lm(y0_form, data = dfi)
  
  dfi$mu0_fit_n <- predict(mu0_glm, newdata = df_mu0n)
  dfi$mu0_fit_y <- predict(mu0_glm, newdata = df_mu0y)
  
  for (s in 1:S){
    
    a <- estimands$a[[s]]
    
    dfi <- dfi %>%
      mutate(
        
        wt0_deno = a * p0_fit + (1 - a) * (1 - p0_fit),
        
        !!sym(paste0("w0_", a)) := as.double(att22==a)/wt0_deno,
        
      )
    
    dfi[dfi$att22 == a, paste0("w0_", a)] <- trimQ(dfi[dfi$att22 == a, paste0("w0_", a)])
    
    dfi <- dfi %>% 
      mutate(
        
        !!sym(paste0("mu0fit_", a)) := a * mu0_fit_y + (1 - a) * mu0_fit_n,
        
        !!sym(paste0("www_", a)) := !!sym(paste0("w0_", a)) * std_cesd_age40,
        !!sym(paste0("iii_", a)) := !!sym(paste0("mu0fit_", a)),
        !!sym(paste0("eif_", a)) := !!sym(paste0("w0_", a)) *
          (std_cesd_age40 - !!sym(paste0("mu0fit_", a))) +
          !!sym(paste0("mu0fit_", a))
        
      )
  }
  
  out_dfi <- dfi %>%
    mutate(www_ate = www_1 - www_0,
           iii_ate = iii_1 - iii_0,
           eif_ate = eif_1 - eif_0) %>%
    summarise_at(vars(contains(c("iii_", "www_", "eif_"))), list(est = ~ wtd.mean(.x),
                                              se = ~ sqrt(wtd.var(.x)/length(.x)))) %>%
    pivot_longer(everything()) %>%
    separate(name, into = c("estimator","estimand", "measure")) %>%
    pivot_wider(names_from = measure, values_from = value)
  
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
  scale_y_continuous("Effects of College Attendance on Depression") +
  coord_flip()

table6_1par <-  out_df %>%
  dplyr::select(-se_naive) %>% 
  mutate_at(c("est", "se", "lower", "upper"), ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", lower, ", ", upper, ")")) %>% 
  mutate(out = paste(est, intv)) %>% 
  dplyr::select(-lower, -upper, -intv) 

write_csv(table6_1par, file = "table6-1par.csv")

save.image(file = "table6-1par.RData")


