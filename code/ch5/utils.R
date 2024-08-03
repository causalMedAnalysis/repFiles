demean <- function(x, w = rep(1, length(x))) x - weighted.mean(x, w, na.rm = TRUE)

trim <- function(x, min = 0.01, max = 1) {
  x[x<min] <- min
  x[x>max] <- max
  x
}

trimQ <- function(x, low = 0.01, high = 0.99) {
  
  min <- quantile(x, low)
  max <- quantile(x, high)
  
  x[x<min] <- min
  x[x>max] <- max
  x
}

# # estimator based on linear models with no D-M interaction
# linmed <- function(df, d, m, y, x){
#   
#   m_form <- paste(m, "~", paste(c(x, d), collapse = " + "))
#   y_form <- paste(y, "~", paste(c(x, d, m), collapse = " + "))
#   
#   m_model <- lm(m_form, data = df)
#   y_model <- lm(y_form, data = df)
#   
#   NDE <- CDE <- y_model$coef[[d]]
#   NIE <- m_model$coef[[d]] * y_model$coef[[m]]
#   ATE <- NDE + NIE
#   
#   point_est <- list(ATE = ATE, NDE = NDE, NIE = NIE)
#   
#   point_est
# }
# 
# # estimator based on linear models with D-M interaction
# linmedx <- function(df, d, m, y, x){
#   
#   # demean x in df
#   df_x <- df[, x, drop = FALSE]
#   for(i in seq_along(x)) df[[x[i]]] <- demean(df_x[[i]])
#   
#   m_form <- paste(m, "~", paste(c(x, d), collapse = " + "))
#   y_form <- paste(y, "~", paste(c(x, d, m, paste(d, "*", m)), collapse = " + "))
#   
#   m_model <- lm(m_form, data = df)
#   y_model <- lm(y_form, data = df)
#   
#   NDE <- y_model$coef[[d]] + y_model$coef[[paste0(d, ":", m)]] * m_model$coef[["(Intercept)"]] # gamma_2 + gamma_4 * beta_0
#   NIE <- m_model$coef[[d]] * (y_model$coef[[m]] + y_model$coef[[paste0(d, ":", m)]]) # beta2 * (gamma_3 + gamma_4)
#   ATE <- NDE + NIE 
#   # CDE0 <- y_model$coef[[d]] # gamma_2
#   
#   point_est <- list(ATE = ATE, NDE = NDE, NIE = NIE)
#   
#   point_est
# }

# estimator based on linear models with no D-M interaction (m can be multivariate)
linmed <- function(df, d, m, y, x){
  
  m_forms <- purrr::map(m, ~ paste(.x, "~", paste(c(x, d), collapse = " + "))) 
  y_form <- paste(y, "~", paste(c(x, d, m), collapse = " + "))
  
  m_models <- purrr::map(m_forms, ~ lm(.x, data = df))
  y_model <- lm(y_form, data = df)
  
  NDE <- CDE <- y_model$coef[[d]]
  
  am_coefs <- purrr::map_dbl(m_models, function(x) x$coef[[d]])
  NIE <- sum(am_coefs * y_model$coef[m])
  ATE <- NDE + NIE
  
  point_est <- list(ATE = ATE, NDE = NDE, NIE = NIE)
  
  point_est
}


# estimator based on linear models with D-M interaction (m can be multivariate)
linmedx <- function(df, d, m, y, x){
  
  # demean x in df
  df_x <- df[, x, drop = FALSE]
  for(i in seq_along(x)) df[[x[i]]] <- demean(df_x[[i]])
  
  m_forms <- purrr::map(m, ~ paste(.x, "~", paste(c(x, d), collapse = " + "))) 
  y_form <- paste(y, "~", paste(c(x, d, m, paste(d, "*", m)), collapse = " + "))
  
  m_models <- purrr::map(m_forms, ~ lm(.x, data = df))
  y_model <- lm(y_form, data = df)
  
  NDE_adj <- purrr::map2_dbl(m, m_models, function(a, b) y_model$coef[[paste0(d, ":", a)]] * b$coef[["(Intercept)"]])
  NIE_adj <- purrr::map2_dbl(m, m_models, function(a, b) b$coef[[d]] * (y_model$coef[[a]] + y_model$coef[[paste0(d, ":", a)]]))
  
  NDE <- y_model$coef[[d]] + sum(NDE_adj) # gamma_2 + gamma_4 * beta_0
  NIE <- sum(NIE_adj) # beta2 * (gamma_3 + gamma_4)
  ATE <- NDE + NIE 
  # CDE0 <- y_model$coef[[d]] # gamma_2
  
  point_est <- list(ATE = ATE, NDE = NDE, NIE = NIE)
  
  point_est
}

# weighting estimator (m can be multivariate)
ipwmed <- function(df, d, m, y, x){
  
  d_form1 <- paste(d, "~", paste(c(x), collapse = " + "))
  d_form2 <- paste(d, "~", paste(c(x, m), collapse = " + "))
  
  d_model1 <- glm(d_form1, data = df, family=binomial("logit"))
  d_model2 <- glm(d_form2, data = df, family=binomial("logit"))
  
  d_ps1 <- predict(d_model1, type = "response") 
  d_ps2 <- predict(d_model2, type = "response")
  
  sw1 <- (1-mean(d_ps1))/(1-d_ps1)
  sw2 <- mean(d_ps1)/d_ps1
  sw3 <- (1-d_ps2)/d_ps2 * mean(d_ps1)/(1-d_ps1)
  
  # truncate weights among obs with nonzero weights. 
  sw1[df[[d]] == 0] <- trimQ(sw1[df[[d]] == 0])
  sw2[df[[d]] == 1] <- trimQ(sw2[df[[d]] == 1])
  sw3[df[[d]] == 1] <- trimQ(sw3[df[[d]] == 1])
  
  Ehat_Y0M0 <- Hmisc::wtd.mean(df[[y]], sw1 * (df[[d]] == 0))
  Ehat_Y1M1 <- Hmisc::wtd.mean(df[[y]], sw2 * (df[[d]] == 1))
  Ehat_Y1M0 <- Hmisc::wtd.mean(df[[y]], sw3 * (df[[d]] == 1))
  
  ATE <- Ehat_Y1M1 - Ehat_Y0M0
  NDE <- Ehat_Y1M0 - Ehat_Y0M0
  NIE <- Ehat_Y1M1 - Ehat_Y1M0
  
  point_est <- list(ATE = ATE, NDE = NDE, NIE = NIE)
  
  point_est
}

# simulation estimator based on the mediation package (m can't be multivariate)
simmed <- function(df, d, m, y, x,
                   m_type = c("continuous", "binary"),
                   y_type = c("continuous", "binary"), sims = 1000){
  
  m_form <- paste(m, "~", paste(c(x, d), collapse = " + "))
  y_form <- paste(y, "~", paste(c(x, d, m), collapse = " + "))
  
  m_type <- match.arg(m_type)
  y_type <- match.arg(y_type)
  
  m_model <- if(m_type == "binary") glm(m_form, data = df, family = binomial("logit")) else lm(m_form, data = df)
  y_model <- if(y_type == "binary") glm(y_form, data = df, family = binomial("logit")) else lm(y_form, data = df)
  
  sim_est <- mediation::mediate(m_model, y_model, treat = d, mediator = m, sims = sims)
  
  point_est <- list(ATE = sim_est$tau.coef,
                    NDE = sim_est$z0,
                    NIE = sim_est$d1)
  
  point_est
}