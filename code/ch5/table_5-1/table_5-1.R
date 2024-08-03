###Table 5.1###

rm(list=ls())

packages<-c("tidyverse", "margins", "mediation", "foreach", "doParallel", "doRNG")

for (package.i in packages) {
  suppressPackageStartupMessages(library(package.i, character.only=TRUE))
}
source("utils.R")

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

nlsy4m1 <- nlsy_raw[complete.cases(nlsy_raw[, c(d, "cesd_age40", m1, x)]),] %>% 
  mutate(std_cesd_age40 = as.numeric(scale(cesd_age40)))

nlsy4m2 <- nlsy_raw[complete.cases(nlsy_raw[, c(d, "cesd_age40", m2, x)]),] %>% 
  mutate(std_cesd_age40 = as.numeric(scale(cesd_age40)))

## linear model estimators

#linear model w/o exposure-mediator interaction

linmed_m1 <- linmed(nlsy4m1, d, m1, y, x)
linmed_m2 <- linmed(nlsy4m2, d, m2, y, x)

#linear model w exposure-mediator interaction

linmedx_m1 <- linmedx(nlsy4m1, d, m1, y, x)
linmedx_m2 <- linmedx(nlsy4m2, d, m2, y, x) 

# # simulation estimator 
# 
# simmed_m1 <- simmed(nlsy4m1, d, m1, y, x, m_type = "binary")
# 
# simmed_m2 <- simmed(nlsy4m2, d, m2, y, x, m_type = "continuous")

# ipw estimator

ipwmed_m1 <- ipwmed(nlsy4m1, d, m1, y, x)

ipwmed_m2 <- ipwmed(nlsy4m2, d, m2, y, x)

# summary of point estimates

est_m1 <- list(linmed_m1, linmedx_m1, ipwmed_m1) %>%
  map(unlist) %>%
  bind_rows() %>%
  t() %>% 
  `colnames<-`( c("^lma", "^lmi", "^ipw"))

est_m2 <- list(linmed_m2, linmedx_m2, ipwmed_m2) %>%
  map(unlist) %>%
  bind_rows() %>%
  t() %>% 
  `colnames<-`( c("^lma", "^lmi", "^ipw"))

# bootstrap SEs for linmed and linmedx and ipwmed

#compute bootstrap estimates for M1

#setup parallel computing cluster
ncores <- parallel::detectCores()-1
my.cluster <- parallel::makeCluster(ncores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster,list("linmed", "linmedx", "ipwmed"), envir=environment())
registerDoRNG(3308004)

nboot <- 2000

m1_boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {
  
  boot_data <- nlsy4m1[sample(nrow(nlsy4m1), nrow(nlsy4m1), replace=TRUE),]
  
  boot_linmed <- linmed(boot_data, d, m1, y, x)
  
  boot_linmedx <- linmedx(boot_data, d, m1, y, x)
  
  boot_ipwmed <- ipwmed(boot_data, d, m1, y, x)
  
  list(linmed = boot_linmed,
       linmedx = boot_linmedx,
       ipwmed = boot_ipwmed)
}

stopCluster(my.cluster)
rm(my.cluster)

clnms <- expand_grid(method = c("linmed", "linmedx", "ipwmed"), estimand = c("ATE", "NDE", "NIE")) %>% unite("") %>% unlist()

m1_boot_mat <- m1_boot %>% 
  map(unlist) %>% 
  unlist() %>% 
  matrix(ncol = 9, byrow = TRUE) %>% 
  `colnames<-`(clnms)

m1_out <- apply(m1_boot_mat, 2, quantile, probs = c(0.025, 0.975)) %>% 
  t() %>% 
  as_tibble() %>%
  bind_cols(est = as.numeric(est_m1), .) %>% 
  mutate_all( ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", `2.5%`, ", ", `97.5%`, ")")) %>% 
  mutate(out = paste(est, intv), name = clnms) %>% 
  dplyr::select(name, est, `2.5%`, `97.5%`, out)

#compute bootstrap estimates for M2

#setup parallel computing cluster
ncores <- parallel::detectCores()-1
my.cluster <- parallel::makeCluster(ncores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("linmed", "linmedx"), envir=environment())
registerDoRNG(3308004)

m2_boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {
  
  boot_data <- nlsy4m2[sample(nrow(nlsy4m2), nrow(nlsy4m2), replace=TRUE),]
  
  boot_linmed <- linmed(boot_data, d, m2, y, x)
  
  boot_linmedx <- linmedx(boot_data, d, m2, y, x)
  
  boot_ipwmed <- ipwmed(boot_data, d, m2, y, x)
  
  list(linmed = boot_linmed,
       linmedx = boot_linmedx,
       ipwmed = boot_ipwmed)
}

stopCluster(my.cluster)
rm(my.cluster)

clnms <- expand_grid(method = c("linmed", "linmedx", "ipwmed"), estimand = c("ATE", "NDE", "NIE")) %>% unite("") %>% unlist()

m2_boot_mat <- m2_boot %>% 
  map(unlist) %>% 
  unlist() %>% 
  matrix(ncol = 9, byrow = TRUE) %>% 
  `colnames<-`(clnms)

m2_out <- apply(m2_boot_mat, 2, quantile, probs = c(0.025, 0.975)) %>% 
  t() %>% 
  as_tibble() %>%
  bind_cols(est = as.numeric(est_m2), .) %>% 
  mutate_all( ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", `2.5%`, ", ", `97.5%`, ")")) %>% 
  mutate(out = paste(est, intv), name = clnms) %>% 
  dplyr::select(name, est, `2.5%`, `97.5%`, out)

write_csv(m1_out, file = "table5-1a.csv")
write_csv(m2_out, file = "table5-1b.csv")

# diagnostic test for M2-M1 given X and D

diag_form <- paste(m2, "~", paste(c(x, d, m1), collapse = " + "))

diag_lm <- lm(diag_form, data=nlsy4m2)

summary(diag_lm)

#sink()
