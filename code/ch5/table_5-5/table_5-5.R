###Table 5.5###

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

## linear model estimator w/o exposure-mediator interaction

linmed_m1 <-  linmed(df, d, m1, y, x)
linmed_m12 <-  linmed(df, d, m12, y, x)

decomp_linmed <- c(ATE = linmed_m1$ATE,
                   AY = linmed_m12$NDE,
                   AM2Y = linmed_m1$NDE - linmed_m12$NDE,
                   AM1Y = linmed_m1$NIE)

## linear model estimator w/ exposure-mediator interactions

linmedx_m1 <-  linmedx(df, d, m1, y, x)
linmedx_m12 <-  linmedx(df, d, m12, y, x)

decomp_linmedx <- c(ATE = linmedx_m1$ATE,
                    AY = linmedx_m12$NDE,
                    AM2Y = linmedx_m1$NDE - linmedx_m12$NDE,
                    AM1Y = linmedx_m1$NIE)

## IPW estimator

ipwmed_m1 <- ipwmed(df, d, m1, y, x)
ipwmed_m12 <- ipwmed(df, d, m12, y, x)

decomp_ipwmed <- c(ATE = ipwmed_m1$ATE,
                   AY = ipwmed_m12$NDE,
                   AM2Y = ipwmed_m1$NDE - ipwmed_m12$NDE,
                   AM1Y = ipwmed_m1$NIE)

# summary of point estimates

est_m <- list(decomp_linmed, decomp_linmedx, decomp_ipwmed) %>%
  map(unlist) %>%
  bind_rows() %>%
  t() %>% 
  `colnames<-`( c("^lma", "^lmi", "^ipw"))

## Bootstrap SEs for linmed, linmed, and ipwmed

# setup parallel computing cluster
ncores <- parallel::detectCores()-1
my.cluster <- parallel::makeCluster(ncores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("linmed", "linmedx", "ipwmed"), envir=environment())
registerDoRNG(3308004)

nboot <- 2000

m_boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {
  
  boot_data <- nlsy[sample(nrow(nlsy), nrow(nlsy), replace=TRUE),]
  
  boot_linmed_m1 <- linmed(boot_data, d, m1, y, x)
  
  boot_linmedx_m1 <- linmedx(boot_data, d, m1, y, x)
  
  boot_ipwmed_m1 <- ipwmed(boot_data, d, m1, y, x)
  
  boot_linmed_m12 <- linmed(boot_data, d, m12, y, x)
  
  boot_linmedx_m12 <- linmedx(boot_data, d, m12, y, x)
  
  boot_ipwmed_m12 <- ipwmed(boot_data, d, m12, y, x)
  
  decomp_linmed <- c(ATE = boot_linmed_m1$ATE,
                     AY = boot_linmed_m12$NDE,
                     AM2Y = boot_linmed_m1$NDE - boot_linmed_m12$NDE,
                     AM1Y = boot_linmed_m1$NIE)
  
  decomp_linmedx <- c(ATE = boot_linmedx_m1$ATE,
                      AY = boot_linmedx_m12$NDE,
                      AM2Y = boot_linmedx_m1$NDE - boot_linmedx_m12$NDE,
                      AM1Y = boot_linmedx_m1$NIE)
  
  decomp_ipwmed <- c(ATE = boot_ipwmed_m1$ATE,
                     AY = boot_ipwmed_m12$NDE,
                     AM2Y = boot_ipwmed_m1$NDE - boot_ipwmed_m12$NDE,
                     AM1Y = boot_ipwmed_m1$NIE)
  
  list(decomp_linmed,
       decomp_linmedx,
       decomp_ipwmed)
}

stopCluster(my.cluster)
rm(my.cluster)

clnms <- expand_grid(method = c("linmed", "linmedx", "ipwmed"), estimand = c("ATE", "AY", "AM2Y", "AM1Y")) %>% unite("") %>% unlist()

m_boot_mat <- m_boot %>% 
  map(unlist) %>% 
  unlist() %>% 
  matrix(ncol = 12, byrow = TRUE) %>% 
  `colnames<-`(clnms)

m_out <- apply(m_boot_mat, 2, quantile, probs = c(0.025, 0.975)) %>% 
  t() %>% 
  as_tibble() %>%
  bind_cols(est = as.numeric(est_m), .) %>% 
  mutate_all( ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", `2.5%`, ", ", `97.5%`, ")")) %>% 
  mutate(out = paste(est, intv), name = clnms) %>% 
  dplyr::select(name, est, `2.5%`, `97.5%`, out)

write_csv(m_out, file = "table5-5.csv")
