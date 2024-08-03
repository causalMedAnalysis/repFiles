###Table 5.2###

rm(list=ls())

packages<-c("tidyverse", "margins", "mediation", "foreach", "doParallel", "doRNG")

for (package.i in packages) {
  suppressPackageStartupMessages(library(package.i, character.only=TRUE))
}
source("utils.R")

##office
datadir <- "../../data/" 
logdir <- "../../code/ch5/_LOGS/"

#sink(paste(logdir, "table_5-2_log.txt", sep=""))

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

## linear model estimators

df <- nlsy
m <- c(m1, m2)

#linear model w/o exposure-mediator interaction

linmed_m <- linmed(df, d, m, y, x)

#linear model w exposure-mediator interaction

linmedx_m <- linmedx(df, d, m, y, x)

# ipw estimator

ipwmed_m <- ipwmed(df, d, m, y, x)

# summary of point estimates

est_m <- list(linmed_m, linmedx_m, ipwmed_m) %>%
  map(unlist) %>%
  bind_rows() %>%
  t() %>% 
  `colnames<-`( c("^lma", "^lmi", "^ipw"))

# bootstrap SEs for linmed and linmedx and ipwmed

#setup parallel computing cluster
ncores <- parallel::detectCores()-1
my.cluster <- parallel::makeCluster(ncores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("linmed", "linmedx", "ipwmed"), envir=environment())
registerDoRNG(3308004)

nboot <- 2000

m_boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {
  
  boot_data <- nlsy[sample(nrow(nlsy), nrow(nlsy), replace=TRUE),]
  
  boot_linmed <- linmed(boot_data, d, m, y, x)
  
  boot_linmedx <- linmedx(boot_data, d, m, y, x)
  
  boot_ipwmed <- ipwmed(boot_data, d, m, y, x)
  
  list(linmed = boot_linmed,
       linmedx = boot_linmedx,
       ipwmed = boot_ipwmed)
}

stopCluster(my.cluster)
rm(my.cluster)

clnms <- expand_grid(method = c("linmed", "linmedx", "ipwmed"), estimand = c("ATE", "NDE", "NIE")) %>% unite("") %>% unlist()

m_boot_mat <- m_boot %>% 
  map(unlist) %>% 
  unlist() %>% 
  matrix(ncol = 9, byrow = TRUE) %>% 
  `colnames<-`(clnms)

m_out <- apply(m_boot_mat, 2, quantile, probs = c(0.025, 0.975)) %>% 
  t() %>% 
  as_tibble() %>%
  bind_cols(est = as.numeric(est_m), .) %>% 
  mutate_all( ~ round(.x, digits = 3)) %>% 
  mutate(intv = paste0("(", `2.5%`, ", ", `97.5%`, ")")) %>% 
  mutate(out = paste(est, intv), name = clnms) %>% 
  dplyr::select(name, est, `2.5%`, `97.5%`, out)

write_csv(m_out, file = "table5-2.csv")

#sink()
