###Table 5.4###

rm(list=ls())

packages<-c("tidyverse", "margins", "mediation", "foreach", "doParallel", "doRNG")

for (package.i in packages) {
  suppressPackageStartupMessages(library(package.i, character.only=TRUE))
}
source("utils.R")

##office
datadir <- "../../data/" 
logdir <- "../../code/ch5/_LOGS/"

##input data
nlsy_raw <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy_raw[complete.cases(nlsy_raw[,c("ever_unemp_age3539", "momedu", "att22", "cesd_age40", "faminc_adj_age3539")]),]

nlsy$momcol <- ifelse(nlsy$momedu>12, 1, 0)

nlsy$incgt50k <- ifelse(nlsy$faminc_adj_age3539>=50000, 1, 0)

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##tabulate data
nptab <- nlsy %>%
  group_by(momcol, att22, incgt50k, ever_unemp_age3539) %>% 
  dplyr::summarize(
    mean = mean(std_cesd_age40),
    n = n(),
    .groups = "drop")

print(nptab)

# GMF function

gmf <- function(d1, d2, d){
  
  mu_y <- nlsy %>% 
    filter(att22 == d) %>% 
    group_by(momcol, ever_unemp_age3539, incgt50k) %>% 
    summarise(
      mu_y = mean(std_cesd_age40),
      n = n(),
      .groups = "drop"
    ) %>% 
    dplyr::select(-n)
  
  pi_m2 <- nlsy %>% 
    filter(att22 == d2) %>% 
    group_by(momcol, ever_unemp_age3539, incgt50k) %>% 
    summarise(
      n = n(),
      .groups = "drop"
    ) %>% 
    group_by(momcol, ever_unemp_age3539) %>% 
    mutate(
      pi_m2 = n/sum(n)
    ) %>% 
    dplyr::select(-n)
  
  pi_m1 <- nlsy %>% 
    filter(att22 == d1) %>% 
    group_by(momcol, ever_unemp_age3539, incgt50k) %>% 
    summarise(
      n = n(),
      .groups = "drop"
    ) %>% 
    ungroup() %>% 
    group_by(momcol, ever_unemp_age3539) %>% 
    mutate(
      n_m1 = sum(n)
    ) %>% 
    ungroup() %>% 
    group_by(momcol) %>% 
    mutate(
      pi_m1 = n_m1/sum(n)
    ) %>% 
    dplyr::select(-n, -n_m1)
  
  pi_c <- nlsy %>% 
    group_by(momcol, ever_unemp_age3539, incgt50k) %>% 
    summarise(
      n = n(),
      .groups = "drop"
    ) %>% 
    ungroup() %>% 
    group_by(momcol) %>% 
    mutate(
      n_c = sum(n)
    ) %>% 
    ungroup() %>% 
    mutate(
      pi_c = n_c/sum(n)
    ) %>% 
    dplyr::select(-n, -n_c)
  
  df <- mu_y %>% 
    left_join(pi_m2, by = c("momcol", "ever_unemp_age3539", "incgt50k")) %>% 
    left_join(pi_m1, by = c("momcol", "ever_unemp_age3539", "incgt50k")) %>% 
    left_join(pi_c, by = c("momcol", "ever_unemp_age3539", "incgt50k"))
  
  with(df, sum(mu_y * pi_m2 * pi_m1 * pi_c))
  
}

ATEhat_npl <- gmf(1,1,1) - gmf(0,0,0)
AYhat_npl <- gmf(0,0,1) - gmf(0,0,0)
AM2Yhat_npl <- gmf(0,1,1) - gmf(0,0,1)
AM1Yhat_npl <- gmf(1,1,1) - gmf(0,1,1)

print(c(ATEhat_npl, AYhat_npl, AM2Yhat_npl, AM1Yhat_npl))