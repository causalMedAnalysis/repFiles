###Table 4.1###

rm(list=ls())

packages<-c("dplyr", "tidyr")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

##office
datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch4/_LOGS/"

##home 
#datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch4/_LOGS/"

sink(paste(logdir, "table_4-1_log.txt", sep=""))

##input data
nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("ever_unemp_age3539", "momedu", "att22", "cesd_age40", "faminc_adj_age3539")]),]

nlsy$momcol <- ifelse(nlsy$momedu>12, 1, 0)

nlsy$incgt50k <- ifelse(nlsy$faminc_adj_age3539>=50000, 1, 0)

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##tabulate data
nptab <-nlsy %>%
		group_by(momcol, att22, incgt50k, ever_unemp_age3539) %>% 
			dplyr::summarize(
				mean = mean(std_cesd_age40),
				n = n(),
				.groups = "drop")

print(nptab)

##compute nonparametric estimates
YhatC0D0L0M0 <- nptab[1,5] 
YhatC0D0L1M0 <- nptab[2,5] 
YhatC0D0L0M1 <- nptab[3,5] 
YhatC0D0L1M1 <- nptab[4,5] 

YhatC0D1L0M0 <- nptab[5,5] 
YhatC0D1L1M0 <- nptab[6,5] 
YhatC0D1L0M1 <- nptab[7,5] 
YhatC0D1L1M1 <- nptab[8,5] 

YhatC1D0L0M0 <- nptab[9,5] 
YhatC1D0L1M0 <- nptab[10,5] 
YhatC1D0L0M1 <- nptab[11,5] 
YhatC1D0L1M1 <- nptab[12,5] 

YhatC1D1L0M0 <- nptab[13,5] 
YhatC1D1L1M0 <- nptab[14,5] 
YhatC1D1L0M1 <- nptab[15,5] 
YhatC1D1L1M1 <- nptab[16,5] 

phatL1_C0D0 <- mean(nlsy[which(nlsy$momcol==0 & nlsy$att22==0), "ever_unemp_age3539"])
phatL1_C0D1 <- mean(nlsy[which(nlsy$momcol==0 & nlsy$att22==1), "ever_unemp_age3539"])
phatL1_C1D0 <- mean(nlsy[which(nlsy$momcol==1 & nlsy$att22==0), "ever_unemp_age3539"])
phatL1_C1D1 <- mean(nlsy[which(nlsy$momcol==1 & nlsy$att22==1), "ever_unemp_age3539"])

phatM1_C0D0 <- mean(nlsy[which(nlsy$momcol==0 & nlsy$att22==0), "incgt50k"])
phatM1_C0D1 <- mean(nlsy[which(nlsy$momcol==0 & nlsy$att22==1), "incgt50k"])
phatM1_C1D0 <- mean(nlsy[which(nlsy$momcol==1 & nlsy$att22==0), "incgt50k"])
phatM1_C1D1 <- mean(nlsy[which(nlsy$momcol==1 & nlsy$att22==1), "incgt50k"])

phatC1 <- mean(nlsy$momcol)

CDEhat_npl <- 
	((YhatC0D1L0M1*(1-phatL1_C0D1) + YhatC0D1L1M1*phatL1_C0D1)*(1-phatC1) + (YhatC1D1L0M1*(1-phatL1_C1D1) + YhatC1D1L1M1*phatL1_C1D1)*(phatC1)) -
	((YhatC0D0L0M1*(1-phatL1_C0D0) + YhatC0D0L1M1*phatL1_C0D0)*(1-phatC1) + (YhatC1D0L0M1*(1-phatL1_C1D0) + YhatC1D0L1M1*phatL1_C1D0)*(phatC1))

IDEhat_npl <- 
	((YhatC0D1L0M1*(1-phatL1_C0D1) + YhatC0D1L1M1*phatL1_C0D1)-(YhatC0D0L0M1*(1-phatL1_C0D0) + YhatC0D0L1M1*phatL1_C0D0))*phatM1_C0D0*(1-phatC1) + 
	((YhatC1D1L0M1*(1-phatL1_C1D1) + YhatC1D1L1M1*phatL1_C1D1)-(YhatC1D0L0M1*(1-phatL1_C1D0) + YhatC1D0L1M1*phatL1_C1D0))*phatM1_C1D0*(phatC1) +
	((YhatC0D1L0M0*(1-phatL1_C0D1) + YhatC0D1L1M0*phatL1_C0D1)-(YhatC0D0L0M0*(1-phatL1_C0D0) + YhatC0D0L1M0*phatL1_C0D0))*(1-phatM1_C0D0)*(1-phatC1) +
	((YhatC1D1L0M0*(1-phatL1_C1D1) + YhatC1D1L1M0*phatL1_C1D1)-(YhatC1D0L0M0*(1-phatL1_C1D0) + YhatC1D0L1M0*phatL1_C1D0))*(1-phatM1_C1D0)*(phatC1)

IIEhat_npl <- 
	(YhatC0D1L0M1*(1-phatL1_C0D1) + YhatC0D1L1M1*phatL1_C0D1)*(1-phatC1)*(phatM1_C0D1-phatM1_C0D0) +
	(YhatC1D1L0M1*(1-phatL1_C1D1) + YhatC1D1L1M1*phatL1_C1D1)*phatC1*(phatM1_C1D1-phatM1_C1D0) +
	(YhatC0D1L0M0*(1-phatL1_C0D1) + YhatC0D1L1M0*phatL1_C0D1)*(1-phatC1)*((1-phatM1_C0D1)-(1-phatM1_C0D0)) +
	(YhatC1D1L0M0*(1-phatL1_C1D1) + YhatC1D1L1M0*phatL1_C1D1)*phatC1*((1-phatM1_C1D1)-(1-phatM1_C1D0))

OEhat_npl <- IDEhat_npl + IIEhat_npl

print(c(CDEhat_npl, IDEhat_npl, IIEhat_npl, OEhat_npl))

sink()
