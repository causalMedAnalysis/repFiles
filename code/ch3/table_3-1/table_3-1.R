###Table 3.1###

rm(list=ls())

packages<-c("dplyr", "tidyr", "margins")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

##office
datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

##home 
#datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

sink(paste(logdir, "table_3-1_log.txt", sep=""))

##input data
nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("ever_unemp_age3539", "momedu", "att22", "cesd_age40")]),]

nlsy$momcol <- ifelse(nlsy$momedu>12, 1, 0)

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##tabulate data
nlsy %>%
	group_by(momcol,att22,ever_unemp_age3539) %>% 
			dplyr::summarize(
				mean = mean(std_cesd_age40),
				n = n(),
				.groups = "drop")

##nonparametric estimators
cat("\n")

#ATEhat^np
m1 <- lm(std_cesd_age40~att22*momcol, data=nlsy)
ATEhat <- mean(marginal_effects(m1, nlsy, variables=c("att22"))$dydx_att22)

#CDEhat^np
m2 <- lm(std_cesd_age40~att22*ever_unemp_age3539*momcol, data=nlsy)

gdata <- nlsy

gdata$att22 <- 1
gdata$ever_unemp_age3539 <- 0
Ehat_Y10 <- predict(m2, gdata)

gdata$att22 <- 0
Ehat_Y00 <- predict(m2, gdata)

gdata$ever_unemp_age3539 <- 1
Ehat_Y01 <- predict(m2, gdata)

gdata$att22 <- 1
Ehat_Y11 <- predict(m2, gdata)

CDE0hat <- mean(Ehat_Y10-Ehat_Y00)
CDE1hat <- mean(Ehat_Y11-Ehat_Y01)

#NDEhat^np and NIEhat^np
gdata$att22 <- 0
gdata$ever_unemp_age3539 <- 0
EhatY_D0M0C <- predict(m2, gdata)

gdata$ever_unemp_age3539 <- 1
EhatY_D0M1C <- predict(m2, gdata)

gdata$att22 <- 1
gdata$ever_unemp_age3539 <- 0
EhatY_D1M0C <- predict(m2, gdata)

gdata$ever_unemp_age3539 <- 1
EhatY_D1M1C <- predict(m2, gdata)

m3 <- lm(ever_unemp_age3539~att22*momcol, data=nlsy)
gdata <- nlsy

gdata$att22 <- 0
PhatM_D0C <- predict(m3, gdata)

gdata$att22 <- 1
PhatM_D1C <- predict(m3, gdata)

NDEhat <- mean((EhatY_D1M0C-EhatY_D0M0C)*(1-PhatM_D0C)+(EhatY_D1M1C-EhatY_D0M1C)*PhatM_D0C)
NIEhat <- mean((EhatY_D1M0C)*((1-PhatM_D1C)-(1-PhatM_D0C))+(EhatY_D1M1C)*((PhatM_D1C)-(PhatM_D0C)))

npest.output <- as.matrix(c(ATEhat, NDEhat, NIEhat, CDE0hat, CDE1hat))

for (i in 1:length(npest.output)) { 
	npest.output[i,1]<-round(npest.output[i,1], digits=3)
	}

npest.output <- data.frame(npest.output, row.names=c("ATEhat^np", "NDEhat^np", "NIEhat^np", "CDE0hat^np", "CDE1hat^np"))
colnames(npest.output) <- c("point.est")

print(npest.output)

sink()
