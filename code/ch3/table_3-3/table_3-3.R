###Table 3.3###

rm(list=ls())

packages<-c("dplyr", "tidyr", "mediation")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

nsim<-2000

set.seed(3308004)

##office
datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

##home
#datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

sink(paste(logdir, "table_3-3_log.txt", sep=""))

##input data
nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("cesd_age40", "ever_unemp_age3539", "att22", "female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")]),]

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##simulation/imputation estimators w/o covariate interactions

Mmodel <- glm(ever_unemp_age3539~att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=nlsy, family=binomial("logit"))

Ymodel <- lm(std_cesd_age40~ever_unemp_age3539+att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=nlsy)

simest <- mediate(Mmodel, Ymodel, treat="att22", mediator="ever_unemp_age3539", sims=nsim)

impcde <- function(data) {
	df <- data

	Ymodel <- lm(std_cesd_age40~ever_unemp_age3539+att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df)

	gdata <- df

	gdata$att22 <- 1
	gdata$ever_unemp_age3539 <- 0
	Yhat10 <- predict(Ymodel, gdata)

	gdata$att22 <- 0
	Yhat00 <- predict(Ymodel, gdata)

	point.est <- mean(Yhat10-Yhat00)

	return(point.est)
	}

impcde.est <- impcde(nlsy)

##simulation estimator w/ covariate interactions

Mmodelx <- glm(ever_unemp_age3539~att22*(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3), data=nlsy, family=binomial("logit"))

Ymodelx <- lm(std_cesd_age40~(ever_unemp_age3539*att22)+(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3)*(ever_unemp_age3539+att22), data=nlsy)

simestx <- mediate(Mmodelx, Ymodelx, treat="att22", mediator="ever_unemp_age3539", sims=nsim)

impcdex <- function(data) {
	df <- data

	Ymodelx <- lm(std_cesd_age40~(ever_unemp_age3539*att22)+(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3)*(ever_unemp_age3539+att22), data=df)

	gdata <- df

	gdata$att22 <- 1
	gdata$ever_unemp_age3539 <- 0
	Yhat10 <- predict(Ymodelx, gdata)

	gdata$att22 <- 0
	Yhat00 <- predict(Ymodelx, gdata)

	point.est <- mean(Yhat10-Yhat00)

	return(point.est)
	}

impcdex.est <- impcdex(nlsy)

##print output
sim.imp.est.output <- matrix(data=NA, nrow=4, ncol=2)

sim.imp.est.output[1,]<-c(round(simest$tau.coef, digits=3), round(simestx$tau.coef, digits=3))
sim.imp.est.output[2,]<-c(round(simest$z0, digits=3), round(simestx$z0, digits=3))
sim.imp.est.output[3,]<-c(round(simest$d1, digits=3), round(simestx$d1, digits=3))
sim.imp.est.output[4,]<-c(round(impcde.est, digits=3), round(impcdex.est, digits=3))

sim.imp.est.output <- data.frame(sim.imp.est.output, row.names=c("ATEhat^sim", "NDEhat^sim", "NIEhat^sim", "CDE0hat^imp"))
colnames(sim.imp.est.output) <- c("Additive Models", "Models w/ interations")

print(sim.imp.est.output)

sink()
