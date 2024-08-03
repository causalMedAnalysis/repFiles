###Table 3.7###

rm(list=ls())

packages<-c("dplyr", "tidyr", "foreach", "doParallel", "doRNG", "boot", "mediation", "Hmisc", "foreign")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

nboot<-2000
nsim<-1000
ncores<-parallel::detectCores()-1
set.seed(3308004)

##office
#datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

##home
datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

sink(paste(logdir, "table_3-7_log.txt", sep=""))

##input data
jobs <- read.dta(paste(datadir, "JOBSII/Jobs-NoMiss-Binary.dta", sep=""))

jobs <- jobs %>% 
	mutate(
		treat = recode(treat, "control" = 0, "exp" = 1),
		work1 = recode(work1, "psyump" = 0, "psyemp" = 1),
		nonwhite = recode(nonwhite, "white0" = 0, "non.white1" = 1),
		educ = recode(educ, "lt-hs" = 1, "highsc" = 2, "somcol" = 3, "bach" = 4, "gradwk" = 5),
		income = recode(income, "lt15k" = 1, "15t24k" = 2, "25t39k" = 3, "40t49k" = 4, "50k+" = 5))

##linear models w/ exposure-mediator interaction

linmedx <- function(data) {
	df <- data

	df <- df %>% 
		mutate(
			econ_hard = econ_hard-mean(econ_hard), 
			sex = sex-mean(sex), 
			age = age-mean(age), 
			nonwhite = nonwhite-mean(nonwhite), 
			educ = educ-mean(educ), 
			income = income-mean(income))

	Mmodel <- lm(job_seek~treat+econ_hard+sex+age+nonwhite+educ+income, data=df)

	Ymodel <- lm(work1~(job_seek*treat)+econ_hard+sex+age+nonwhite+educ+income, data=df)

	NDE <- (Ymodel$coefficients["treat"] + Ymodel$coefficients["job_seek:treat"]*Mmodel$coefficients["(Intercept)"])
	NIE <- Mmodel$coefficients["treat"]*(Ymodel$coefficients["job_seek"] + Ymodel$coefficients["job_seek:treat"])
	ATE <- NDE+NIE
	CDE4 <- Ymodel$coefficients["treat"] + Ymodel$coefficients["job_seek:treat"]*4
	
	point.est <- list(ATE, NDE, NIE, CDE4)

	return(point.est)
	}

linmedx.est <- linmedx(jobs)
linmedx.est <- matrix(unlist(linmedx.est), ncol=4, byrow=TRUE)

my.cluster<-parallel::makeCluster(ncores, type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
invisible(clusterEvalQ(cl=my.cluster, library(dplyr)))
clusterExport(cl=my.cluster, list("linmedx"), envir=environment())
registerDoRNG(3308004)

linmedx.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {
	boot.data <- jobs[sample(nrow(jobs), nrow(jobs), replace=TRUE),]

	boot.est <- linmedx(boot.data)

	return(boot.est)
	}

stopCluster(my.cluster)
rm(my.cluster)

linmedx.boot <- matrix(unlist(linmedx.boot), ncol=4, byrow=TRUE)

linmedx.output <- matrix(data=NA, nrow=4, ncol=3)

for (i in 1:4) { 
	linmedx.output[i,1]<-round(linmedx.est[i], digits=3)
	linmedx.output[i,2]<-round(quantile(linmedx.boot[,i], prob=0.025), digits=3)
	linmedx.output[i,3]<-round(quantile(linmedx.boot[,i], prob=0.975), digits=3)
	}

linmedx.output <- data.frame(linmedx.output, row.names=c("ATEhat^lmi", "NDEhat^lmi", "NIEhat^lmi", "CDE4hat^lmi"))
colnames(linmedx.output) <- c("point.est", "ll.95ci", "ul.95ci")

##simulation estimators

Mmodel <- lm(job_seek~treat+econ_hard+sex+age+nonwhite+educ+income, data=jobs)

Ymodel <- glm(work1~(job_seek*treat)+econ_hard+sex+age+nonwhite+educ+income, data=jobs, family=binomial("logit"))

simest <- mediate(Mmodel, Ymodel, boot=TRUE, treat="treat", mediator="job_seek", sims=nsim, parallel="multicore")

simcde <- function(data, indices) {
	df <- data[indices,]

	Ymodel <- glm(work1~(job_seek*treat)+econ_hard+sex+age+nonwhite+educ+income, data=df, family=binomial("logit"))

	gdata <- df

	gdata$treat <- 1
	gdata$job_seek <- 4
	Ehat_Y14 <- predict(Ymodel, gdata, type = "response")

	gdata$treat <- 0
	Ehat_Y04 <- predict(Ymodel, gdata, type = "response")

	point.est <- mean(Ehat_Y14-Ehat_Y04)

	return(point.est)
	}

simcde.boot <- boot(data=jobs, statistic=simcde, R=nboot)

simest.output <- matrix(data=NA, nrow=4, ncol=3)
simest.output[1,]<-c(round(simest$tau.coef, digits=3), round(simest$tau.ci[1], digits=3), round(simest$tau.ci[2], digits=3))
simest.output[2,]<-c(round(simest$z0, digits=3), round(simest$z0.ci[1], digits=3), round(simest$z0.ci[2], digits=3))
simest.output[3,]<-c(round(simest$d1, digits=3), round(simest$d1.ci[1], digits=3),round(simest$d1.ci[2], digits=3))
simest.output[4,]<-c(round(simcde.boot$t0, digits=3), round(quantile(simcde.boot$t, prob=0.025), digits=3), round(quantile(simcde.boot$t, prob=0.975), digits=3))

simest.output <- data.frame(simest.output, row.names=c("ATEhat^sim", "NDEhat^sim", "NIEhat^sim", "CDE4hat^sim"))
colnames(simest.output) <- c("point.est", "ll.95ci", "ul.95ci")

##inverse probability weighting estimators

ipwmed <- function(data) {
	df <- data

	Dmodel_1 <- glm(treat~econ_hard+sex+age+nonwhite+educ+income, data=df, family=binomial("logit"))

	Dmodel_2 <- glm(treat~job_seek+econ_hard+sex+age+nonwhite+educ+income, data=df, family=binomial("logit"))

	df$phat_d_C <- predict(Dmodel_1, type = "response")
	df$phat_dstar_C <- 1-df$phat_d_C

	df$phat_d_CM <- predict(Dmodel_2, type = "response")
	df$phat_dstar_CM <- 1-df$phat_d_CM

	df$phat_d <- mean(df$treat)
	df$phat_dstar <- 1-df$phat_d

	df$sw1[df$treat==0] <- df$phat_dstar[df$treat==0]/df$phat_dstar_C[df$treat==0]
	df$sw2[df$treat==1] <- df$phat_d[df$treat==1]/df$phat_d_C[df$treat==1]
	df$sw3[df$treat==1] <- (df$phat_dstar_CM[df$treat==1]*df$phat_d[df$treat==1])/(df$phat_d_CM[df$treat==1]*df$phat_dstar_C[df$treat==1])

	for (i in c("sw1", "sw2", "sw3")) {
		df[which(df[,i] > quantile(df[,i], probs=0.99, na.rm=T)), i] <- quantile(df[,i], probs=0.99, na.rm=T)
		df[which(df[,i] < quantile(df[,i], probs=0.01, na.rm=T)), i] <- quantile(df[,i], probs=0.01, na.rm=T)
		}

	Ehat_Y0M0 <- weighted.mean(df$work1[df$treat==0], df$sw1[df$treat==0])
	Ehat_Y1M1 <- weighted.mean(df$work1[df$treat==1], df$sw2[df$treat==1])
	Ehat_Y1M0 <- weighted.mean(df$work1[df$treat==1], df$sw3[df$treat==1])

	ATE <- Ehat_Y1M1-Ehat_Y0M0
	NDE <- Ehat_Y1M0-Ehat_Y0M0
	NIE <- Ehat_Y1M1-Ehat_Y1M0

	point.est <- list(ATE, NDE, NIE)

	return(point.est)
	}

ipwmed.est <- ipwmed(jobs)
ipwmed.est <- matrix(unlist(ipwmed.est), ncol=3, byrow=TRUE)

my.cluster<-parallel::makeCluster(ncores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("ipwmed"), envir=environment())
registerDoRNG(3308004)

ipwmed.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {
	boot.data <- jobs[sample(nrow(jobs), nrow(jobs), replace=TRUE),]

	boot.est <- ipwmed(boot.data)

	return(boot.est)
	}

stopCluster(my.cluster)
rm(my.cluster)

ipwmed.boot <- matrix(unlist(ipwmed.boot), ncol=3, byrow=TRUE)

ipwcde <- function(data, indices) {

	df <- data[indices,]

	Dmodel <- glm(treat~econ_hard+sex+age+nonwhite+educ+income, data=df, family=binomial("logit"))

	Mmodel_1 <- lm(job_seek~treat+econ_hard+sex+age+nonwhite+educ+income, data=df)

	df$phat_d_C <- predict(Dmodel, type = "response")
	df$phat_D_denom <- (df$treat*df$phat_d_C)+((1-df$treat)*(1-df$phat_d_C))

	df$Ehat_m_CD <- predict(Mmodel_1)
	df$sighat_m_CD <- sqrt(mean(Mmodel_1$residuals^2))
	df$phat_M_denom <- dnorm(df$job_seek, mean=df$Ehat_m_CD, sd=df$sighat_m_CD)

	df$phat_d <- mean(df$treat)
	df$phat_D_num <- (df$treat*df$phat_d)+((1-df$treat)*(1-df$phat_d))

	Mmodel_2 <- lm(job_seek~treat, data=df)

	df$Ehat_m_D <- predict(Mmodel_2)
	df$sighat_m_D <- sqrt(mean(Mmodel_2$residuals^2))
	df$phat_M_num <- dnorm(df$job_seek, mean=df$Ehat_m_D, sd=df$sighat_m_D)

	df$sw4 <- (df$phat_M_num*df$phat_D_num)/(df$phat_M_denom*df$phat_D_denom)

	df[which(df[,"sw4"] > quantile(df[,"sw4"], probs=0.99, na.rm=T)), "sw4"] <- quantile(df[, "sw4"], probs=0.99, na.rm=T)
	df[which(df[,"sw4"] < quantile(df[,"sw4"], probs=0.01, na.rm=T)), "sw4"] <- quantile(df[, "sw4"], probs=0.01, na.rm=T)

	Ymodel <- lm(work1~treat*job_seek, data=df, weights=sw4)

	point.est <- Ymodel$coefficients["treat"] + Ymodel$coefficients["treat:job_seek"]*4

	return(point.est)
	}

ipwcde.boot <- boot(data=jobs, statistic=ipwcde, R=nboot)

ipwest.output <- matrix(data=NA, nrow=4, ncol=3)

for (i in 1:3) {
	ipwest.output[i,]<-c(round(ipwmed.est[,i], digits=3), round(quantile(ipwmed.boot[,i], prob=0.025), digits=3), round(quantile(ipwmed.boot[,i], prob=0.975), digits=3))
	}

ipwest.output[4,]<-c(round(ipwcde.boot$t0, digits=3), round(quantile(ipwcde.boot$t, prob=0.025), digits=3), round(quantile(ipwcde.boot$t, prob=0.975), digits=3))

ipwest.output <- data.frame(ipwest.output, row.names=c("ATEhat^ipw", "NDEhat^ipw", "NIEhat^ipw", "CDE4hat^ipw"))
colnames(ipwest.output) <- c("point.est", "ll.95ci", "ul.95ci")

##print results
print(linmedx.output)
cat("\n")
print(simest.output)
cat("\n")
print(ipwest.output)

sink()
