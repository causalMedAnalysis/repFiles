###Table 3.5###

rm(list=ls())

packages<-c("dplyr", "tidyr", "foreach", "doParallel", "doRNG")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

nboot <- 2000

ncores <- parallel::detectCores()-1

##office
datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

##home
#datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

sink(paste(logdir, "table_3-5_log.txt", sep=""))

##input data 
nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("cesd_age40", "ever_unemp_age3539", "att22", "female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")]),]

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##inverse probability weighting estimators

#define functions for estimators

ipwmed <- function(data) {
	
	df <- data

	Dmodel_1 <- glm(att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df, family=binomial("logit"))

	Dmodel_2 <- glm(att22~ever_unemp_age3539+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df, family=binomial("logit"))

	df$phat_d_C <- predict(Dmodel_1, type = "response")
	df$phat_dstar_C <- 1-df$phat_d_C

	df$phat_d_CM <- predict(Dmodel_2, type = "response")
	df$phat_dstar_CM <- 1-df$phat_d_CM

	df$phat_d <- mean(df$att22)
	df$phat_dstar <- 1-df$phat_d

	df$sw1[df$att22==0] <- df$phat_dstar[df$att22==0]/df$phat_dstar_C[df$att22==0]
	df$sw2[df$att22==1] <- df$phat_d[df$att22==1]/df$phat_d_C[df$att22==1]
	df$sw3[df$att22==1] <- (df$phat_dstar_CM[df$att22==1]*df$phat_d[df$att22==1])/(df$phat_d_CM[df$att22==1]*df$phat_dstar_C[df$att22==1])

	for (i in c("sw1", "sw2", "sw3")) {
		df[which(df[,i] > quantile(df[,i], probs=0.99, na.rm=T)), i] <- quantile(df[,i], probs=0.99, na.rm=T)
		df[which(df[,i] < quantile(df[,i], probs=0.01, na.rm=T)), i] <- quantile(df[,i], probs=0.01, na.rm=T)
		}

	Ehat_Y0M0 <- weighted.mean(df$std_cesd_age40[df$att22==0], df$sw1[df$att22==0])
	Ehat_Y1M1 <- weighted.mean(df$std_cesd_age40[df$att22==1], df$sw2[df$att22==1])
	Ehat_Y1M0 <- weighted.mean(df$std_cesd_age40[df$att22==1], df$sw3[df$att22==1])

	ATE <- Ehat_Y1M1-Ehat_Y0M0
	NDE <- Ehat_Y1M0-Ehat_Y0M0
	NIE <- Ehat_Y1M1-Ehat_Y1M0

	point.est <- list(ATE, NDE, NIE)

	return(point.est)
	}

ipwcde <- function(data) {

	df <- data

	Dmodel <- glm(att22~female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df, family=binomial("logit"))

	Mmodel_1 <- glm(ever_unemp_age3539~att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df, family=binomial("logit"))

	df$phat_d_C <- predict(Dmodel, type = "response")
	df$phat_D_denom <- (df$att22*df$phat_d_C)+((1-df$att22)*(1-df$phat_d_C))

	df$phat_m_CD <- predict(Mmodel_1, type = "response")
	df$phat_M_denom <- (df$ever_unemp_age3539*df$phat_m_CD)+((1-df$ever_unemp_age3539)*(1-df$phat_m_CD))

	df$phat_d <- mean(df$att22)
	df$phat_D_num <- (df$att22*df$phat_d)+((1-df$att22)*(1-df$phat_d))

	Mmodel_2 <- glm(ever_unemp_age3539~att22, data=df, family=binomial("logit"))

	df$phat_m_D <- predict(Mmodel_2, type = "response")
	df$phat_M_num <- (df$ever_unemp_age3539*df$phat_m_D)+((1-df$ever_unemp_age3539)*(1-df$phat_m_D))

	df$sw4 <- (df$phat_M_num*df$phat_D_num)/(df$phat_M_denom*df$phat_D_denom)

	df[which(df[,"sw4"] > quantile(df[,"sw4"], probs=0.99, na.rm=T)), "sw4"] <- quantile(df[, "sw4"], probs=0.99, na.rm=T)
	df[which(df[,"sw4"] < quantile(df[,"sw4"], probs=0.01, na.rm=T)), "sw4"] <- quantile(df[, "sw4"], probs=0.01, na.rm=T)

	Ymodel <- lm(std_cesd_age40~att22*ever_unemp_age3539, data=df, weights=sw4)

	point.est <- Ymodel$coefficients["att22"]

	return(point.est)
	}

#setup parallel computing cluster
my.cluster <- parallel::makeCluster(ncores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("ipwmed", "ipwcde"), envir=environment())
registerDoRNG(3308004)

#compute bootstrap estimates
ipw.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {

	boot.data <- nlsy[sample(nrow(nlsy), nrow(nlsy), replace=TRUE),]

	boot.est <- ipwmed(boot.data)
	
	boot.est.cde <- ipwcde(boot.data)

	return(c(boot.est, boot.est.cde))
	}

stopCluster(my.cluster)
rm(my.cluster)

ipw.boot <- matrix(unlist(ipw.boot), ncol=4, byrow=TRUE)

#compute bootstrap 95% CIs
ipwest.output <- matrix(data=NA, nrow=4, ncol=3)

for (i in 1:nrow(ipwest.output)) {
	ipwest.output[i,1] <- round(quantile(ipw.boot[,i], prob=0.025), digits=3)
	ipwest.output[i,2] <- round(quantile(ipw.boot[,i], prob=0.975), digits=3)
	}

#compute bootstrap pvalues for null of no effect
ATE_null <- NDE_null <- NIE_null <- CDE_null <- 0

ipw.boot <- as.data.frame(ipw.boot)

ipw.boot <- 
	ipw.boot %>% 
		mutate(
			ATEpval = 2*min(mean(ifelse(V1>ATE_null, 1, 0)), mean(ifelse(V1<ATE_null, 1, 0))),
			NDEpval = 2*min(mean(ifelse(V2>NDE_null, 1, 0)), mean(ifelse(V2<NDE_null, 1, 0))),
			NIEpval = 2*min(mean(ifelse(V3>NIE_null, 1, 0)), mean(ifelse(V3<NIE_null, 1, 0))),
			CDEpval = 2*min(mean(ifelse(V4>CDE_null, 1, 0)), mean(ifelse(V4<CDE_null, 1, 0))))

for (i in 1:nrow(ipwest.output)) {
	ipwest.output[i,3] <- round(ipw.boot[1,i+4], digits=3)
	}

ipwest.output <- data.frame(ipwest.output, row.names=c("ATEhat^ipw", "NDEhat^ipw", "NIEhat^ipw", "CDE0hat^ipw"))
colnames(ipwest.output) <- c("ll.95ci", "ul.95ci", "pvalue")

print(ipwest.output)

sink()
