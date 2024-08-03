###Table 3.4###

rm(list=ls())

packages<-c("dplyr", "tidyr")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

##office
#datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

##home
datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"

sink(paste(logdir, "table_3-4_log.txt", sep=""))

##input data

nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("cesd_age40", "ever_unemp_age3539", "att22", "female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")]),]

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##inverse probability weighting estimators

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

ipwmed.est <- ipwmed(nlsy)
ipwmed.est <- matrix(unlist(ipwmed.est), ncol=3, byrow=TRUE)

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

ipwcde.est <- ipwcde(nlsy)

ipwest.output <- matrix(data=NA, nrow=4, ncol=1)

ipwest.output[,1] <- t(c(round(ipwmed.est[,1], digits=3), round(ipwmed.est[,2], digits=3), round(ipwmed.est[,3], digits=3), round(ipwcde.est, digits=3)))

ipwest.output <- data.frame(ipwest.output, row.names=c("ATEhat^ipw", "NDEhat^ipw", "NIEhat^ipw", "CDE0hat^ipw"))
colnames(ipwest.output) <- c("point.est")

print(ipwest.output)

##note: minor differences compared with stata due to differences in the numerical precision with which predicted probabilities are stored

sink()
