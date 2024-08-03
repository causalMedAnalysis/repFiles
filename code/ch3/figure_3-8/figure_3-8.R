###Table 3.5###

rm(list=ls())

packages<-c("dplyr", "tidyr", "foreach", "doParallel", "doRNG", "ggplot2")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

nboot <- 2000

ncores <- parallel::detectCores()-1

##office
#datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"
#figdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/figures/ch3/"

##home
datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"
figdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/figures/ch3/"

sink(paste(logdir, "figure_3-8_log.txt", sep=""))

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

#setup parallel computing cluster
my.cluster <- parallel::makeCluster(ncores,type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("ipwmed"), envir=environment())
registerDoRNG(3308004)

#compute bootstrap estimates
ipw.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {
	boot.data <- nlsy[sample(nrow(nlsy), nrow(nlsy), replace=TRUE),]

	boot.est <- ipwmed(boot.data)
	
	return(boot.est)
	}

stopCluster(my.cluster)
rm(my.cluster)

ipw.boot <- as.data.frame(matrix(unlist(ipw.boot), ncol=3, byrow=TRUE))

dens <- ggplot(ipw.boot, aes(x=V3)) +
		geom_density(adjust=0.9, alpha=0.1, fill="black") +
	 	theme_bw() +
		scale_y_continuous(name="Density", limits=c(0, 60), breaks=seq(0, 60, 10)) +
		scale_x_continuous(name=expression(hat(theta)[b]), limits=c(-0.05, 0.02), breaks=round(seq(-0.05, 0.02, 0.01), 2)) 

dens.data <- ggplot_build(dens)$data[[1]]

dens <- dens + geom_area(data=subset(dens.data, x>0), aes(x=x, y=y), fill="grey30") 

print(dens)

ggsave(paste(figdir, "figure_3-8.png", sep=""), height=4.5, width=4.5, units="in", dpi=600)

sink()
