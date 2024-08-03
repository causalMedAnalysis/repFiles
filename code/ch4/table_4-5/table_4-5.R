###Table 4.5###

rm(list=ls())

packages<-c("dplyr", "tidyr", "devtools", "doParallel", "doRNG")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

#install_github("xiangzhou09/rwrmed")
library(rwrmed)

nboot <- 2000
ncores <- parallel::detectCores()-2

##office
datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch4/_LOGS/"

##home
#datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch4/_LOGS/"

set.seed(60657)

##input data
nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("cesd_age40", "ever_unemp_age3539", "faminc_adj_age3539", "log_faminc_adj_age3539", 
	"att22", "female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")]),]

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##RWR estimates

#specify form of models for L, M, and Y
Lmodel.x <- lm(ever_unemp_age3539~att22*(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3), data=nlsy)

Mform.x <- log_faminc_adj_age3539~att22*(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3)

Yform.x <- std_cesd_age40~(log_faminc_adj_age3539*att22)+(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3)*(log_faminc_adj_age3539+att22)+(ever_unemp_age3539*log_faminc_adj_age3539)

#fit RWR models
rwr.fit <- rwrmed(
	treatment="att22", 
	pre_cov=c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"), 
	zmodels=list(Lmodel.x),
	m_form=Mform.x,
	y_form=Yform.x,
	data=nlsy)

#compute point estimates 
rwr.decomp <- decomp(rwr.fit, a0=0, a1=1, m=log(50000), bootstrap=F)

#setup parallel computing cluster
my.cluster <- parallel::makeCluster(ncores, type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("Lmodel.x", "Mform.x", "Yform.x"), envir=environment())
clusterEvalQ(cl=my.cluster, library(rwrmed))
registerDoRNG(3308004)

#compute bootstrap estimates
rwrmed.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {

	boot.data <- nlsy[sample(nrow(nlsy), nrow(nlsy), replace=TRUE),]

	Lmodel.boot <- lm(ever_unemp_age3539~att22*(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3), data=boot.data)

	boot.rwr.fit <- rwrmed(
				treatment="att22", 
				pre_cov=c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"), 
				zmodels=list(Lmodel.boot),
				m_form=Mform.x,
				y_form=Yform.x,
				data=boot.data)

	boot.rwr.decomp <- decomp(boot.rwr.fit, a0=0, a1=1, m=log(50000), bootstrap=F)

	IDE <- boot.rwr.decomp$twocomp[1,1]
	IIE <- boot.rwr.decomp$twocomp[2,1]
	OE <- boot.rwr.decomp$twocomp[3,1]
	CDE <- boot.rwr.decomp$fourcomp[1,1]

	return(list(IDE, IIE, OE, CDE))
	}

stopCluster(my.cluster)
rm(my.cluster)

rwrmed.boot <- matrix(unlist(rwrmed.boot), ncol=4, byrow=TRUE)

#compute 95% CIs
rwr.output <- matrix(data=NA, nrow=4, ncol=4)

rwr.output[,1] <- c(rwr.decomp$twocomp[,1], rwr.decomp$fourcomp[1,1])

for (i in 1:nrow(rwr.output)) {

	rwr.output[i,1] <- round(rwr.output[i,1], digits=3)

	rwr.output[i,2] <- round(quantile(rwrmed.boot[,i], prob=0.025), digits=3)
	rwr.output[i,3] <- round(quantile(rwrmed.boot[,i], prob=0.975), digits=3)
	}

#compute bootstrap pvalues for null of no effect
IDE_null <- IIE_null <- OE_null <- CDE_null <- 0

rwrmed.boot <- as.data.frame(rwrmed.boot)

rwrmed.boot <- 
	rwrmed.boot %>% 
		mutate(
			IDEpval = 2*min(mean(ifelse(V1>IDE_null, 1, 0)), mean(ifelse(V1<IDE_null, 1, 0))),
			IIEpval = 2*min(mean(ifelse(V2>IIE_null, 1, 0)), mean(ifelse(V2<IIE_null, 1, 0))),
			OEpval = 2*min(mean(ifelse(V3>OE_null, 1, 0)), mean(ifelse(V3<OE_null, 1, 0))),
			CDEpval = 2*min(mean(ifelse(V4>CDE_null, 1, 0)), mean(ifelse(V4<CDE_null, 1, 0))))

for (i in 1:nrow(rwr.output)) {
	rwr.output[i,4] <- round(rwrmed.boot[1,i+4], digits=3)
	}

rownames(rwr.output)<-c("IDE", "IIE", "OE", "CDE")
colnames(rwr.output) <- c("point.est", "ll.95ci", "ul.95ci", "pval")

##simulation estimates

#define function for simulation estimators
medsim <- function(data, num.sim=2000, Lform, Mform, Yform) {
	
	df <- data

	Lmodel <- glm(Lform, data=df, family=binomial("logit"))
  
	Mmodel <- lm(Mform, data=df)

	Ymodel <- lm(Yform, data=df)

	idata <- df

	Y0L0M0 <- Y1L1M1 <- Y1L1M0 <- Y0L0m <- Y1L1m <- numeric(num.sim)

	for (i in 1:num.sim) {

		idata$att22 <- 0

		phat_L0 <- predict(Lmodel, newdata=idata, type="response")
		yhat_M0 <- predict(Mmodel, newdata=idata, type="response")

		L0 <- rbinom(nrow(idata), size=1, prob=phat_L0)
		M0 <- rnorm(nrow(idata), yhat_M0, sd=sigma(Mmodel))

		idata$att22 <- 1

		phat_L1 <- predict(Lmodel, newdata=idata, type="response")
		yhat_M1 <- predict(Mmodel, newdata=idata, type="response")

		L1 <- rbinom(nrow(idata), size=1, prob=phat_L1)
		M1 <- rnorm(nrow(idata), yhat_M1, sd=sigma(Mmodel))

		idata$ever_unemp_age3539 <- L1
		idata$log_faminc_adj_age3539 <- M1

		yhat_Y1L1M1 <- predict(Ymodel, newdata=idata, type="response")
		Y1L1M1[i] <- mean(rnorm(nrow(idata), yhat_Y1L1M1, sd=sigma(Ymodel)))

		idata$log_faminc_adj_age3539 <- M0

		yhat_Y1L1M0 <- predict(Ymodel, newdata=idata, type="response")
		Y1L1M0[i] <- mean(rnorm(nrow(idata), yhat_Y1L1M0, sd=sigma(Ymodel)))

		idata$att22 <- 0
		idata$ever_unemp_age3539 <- L0

		yhat_Y0L0M0 <- predict(Ymodel, newdata=idata, type="response")
		Y0L0M0[i] <- mean(rnorm(nrow(idata), yhat_Y0L0M0, sd=sigma(Ymodel)))

		idata$log_faminc_adj_age3539 <- log(50000)

		yhat_Y0L0m <- predict(Ymodel, newdata=idata, type="response")
		Y0L0m[i] <- mean(rnorm(nrow(idata), yhat_Y0L0m, sd=sigma(Ymodel)))

		idata$att22 <- 1
		idata$ever_unemp_age3539 <- L1

		yhat_Y1L1m <- predict(Ymodel, newdata=idata, type="response")
		Y1L1m[i] <- mean(rnorm(nrow(idata), yhat_Y1L1m, sd=sigma(Ymodel)))
		}

	IDE <- mean(Y1L1M0) - mean(Y0L0M0)  
	IIE <- mean(Y1L1M1) - mean(Y1L1M0)
	OE <- mean(Y1L1M1) - mean(Y0L0M0)
	CDE <- mean(Y1L1m) - mean(Y0L0m)

	point.est <- list(IDE, IIE, OE, CDE)

	return(point.est)
	}

#specify form of models for L, M, and Y
Lform.x <- ever_unemp_age3539 ~ att22 * (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

Mform.x <- log_faminc_adj_age3539 ~ att22 * (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

Yform.x <- std_cesd_age40 ~ (log_faminc_adj_age3539 * att22) + (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3) * (log_faminc_adj_age3539 + att22) + (ever_unemp_age3539 * log_faminc_adj_age3539)

#compute point estimates
medsim.est <- medsim(data=nlsy, num.sim=2000, Lform=Lform.x, Mform=Mform.x, Yform=Yform.x)
medsim.est <- matrix(unlist(medsim.est), ncol=4, byrow=TRUE)

#setup parallel computing cluster
my.cluster <- parallel::makeCluster(ncores, type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("medsim", "Lform.x", "Mform.x", "Yform.x"), envir=environment())
registerDoRNG(3308004)

#compute bootstrap estimates
medsim.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {

	boot.data <- nlsy[sample(nrow(nlsy), nrow(nlsy), replace=TRUE),]

	boot.est <- medsim(boot.data, num.sim=2000, Lform=Lform.x, Mform=Mform.x, Yform=Yform.x)
	
	return(boot.est)
	}

stopCluster(my.cluster)
rm(my.cluster)

medsim.boot <- matrix(unlist(medsim.boot), ncol=4, byrow=TRUE)

#compute 95% CIs
medsim.output <- matrix(data=NA, nrow=4, ncol=4)

for (i in 1:nrow(medsim.output)) {

	medsim.output[i,1] <- round(medsim.est[i], digits=3)

	medsim.output[i,2] <- round(quantile(medsim.boot[,i], prob=0.025), digits=3)
	medsim.output[i,3] <- round(quantile(medsim.boot[,i], prob=0.975), digits=3)
	}

#compute bootstrap pvalues for null of no effect
IDE_null <- IIE_null <- OE_null <- CDE_null <- 0

medsim.boot <- as.data.frame(medsim.boot)

medsim.boot <- 
	medsim.boot %>% 
		mutate(
			IDEpval = 2*min(mean(ifelse(V1>IDE_null, 1, 0)), mean(ifelse(V1<IDE_null, 1, 0))),
			IIEpval = 2*min(mean(ifelse(V2>IIE_null, 1, 0)), mean(ifelse(V2<IIE_null, 1, 0))),
			OEpval = 2*min(mean(ifelse(V3>OE_null, 1, 0)), mean(ifelse(V3<OE_null, 1, 0))),
			CDEpval = 2*min(mean(ifelse(V4>CDE_null, 1, 0)), mean(ifelse(V4<CDE_null, 1, 0))))

for (i in 1:nrow(medsim.output)) {
	medsim.output[i,4] <- round(medsim.boot[1,i+4], digits=3)
	}

medsim.output <- data.frame(medsim.output, row.names=c("IDE", "IIE", "OE", "CDE"))
colnames(medsim.output) <- c("point.est", "ll.95ci", "ul.95ci", "pval")

##ipw estimates

##weighting estimators

#define function for weighting estimators
ipwmed <- function(data, Dform, Lform, Mform) {

	df <- data

	df$id <- 1:nrow(df)

	Dmodel <- glm(Dform, data=df, family=binomial("logit"))

	Lmodel <- glm(Lform, data=df, family=binomial("logit"))
  
	Mmodel <- lm(Mform, data=df)

	idataD0 <- df %>% mutate(att22 = 0)
	idataD1 <- df %>% mutate(att22 = 1)

	idataD0L0 <- df %>% mutate(att22 = 0, ever_unemp_age3539 = 0)
	idataD1L0 <- df %>% mutate(att22 = 1, ever_unemp_age3539 = 0)
	idataD0L1 <- df %>% mutate(att22 = 0, ever_unemp_age3539 = 1)
	idataD1L1 <- df %>% mutate(att22 = 1, ever_unemp_age3539 = 1)

	phatD_C <- df %>% 
		mutate(
			pD1_C = predict(Dmodel, newdata=df, type = "response"),
			pD0_C = 1-pD1_C,
			pD1 = mean(att22),
			pD0 = 1-pD1) %>%
		select(id, pD1_C, pD0_C, pD1, pD0)

	phatL_D1C <- idataD1 %>% 
		mutate(
			pL1_D1C = predict(Lmodel, newdata=idataD1, type = "response"),
			pL0_D1C = 1-pL1_D1C) %>%
		select(id, pL1_D1C, pL0_D1C)

	phatL_D0C <- idataD0 %>% 
		mutate(
			pL1_D0C = predict(Lmodel, newdata=idataD0, type = "response"),
			pL0_D0C = 1-pL1_D0C) %>%
		select(id, pL1_D0C, pL0_D0C)

	phatM_D <- df %>%
		mutate(
			pM_D1 = case_when(att22 == 1 ~ dnorm(log_faminc_adj_age3539, mean(df$log_faminc_adj_age3539[df$att==1]), sd=sd(df$log_faminc_adj_age3539-(mean(df$log_faminc_adj_age3539[df$att==1])*df$att22)-(mean(df$log_faminc_adj_age3539[df$att==0])*(1-df$att22))))),
			pM_D0 = case_when(att22 == 0 ~ dnorm(log_faminc_adj_age3539, mean(df$log_faminc_adj_age3539[df$att==0]), sd=sd(df$log_faminc_adj_age3539-(mean(df$log_faminc_adj_age3539[df$att==1])*df$att22)-(mean(df$log_faminc_adj_age3539[df$att==0])*(1-df$att22)))))) %>%
		select(id, pM_D1, pM_D0) 

	phatM_D1LC <- idataD1 %>% 
		mutate(
			pM_D1LC = dnorm(log_faminc_adj_age3539, predict(Mmodel, newdata=idataD1, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D1LC)

	phatM_D0LC <- idataD1 %>% 
		mutate(pM_D0LC = dnorm(log_faminc_adj_age3539, predict(Mmodel, newdata=idataD0, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D0LC)

	phatM_D0L0C <- idataD0L0 %>% 
		mutate(pM_D0L0C = dnorm(log_faminc_adj_age3539, predict(Mmodel, newdata=idataD0L0, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D0L0C)

	phatM_D1L0C <- idataD1L0 %>% 
		mutate(pM_D1L0C = dnorm(log_faminc_adj_age3539, predict(Mmodel, newdata=idataD1L0, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D1L0C)

	phatM_D0L1C <- idataD0L1 %>% 
		mutate(pM_D0L1C = dnorm(log_faminc_adj_age3539, predict(Mmodel, newdata=idataD0L1, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D0L1C)

	phatM_D1L1C <- idataD1L1 %>% 
		mutate(pM_D1L1C = dnorm(log_faminc_adj_age3539, predict(Mmodel, newdata=idataD1L1, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D1L1C)

	df.wts <- df %>%
	    full_join(phatD_C, by = "id") %>%
	    full_join(phatL_D1C, by = "id") %>%
	    full_join(phatL_D0C, by = "id") %>%
	    full_join(phatM_D, by = "id") %>%
	    full_join(phatM_D1LC, by = "id") %>%
	    full_join(phatM_D0LC, by = "id") %>%
	    full_join(phatM_D0L0C, by = "id") %>%
	    full_join(phatM_D1L0C, by = "id") %>%
	    full_join(phatM_D0L1C, by = "id") %>%
	    full_join(phatM_D1L1C, by = "id")

	df.wts <- df.wts %>%
		mutate(
			sw1 = case_when(att22 == 0 ~ (pD0/pD0_C) * (1/pM_D0LC) * ((pM_D0L0C * pL0_D0C) + (pM_D0L1C * pL1_D0C))),
			sw2 = case_when(att22 == 1 ~ (pD1/pD1_C) * (1/pM_D1LC) * ((pM_D1L0C * pL0_D1C) + (pM_D1L1C * pL1_D1C))),
			sw3 = case_when(att22 == 1 ~ (pD1/pD1_C) * (1/pM_D1LC) * ((pM_D0L0C * pL0_D0C) + (pM_D0L1C * pL1_D0C))),
			sw4 = case_when(
				att22 == 1 ~ (pD1/pD1_C) * (pM_D1/pM_D1LC),
				att22 == 0 ~ (pD0/pD0_C) * (pM_D0/pM_D0LC))) %>%
		mutate(
			across(c(sw1, sw2, sw3, sw4), 
			~ifelse(. < quantile(., 0.01, na.rm = TRUE), quantile(., 0.01, na.rm = TRUE),
                   ifelse(. > quantile(., 0.99, na.rm = TRUE), quantile(., 0.99, na.rm = TRUE), .))))

	Ehat_Y0M0 <- weighted.mean(df.wts$std_cesd_age40[df.wts$att22==0], df.wts$sw1[df$att22==0])
	Ehat_Y1M1 <- weighted.mean(df.wts$std_cesd_age40[df.wts$att22==1], df.wts$sw2[df$att22==1])
	Ehat_Y1M0 <- weighted.mean(df.wts$std_cesd_age40[df.wts$att22==1], df.wts$sw3[df$att22==1])

	IDE <- Ehat_Y1M0-Ehat_Y0M0
	IIE <- Ehat_Y1M1-Ehat_Y1M0
	OE <- Ehat_Y1M1-Ehat_Y0M0

	Ymodel.wtd <- lm(std_cesd_age40 ~ att22 * log_faminc_adj_age3539, data=df.wts, weights=sw4)

	CDE <- Ymodel.wtd$coefficients["att22"] + log(50000)*Ymodel.wtd$coefficients["att22:log_faminc_adj_age3539"]

	##### NEED TO ADD CDE ESTIMATES #####

	point.est <- list(IDE, IIE, OE, CDE)

	return(point.est)
	}

#specify form of models for D, L, and M
Dform.x <- att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

Lform.x <- ever_unemp_age3539 ~ att22 * (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

Mform.x <- log_faminc_adj_age3539 ~ att22 * (ever_unemp_age3539 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

#compute point estimates
ipwmed.est <- ipwmed(data=nlsy, Dform=Dform.x, Lform=Lform.x, Mform=Mform.x)
ipwmed.est <- matrix(unlist(ipwmed.est), ncol=4, byrow=TRUE)

#setup parallel computing cluster
my.cluster <- parallel::makeCluster(ncores, type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("medsim", "Lform.x", "Mform.x", "Yform.x"), envir=environment())
clusterEvalQ(cl=my.cluster, library(dplyr))
registerDoRNG(3308004)

#compute bootstrap estimates
ipwmed.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {

	boot.data <- nlsy[sample(nrow(nlsy), nrow(nlsy), replace=TRUE),]

	boot.est <- ipwmed(data=boot.data, Dform=Dform.x, Lform=Lform.x, Mform=Mform.x)
	
	return(boot.est)
	}

stopCluster(my.cluster)
rm(my.cluster)

ipwmed.boot <- matrix(unlist(ipwmed.boot), ncol=4, byrow=TRUE)

#compute 95% CIs
ipwmed.output <- matrix(data=NA, nrow=4, ncol=4)

for (i in 1:nrow(ipwmed.output)) {

	ipwmed.output[i,1] <- round(ipwmed.est[i], digits=3)

	ipwmed.output[i,2] <- round(quantile(ipwmed.boot[,i], prob=0.025), digits=3)
	ipwmed.output[i,3] <- round(quantile(ipwmed.boot[,i], prob=0.975), digits=3)
	}

#compute bootstrap pvalues for null of no effect
IDE_null <- IIE_null <- OE_null <- CDE_null <- 0

ipwmed.boot <- as.data.frame(ipwmed.boot)

ipwmed.boot <- 
	ipwmed.boot %>% 
		mutate(
			IDEpval = 2*min(mean(ifelse(V1>IDE_null, 1, 0)), mean(ifelse(V1<IDE_null, 1, 0))),
			IIEpval = 2*min(mean(ifelse(V2>IIE_null, 1, 0)), mean(ifelse(V2<IIE_null, 1, 0))),
			OEpval = 2*min(mean(ifelse(V3>OE_null, 1, 0)), mean(ifelse(V3<OE_null, 1, 0))),
			CDEpval = 2*min(mean(ifelse(V4>CDE_null, 1, 0)), mean(ifelse(V4<CDE_null, 1, 0))))

for (i in 1:nrow(ipwmed.output)) {
	ipwmed.output[i,4] <- round(ipwmed.boot[1,i+4], digits=3)
	}

ipwmed.output <- data.frame(ipwmed.output, row.names=c("IDE", "IIE", "OE", "CDE"))
colnames(ipwmed.output) <- c("point.est", "ll.95ci", "ul.95ci", "pval")

##print output
sink(paste(logdir, "table_4-5_log.txt", sep=""))

cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("RWR Estimates\n")
print(rwr.output)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Simulation Estimates\n")
print(medsim.output)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("IPW Estimates\n")
print(ipwmed.output)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

sink()
