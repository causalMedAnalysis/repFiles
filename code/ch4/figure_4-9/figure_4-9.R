###Figure 4.9###

rm(list=ls())

packages<-c("dplyr", "tidyr", "foreign", "VGAM", "devtools", "doParallel", "doRNG", "ggplot2")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

#install_github("xiangzhou09/rwrmed")

library(rwrmed)

nboot <- 2000
nsim <- 2000
ncores <- parallel::detectCores()-2

##office
datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch4/_LOGS/"
figdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/figures/ch4/"

##home
#datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch4/_LOGS/"
#figdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/figures/ch4/"

set.seed(60657)

##input data
plowData <- read.dta(paste(datadir, "plowUse/plowUse.dta", sep=""))

allvars <- c("women_politics", "plow", "ln_income", "agricultural_suitability", 
             "tropical_climate", "large_animals", "rugged", "polity2_2000")

plowData <- na.omit(plowData[, c('isocode', allvars)])

plowData$plow <- round(plowData$plow)

plowData$women_politics <- plowData$women_politics/100

plowData$authGovCat <- ifelse(plowData$polity2_2000 >= -10 & plowData$polity2_2000 <= -6, 1,
                          ifelse(plowData$polity2_2000 >= -5 & plowData$polity2_2000 <= 5, 2,
                                 ifelse(plowData$polity2_2000 >= 6 & plowData$polity2_2000 <= 10, 3, NA)))

##RWR estimates

#define models
Lmodel.rwr <- lm(authGovCat~plow+agricultural_suitability+tropical_climate+large_animals+rugged, data=plowData)

Mform.rwr <- ln_income~plow+agricultural_suitability+tropical_climate+large_animals+rugged

Yform.rwr <- women_politics~(plow*ln_income)+authGovCat+agricultural_suitability+tropical_climate+large_animals+rugged

#fit models
rwr.fit <- rwrmed(
	treatment="plow", 
	pre_cov=c("agricultural_suitability", "tropical_climate", "large_animals", "rugged"), 
	zmodels=list(Lmodel.rwr),
	m_form=Mform.rwr,
	y_form=Yform.rwr,
	data=plowData)

#compute point and interval estimates 
rwr.decomp <- decomp(rwr.fit, a0=0, a1=1, m=7.5, bootstrap=T, rep=nboot)

rwr.output <- data.frame(
	point.est = c(rwr.decomp$twocomp[,1], rwr.decomp$fourcomp[1,1]),
	ll.95ci = c(rwr.decomp$twocomp[,3], rwr.decomp$fourcomp[1,3]),
	ul.95ci = c(rwr.decomp$twocomp[,4], rwr.decomp$fourcomp[1,4]))

rwr.output <- data.frame(lapply(rwr.output, function(x) round(x, 3)))

rwr.output <- data.frame(estimand = c("IDE(1,0)", "IIE(1,0)", "OE(1,0)", "CDE(1,0,1.8K)"), rwr.output)

##simulation estimates

#define function for simulation estimators
medsim <- function(data, num.sim=1000, Lform, Mform, Yform) {
	
	df <- data

      Lmodel <- vgam(Lform, family=cumulative(parallel=TRUE, reverse=FALSE), data=df)

	Mmodel <- lm(Mform, data=df)

	Ymodel <- glm(Yform, data=df, family=quasibinomial(link="log"))

	idata <- df

	Y0L0M0 <- Y1L1M1 <- Y1L1M0 <- Y0L0m <- Y1L1m <- numeric(num.sim)

	for (i in 1:num.sim) {

		idata$plow <- 0

		phat_L0 <- predict(Lmodel, newdata=idata, type="response")

		phat1_L0 <- phat_L0[,1]
		phat2_L0 <- phat_L0[,2] 
		phat3_L0 <- phat_L0[,3]

		runif <- runif(nrow(idata))
		L0 <- ifelse(runif <= phat1_L0, 1,
				ifelse(runif > phat1_L0 & runif <= (phat1_L0 + phat2_L0), 2,
					ifelse(runif > (phat1_L0 + phat2_L0), 3, NA)))

		yhat_M0 <- predict(Mmodel, newdata=idata, type="response")

		M0 <- rnorm(nrow(idata), yhat_M0, sd=sigma(Mmodel))

		idata$plow <- 1

		phat_L1 <- predict(Lmodel, newdata=idata, type="response")

		phat1_L1 <- phat_L1[,1]
		phat2_L1 <- phat_L1[,2] 
		phat3_L1 <- phat_L1[,3]

		runif <- runif(nrow(idata))
		L1 <- ifelse(runif <= phat1_L1, 1,
				ifelse(runif > phat1_L1 & runif <= (phat1_L1 + phat2_L1), 2,
					ifelse(runif > (phat1_L1 + phat2_L1), 3, NA)))

		yhat_M1 <- predict(Mmodel, newdata=idata, type="response")

		M1 <- rnorm(nrow(idata), yhat_M1, sd=sigma(Mmodel))

		idata$authGovCat <- L1
		idata$ln_income <- M1

		phat_Y1L1M1 <- predict(Ymodel, newdata=idata, type="response")
		Y1L1M1[i] <- mean(rbinom(n=nrow(idata), size=100, prob=phat_Y1L1M1)/100)

		idata$ln_income <- M0

		phat_Y1L1M0 <- predict(Ymodel, newdata=idata, type="response")
		Y1L1M0[i] <- mean(rbinom(n=nrow(idata), size=100, prob=phat_Y1L1M0)/100)

		idata$plow <- 0
		idata$authGovCat <- L0

		phat_Y0L0M0 <- predict(Ymodel, newdata=idata, type="response")
		Y0L0M0[i] <- mean(rbinom(n=nrow(idata), size=100, prob=phat_Y0L0M0)/100)

		idata$ln_income <- 7.5

		phat_Y0L0m <- predict(Ymodel, newdata=idata, type="response")
		Y0L0m[i] <- mean(rbinom(n=nrow(idata), size=100, prob=phat_Y0L0m)/100)

		idata$plow <- 1
		idata$authGovCat <- L1

		phat_Y1L1m <- predict(Ymodel, newdata=idata, type="response")
		Y1L1m[i] <- mean(rbinom(n=nrow(idata), size=100, prob=phat_Y1L1m)/100)
		}

	IDE <- mean(Y1L1M0) - mean(Y0L0M0)  
	IIE <- mean(Y1L1M1) - mean(Y1L1M0)
	OE <- mean(Y1L1M1) - mean(Y0L0M0)
	CDE <- mean(Y1L1m) - mean(Y0L0m)

	point.est <- list(IDE, IIE, OE, CDE)

	return(point.est)
	}

#specify form of models for L, M, and Y
Lform.sim <- authGovCat~plow+agricultural_suitability+tropical_climate+large_animals+rugged

Mform.sim <- ln_income~plow+agricultural_suitability+tropical_climate+large_animals+rugged

Yform.sim <- women_politics~(plow*ln_income)+authGovCat+agricultural_suitability+tropical_climate+large_animals+rugged

#compute point estimates
medsim.est <- medsim(data=plowData, num.sim=nsim, Lform=Lform.sim, Mform=Mform.sim, Yform=Yform.sim)
medsim.est <- matrix(unlist(medsim.est), ncol=4, byrow=TRUE)

#setup parallel computing cluster
my.cluster <- parallel::makeCluster(ncores, type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("medsim", "Lform.sim", "Mform.sim", "Yform.sim"), envir=environment())
clusterEvalQ(my.cluster, {
  library(VGAM)
  library(dplyr)
  library(tidyr)
  })
registerDoRNG(3308004)

#compute bootstrap estimates
medsim.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {

	boot.data <- plowData[sample(nrow(plowData), nrow(plowData), replace=TRUE),]

	boot.est <- medsim(boot.data, num.sim=nsim, Lform=Lform.sim, Mform=Mform.sim, Yform=Yform.sim)
	
	return(boot.est)
	}

stopCluster(my.cluster)
rm(my.cluster)

medsim.boot <- matrix(unlist(medsim.boot), ncol=4, byrow=TRUE)

#compute 95% CIs
medsim.output <- matrix(data=NA, nrow=4, ncol=3)

for (i in 1:nrow(medsim.output)) {

	medsim.output[i,1] <- round(medsim.est[i], digits=3)
	medsim.output[i,2] <- round(quantile(medsim.boot[,i], prob=0.025, na.rm=T), digits=3)
	medsim.output[i,3] <- round(quantile(medsim.boot[,i], prob=0.975, na.rm=T), digits=3)
	}

medsim.output <- data.frame(estimand = c("IDE(1,0)", "IIE(1,0)", "OE(1,0)", "CDE(1,0,1.8K)"), medsim.output)
colnames(medsim.output) <- c("estimand", "point.est", "ll.95ci", "ul.95ci")

##ipw estimates

##weighting estimators

#define function for weighting estimators
ipwmed <- function(data, Dform, Lform, Mform) {

	df <- data

	df$id <- 1:nrow(df)

	Dmodel <- glm(Dform, data=df, family=binomial("logit"))

      Lmodel <- vgam(Lform, family=cumulative(parallel=TRUE, reverse=FALSE), data=df)

	Mmodel <- lm(Mform, data=df)

	idataD0 <- df %>% mutate(plow = 0)
	idataD1 <- df %>% mutate(plow = 1)

	idataD0L1 <- df %>% mutate(plow = 0, authGovCat = 1)
	idataD1L1 <- df %>% mutate(plow = 1, authGovCat = 1)
	idataD0L2 <- df %>% mutate(plow = 0, authGovCat = 2)
	idataD1L2 <- df %>% mutate(plow = 1, authGovCat = 2)
	idataD0L3 <- df %>% mutate(plow = 0, authGovCat = 3)
	idataD1L3 <- df %>% mutate(plow = 1, authGovCat = 3)

	phatD_C <- df %>% 
		mutate(
			pD1_C = predict(Dmodel, newdata=df, type = "response"),
			pD0_C = 1-pD1_C,
			pD1 = mean(plow),
			pD0 = 1-pD1) %>%
		select(id, pD1_C, pD0_C, pD1, pD0)

	phatL_D1C <- idataD1 %>% 
		mutate(
			pL1_D1C = predict(Lmodel, newdata=idataD1, type = "response")[,1],
			pL2_D1C = predict(Lmodel, newdata=idataD1, type = "response")[,2],
			pL3_D1C = predict(Lmodel, newdata=idataD1, type = "response")[,3]) %>%
		select(id, pL1_D1C, pL2_D1C, pL3_D1C)

	phatL_D0C <- idataD0 %>% 
		mutate(
			pL1_D0C = predict(Lmodel, newdata=idataD1, type = "response")[,1],
			pL2_D0C = predict(Lmodel, newdata=idataD1, type = "response")[,2],
			pL3_D0C = predict(Lmodel, newdata=idataD1, type = "response")[,3]) %>%
		select(id, pL1_D0C, pL2_D0C, pL3_D0C)

	phatM_D <- df %>%
		mutate(
			pM_D1 = case_when(plow == 1 ~ dnorm(ln_income, mean(df$ln_income[df$plow==1]), sd=sd(df$ln_income-(mean(df$ln_income[df$plow==1])*df$plow)-(mean(df$ln_income[df$plow==0])*(1-df$plow))))),
			pM_D0 = case_when(plow == 0 ~ dnorm(ln_income, mean(df$ln_income[df$plow==0]), sd=sd(df$ln_income-(mean(df$ln_income[df$plow==1])*df$plow)-(mean(df$ln_income[df$plow==0])*(1-df$plow)))))) %>%
		select(id, pM_D1, pM_D0) 

	phatM_D1LC <- idataD1 %>% 
		mutate(
			pM_D1LC = dnorm(ln_income, predict(Mmodel, newdata=idataD1, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D1LC)

	phatM_D0LC <- idataD1 %>% 
		mutate(pM_D0LC = dnorm(ln_income, predict(Mmodel, newdata=idataD0, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D0LC)

	phatM_D0L1C <- idataD0L1 %>% 
		mutate(pM_D0L1C = dnorm(ln_income, predict(Mmodel, newdata=idataD0L1, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D0L1C)

	phatM_D0L2C <- idataD0L2 %>% 
		mutate(pM_D0L2C = dnorm(ln_income, predict(Mmodel, newdata=idataD0L2, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D0L2C)

	phatM_D0L3C <- idataD0L3 %>% 
		mutate(pM_D0L3C = dnorm(ln_income, predict(Mmodel, newdata=idataD0L3, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D0L3C)

	phatM_D1L1C <- idataD1L1 %>% 
		mutate(pM_D1L1C = dnorm(ln_income, predict(Mmodel, newdata=idataD1L1, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D1L1C)

	phatM_D1L2C <- idataD1L2 %>% 
		mutate(pM_D1L2C = dnorm(ln_income, predict(Mmodel, newdata=idataD1L2, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D1L2C)

	phatM_D1L3C <- idataD1L3 %>% 
		mutate(pM_D1L3C = dnorm(ln_income, predict(Mmodel, newdata=idataD1L3, type = "response"), sd=sigma(Mmodel))) %>%
		select(id, pM_D1L3C)

	df.wts <- df %>%
	    full_join(phatD_C, by = "id") %>%
	    full_join(phatL_D1C, by = "id") %>%
	    full_join(phatL_D0C, by = "id") %>%
	    full_join(phatM_D, by = "id") %>%
	    full_join(phatM_D1LC, by = "id") %>%
	    full_join(phatM_D0LC, by = "id") %>%
	    full_join(phatM_D0L1C, by = "id") %>%
	    full_join(phatM_D0L2C, by = "id") %>%
	    full_join(phatM_D0L3C, by = "id") %>%
	    full_join(phatM_D1L1C, by = "id") %>%
	    full_join(phatM_D1L2C, by = "id") %>%
	    full_join(phatM_D1L3C, by = "id")

	df.wts <- df.wts %>%
		mutate(
			sw1 = case_when(plow == 0 ~ (pD0/pD0_C) * (1/pM_D0LC) * ((pM_D0L1C * pL1_D0C) + (pM_D0L2C * pL2_D0C) + (pM_D0L3C * pL3_D0C))),
			sw2 = case_when(plow == 1 ~ (pD1/pD1_C) * (1/pM_D1LC) * ((pM_D1L1C * pL1_D1C) + (pM_D1L2C * pL2_D1C) + (pM_D1L3C * pL3_D1C))),
			sw3 = case_when(plow == 1 ~ (pD1/pD1_C) * (1/pM_D1LC) * ((pM_D0L1C * pL1_D0C) + (pM_D0L2C * pL2_D0C) + (pM_D0L3C * pL3_D0C))),
			sw4 = case_when(
				plow == 1 ~ (pD1/pD1_C) * (pM_D1/pM_D1LC),
				plow == 0 ~ (pD0/pD0_C) * (pM_D0/pM_D0LC))) %>%
		mutate(
			across(c(sw1, sw2, sw3, sw4), 
			~ifelse(. < quantile(., 0.02, na.rm = TRUE), quantile(., 0.02, na.rm = TRUE),
                   ifelse(. > quantile(., 0.98, na.rm = TRUE), quantile(., 0.98, na.rm = TRUE), .))))

	Ehat_Y0M0 <- weighted.mean(df.wts$women_politics[df.wts$plow==0], df.wts$sw1[df$plow==0])
	Ehat_Y1M1 <- weighted.mean(df.wts$women_politics[df.wts$plow==1], df.wts$sw2[df$plow==1])
	Ehat_Y1M0 <- weighted.mean(df.wts$women_politics[df.wts$plow==1], df.wts$sw3[df$plow==1])

	IDE <- Ehat_Y1M0-Ehat_Y0M0
	IIE <- Ehat_Y1M1-Ehat_Y1M0
	OE <- Ehat_Y1M1-Ehat_Y0M0

	Ymodel.wtd <- lm(women_politics ~ plow * ln_income, data=df.wts, weights=sw4)

	CDE <- Ymodel.wtd$coefficients["plow"] + 7.5 * Ymodel.wtd$coefficients["plow:ln_income"]

	point.est <- list(IDE, IIE, OE, CDE)

	return(point.est)
	}

#specify form of models for D, L, and M
Dform.ipw <- plow~agricultural_suitability+tropical_climate+large_animals+rugged

Lform.ipw <- authGovCat~plow+agricultural_suitability+tropical_climate+large_animals+rugged

Mform.ipw <- ln_income~authGovCat+plow+agricultural_suitability+tropical_climate+large_animals+rugged

#compute point estimates
ipwmed.est <- ipwmed(data=plowData, Dform=Dform.ipw, Lform=Lform.ipw, Mform=Mform.ipw)
ipwmed.est <- matrix(unlist(ipwmed.est), ncol=4, byrow=TRUE)

#setup parallel computing cluster
my.cluster <- parallel::makeCluster(ncores, type="PSOCK")
doParallel::registerDoParallel(cl=my.cluster)
clusterExport(cl=my.cluster, list("ipwmed", "Dform.ipw", "Lform.ipw", "Mform.ipw"), envir=environment())
clusterEvalQ(my.cluster, {
  library(VGAM)
  library(dplyr)
  library(tidyr)
  })
registerDoRNG(3308004)

#compute bootstrap estimates
ipwmed.boot <- foreach(i=1:nboot, .combine=cbind) %dopar% {

	boot.data <- plowData[sample(nrow(plowData), nrow(plowData), replace=TRUE),]

	boot.est <- ipwmed(data=boot.data, Dform=Dform.ipw, Lform=Lform.ipw, Mform=Mform.ipw)
	
	return(boot.est)
	}

stopCluster(my.cluster)
rm(my.cluster)

ipwmed.boot <- matrix(unlist(ipwmed.boot), ncol=4, byrow=TRUE)

#compute 95% CIs
ipwmed.output <- matrix(data=NA, nrow=4, ncol=3)

for (i in 1:nrow(ipwmed.output)) {

	ipwmed.output[i,1] <- round(ipwmed.est[i], digits=3)
	ipwmed.output[i,2] <- round(quantile(ipwmed.boot[,i], prob=0.025, na.rm=T), digits=3)
	ipwmed.output[i,3] <- round(quantile(ipwmed.boot[,i], prob=0.975, na.rm=T), digits=3)
	}

ipwmed.output <- data.frame(estimand = c("IDE(1,0)", "IIE(1,0)", "OE(1,0)", "CDE(1,0,1.8K)"), ipwmed.output)
colnames(ipwmed.output) <- c("estimand", "point.est", "ll.95ci", "ul.95ci")

##print output
sink(paste(logdir, "figure_4-9_log.txt", sep=""))

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

##display results in dot-whisket plot

rwr.output$method <- "RWR Estimates"
medsim.output$method <- "Simulation Estimates"
ipwmed.output$method <- "IPW Estimates"

combined.output <- rbind(rwr.output, medsim.output, ipwmed.output)

ggplot(data = combined.output, aes(y = estimand, x = point.est)) +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = ll.95ci, xmax = ul.95ci), height = 0.1) +
  facet_wrap(~ factor(method, levels = c("RWR Estimates", "Simulation Estimates", "IPW Estimates")), ncol = 3) +
  geom_vline(xintercept = 0, color = "black", size = 0.4) +
  scale_x_continuous(breaks = seq(-0.12, 0.08, by = 0.02)) +
  ylab("Estimand") +
  xlab("Difference in Proportions") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(panel.spacing = grid::unit(1, "lines"))

ggsave(paste(figdir, "figure_4-9.png", sep=""), height=4.5, width=9, units="in", dpi=600)


