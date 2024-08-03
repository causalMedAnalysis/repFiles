###Table 4.3###

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

set.seed(3308004)

##input data 
nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("cesd_age40", "ever_unemp_age3539", "faminc_adj_age3539", "log_faminc_adj_age3539", 
	"att22", "female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")]),]

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##simulation estimators

#define estimation function 
medsim <- function(data, num.sim=100, Lform, Mform, Yform) {
	
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

	OE <- mean(Y1L1M1) - mean(Y0L0M0)
	IDE <- mean(Y1L1M0) - mean(Y0L0M0)  
	IIE <- mean(Y1L1M1) - mean(Y1L1M0)
	CDE <- mean(Y1L1m) - mean(Y0L0m)

	point.est <- list(OE, IDE, IIE, CDE)

	return(point.est)
	}

#compute simulation estimates w/ exposure-mediator interaction
Lform.x <- ever_unemp_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

Mform.x <- log_faminc_adj_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

Yform.x <- std_cesd_age40~(log_faminc_adj_age3539 * att22) + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + ever_unemp_age3539

medsim.est.x <- medsim(data=nlsy, num.sim=2000, Lform=Lform.x, Mform=Mform.x, Yform=Yform.x)

#compute simulation estimates w/ all two-way interactions
Lform.xx <- ever_unemp_age3539 ~ att22 * (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

Mform.xx <- log_faminc_adj_age3539 ~ att22 * (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

Yform.xx <- std_cesd_age40 ~ (log_faminc_adj_age3539 * att22) + 
	(female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3) * 
	(log_faminc_adj_age3539 + att22) + (ever_unemp_age3539 * log_faminc_adj_age3539)

medsim.est.xx <- medsim(data=nlsy, num.sim=2000, Lform=Lform.xx, Mform=Mform.xx, Yform=Yform.xx)

#print output
medsim.est.x <- matrix(unlist(medsim.est.x), ncol=4, byrow=TRUE)
medsim.est.xx <- matrix(unlist(medsim.est.xx), ncol=4, byrow=TRUE)

medsim.x.out <- data.frame(est = round(c(medsim.est.x), digits=3), row.names=c("OEhat^sim", "IDEhat^sim", "IIEhat^sim", "CDEhat^sim"))
medsim.xx.out <- data.frame(est = round(c(medsim.est.xx), digits=3), row.names=c("OEhat^sim", "IDEhat^sim", "IIEhat^sim", "CDEhat^sim"))

sink(paste(logdir, "table_4-3_log.txt", sep=""))

cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("simEst w/ ln(M) and ln(M)xD Interaction\n")
print(medsim.x.out)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("simEst w/ ln(M) and All Two-way Interactions\n")
print(medsim.xx.out)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

sink()
