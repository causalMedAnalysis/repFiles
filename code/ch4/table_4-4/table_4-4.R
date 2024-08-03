###Table 4.4###

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

##input data 
nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("cesd_age40", "ever_unemp_age3539", "faminc_adj_age3539", "log_faminc_adj_age3539", 
	"att22", "female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")]),]

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##weighting estimators

#define estimation function 
ipwmed <- function(data, Dform, Lform, Mform) {

	df <- data

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

	OE <- Ehat_Y1M1-Ehat_Y0M0
	IDE <- Ehat_Y1M0-Ehat_Y0M0
	IIE <- Ehat_Y1M1-Ehat_Y1M0

	Ymodel.wtd <- lm(std_cesd_age40 ~ att22 * log_faminc_adj_age3539, data=df.wts, weights=sw4)

	CDE <- Ymodel.wtd$coefficients["att22"] + log(50000)*Ymodel.wtd$coefficients["att22:log_faminc_adj_age3539"]

	##### NEED TO ADD CDE ESTIMATES #####

	point.est <- list(OE, IDE, IIE, CDE)

	return(point.est)
	}

#compute ipw estimates w/ additive models
Dform.a <- att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

Lform.a <- ever_unemp_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

Mform.a <- log_faminc_adj_age3539 ~ att22 + ever_unemp_age3539 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

ipwmed.est.a <- ipwmed(data=nlsy, Dform=Dform.a, Lform=Lform.a, Mform=Mform.a)

#compute ipw estimates w/ interactive models
Dform.x <- att22 ~ female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

Lform.x <- ever_unemp_age3539 ~ att22 * (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

Mform.x <- log_faminc_adj_age3539 ~ att22 * (ever_unemp_age3539 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

ipwmed.est.x <- ipwmed(data=nlsy, Dform=Dform.x, Lform=Lform.x, Mform=Mform.x)

#print output
ipwmed.est.a <- matrix(unlist(ipwmed.est.a), ncol=4, byrow=TRUE)
ipwmed.est.x <- matrix(unlist(ipwmed.est.x), ncol=4, byrow=TRUE)

ipwmed.a.out <- data.frame(est = round(c(ipwmed.est.a), digits=3), row.names=c("OEhat^ipw", "IDEhat^ipw", "IIEhat^ipw", "CDEhat^ipw"))
ipwmed.x.out <- data.frame(est = round(c(ipwmed.est.x), digits=3), row.names=c("OEhat^ipw", "IDEhat^ipw", "IIEhat^ipw", "CDEhat^ipw"))

sink(paste(logdir, "table_4-4_log.txt", sep=""))

cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("ipwEst w/ additive models for D, L, and ln(M)\n")
print(ipwmed.a.out)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("ipwEst w/ interactive models for L and ln(M)\n")
print(ipwmed.x.out)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

sink()
