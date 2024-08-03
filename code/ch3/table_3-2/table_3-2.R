###Table 3.2###

rm(list=ls())

packages<-c("dplyr", "tidyr")

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

sink(paste(logdir, "table_3-2_log.txt", sep=""))

##input data
nlsy <- as.data.frame(readRDS(paste(datadir, "NLSY79/nlsy79BK_ed2.RDS", sep="")))

nlsy <- nlsy[complete.cases(nlsy[,c("cesd_age40", "ever_unemp_age3539", "att22", "female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3")]),]

nlsy$std_cesd_age40 <- (nlsy$cesd_age40-mean(nlsy$cesd_age40))/sd(nlsy$cesd_age40)

##linear model estimators

#linear model w/o exposure-mediator interaction

linmed <- function(data) {
	df <- data

	Mmodel <- lm(ever_unemp_age3539~att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df)

	Ymodel <- lm(std_cesd_age40~ever_unemp_age3539+att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df)

	NDE <- CDE <- Ymodel$coefficients["att22"]
	NIE <- Mmodel$coefficients["att22"]*Ymodel$coefficients["ever_unemp_age3539"]
	ATE <- NDE+NIE
	
	point.est <- list(ATE, NDE, NIE, CDE)

	return(point.est)
	}

linmed.est <- linmed(nlsy)
linmed.est <- matrix(unlist(linmed.est), ncol=4, byrow=TRUE)

#linear model w/ exposure-mediator interaction

linmedx <- function(data) {
	df <- data

	df <- df %>% 
			mutate(
				female = female-mean(female), 
				black = black-mean(black), 
				hispan = hispan-mean(hispan), 
				paredu = paredu-mean(paredu), 
				parprof = parprof-mean(parprof), 
				parinc_prank = parinc_prank-mean(parinc_prank), 
				famsize = famsize-mean(famsize), 
				afqt3 = afqt3-mean(afqt3))

	Mmodel <- lm(ever_unemp_age3539~att22+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df)

	Ymodel <- lm(std_cesd_age40~(ever_unemp_age3539*att22)+female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3, data=df)

	NDE <- (Ymodel$coefficients["att22"] + Ymodel$coefficients["ever_unemp_age3539:att22"]*Mmodel$coefficients["(Intercept)"])
	NIE <- Mmodel$coefficients["att22"]*(Ymodel$coefficients["ever_unemp_age3539"] + Ymodel$coefficients["ever_unemp_age3539:att22"])
	ATE <- NDE+NIE
	CDE0 <- Ymodel$coefficients["att22"]
	
	point.est <- list(ATE, NDE, NIE, CDE0)

	return(point.est)
	}

linmedx.est <- linmedx(nlsy)
linmedx.est <- matrix(unlist(linmedx.est), ncol=4, byrow=TRUE)

#linear model w/ exposure-mediator interaction and covariate-exposure-mediator interactions

linmedxx <- function(data) {
	df <- data

	df <- df %>% 
			mutate(
				female = female-mean(female), 
				black = black-mean(black), 
				hispan = hispan-mean(hispan), 
				paredu = paredu-mean(paredu), 
				parprof = parprof-mean(parprof), 
				parinc_prank = parinc_prank-mean(parinc_prank), 
				famsize = famsize-mean(famsize), 
				afqt3 = afqt3-mean(afqt3))

	Mmodel <- lm(ever_unemp_age3539~att22*(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3), data=df)

	Ymodel <- lm(std_cesd_age40~(ever_unemp_age3539*att22)+(female+black+hispan+paredu+parprof+parinc_prank+famsize+afqt3)*(ever_unemp_age3539+att22), data=df)

	NDE <- (Ymodel$coefficients["att22"] + Ymodel$coefficients["ever_unemp_age3539:att22"]*Mmodel$coefficients["(Intercept)"])
	NIE <- Mmodel$coefficients["att22"]*(Ymodel$coefficients["ever_unemp_age3539"] + Ymodel$coefficients["ever_unemp_age3539:att22"])
	ATE <- NDE+NIE
	CDE0 <- Ymodel$coefficients["att22"]
	
	point.est <- list(ATE, NDE, NIE, CDE0)

	return(point.est)
	}

linmedxx.est <- linmedxx(nlsy)
linmedxx.est <- matrix(unlist(linmedxx.est), ncol=4, byrow=TRUE)

##print output

output <- t(rbind(linmed.est, linmedx.est, linmedxx.est))

output <- data.frame(output, row.names=c("ATEhat", "NDEhat", "NIEhat", "CDE0hat"))
colnames(output) <- c("^lma", "^lmi", "^lmi+")

print(output)

sink()
