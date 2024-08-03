###Figure 3.10###

rm(list=ls())

packages<-c("dplyr", "tidyr", "ggplot2", "gridExtra", "metR", "foreign")
	
#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

##office
#datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"
#figdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/figures/ch3/"

##home
datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch3/_LOGS/"
figdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/figures/ch3/"

sink(paste(logdir, "figure_3-10_log.txt", sep=""))

##input data
jobs <- read.dta(paste(datadir, "JOBSII/Jobs-NoMiss-Binary.dta", sep=""))

jobs <- jobs %>% 
	mutate(
		treat = recode(treat, "control" = 0, "exp" = 1),
		work1 = recode(work1, "psyump" = 0, "psyemp" = 1),
		nonwhite = recode(nonwhite, "white0" = 0, "non.white1" = 1),
		educ = recode(educ, "lt-hs" = 1, "highsc" = 2, "somcol" = 3, "bach" = 4, "gradwk" = 5),
		income = recode(income, "lt15k" = 1, "15t24k" = 2, "25t39k" = 3, "40t49k" = 4, "50k+" = 5))


##compute point estimates using linear models w/ exposure-mediator interaction

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

##specify range for sensitivity parameters under M-Y confounding
sens.grid <- expand.grid(delta_UYgivCDM=seq(-0.1, 0.1, 0.01), delta_DUgivCM=seq(-0.1, 0.1, 0.01))

##compute bias-adjusted estimates
adj.grid <- cbind(sens.grid,
	nde.adj=linmedx.est[1,2]-(sens.grid$delta_UYgivCDM*sens.grid$delta_DUgivCM),
	nie.adj=linmedx.est[1,3]+(sens.grid$delta_UYgivCDM*sens.grid$delta_DUgivCM))

##create contour plots of bias-adjusted estimates against sensitivity parameters
nde.plot <- ggplot(adj.grid, aes(x=delta_UYgivCDM, y=delta_DUgivCM, z=nde.adj, colour=stat(level))) +  
		geom_contour(breaks=seq(round(min(adj.grid$nde.adj), 3), round(max(adj.grid$nde.adj), 3), 0.002), show.legend=FALSE) +
		scale_colour_distiller(palette="Greys", direction=1) +
		xlab(expression(delta["UY|C,D,M"])) +
		ylab(expression(delta["DU|C,M"])) +
		scale_x_continuous(breaks=seq(-0.1, 0.1, 0.02)) +
		scale_y_continuous(breaks=seq(-0.1, 0.1, 0.02)) +
		ggtitle("A. Bias-adjusted NDE(1,0) Estimates") +
		theme_bw(base_size=11) +
		theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
		geom_text_contour(
			breaks=seq(round(min(adj.grid$nde.adj), 3), round(max(adj.grid$nde.adj), 3), 0.002),
			stroke=0.3,
			size=3,
			skip=0,
			color="black")

nie.plot <- ggplot(adj.grid, aes(x=delta_UYgivCDM,y=delta_DUgivCM,z=nie.adj,colour=stat(level))) +  
		geom_contour(breaks=seq(round(min(adj.grid$nie.adj),3),round(max(adj.grid$nie.adj),3),0.002),show.legend=FALSE) +
		scale_colour_distiller(palette="Greys",direction=1) +
		xlab(expression(delta["UY|C,D,M"])) +
		ylab(expression(delta["DU|C,M"])) +
		scale_x_continuous(breaks=seq(-0.1,0.1,0.02)) +
		scale_y_continuous(breaks=seq(-0.1,0.1,0.02)) +
		ggtitle("B. Bias-adjusted NIE(1,0) Estimates") +
		theme_bw(base_size=11) +
		theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
		geom_text_contour(
			breaks=seq(round(min(adj.grid$nie.adj),3),round(max(adj.grid$nie.adj),3),0.002),
			stroke=0.3,
			size=3,
			skip=0,
			color="black")

comb.plot <- grid.arrange(nde.plot, nie.plot, ncol=1)

ggsave(paste(figdir, "figure_3-10.png", sep=""), plot=comb.plot, height=8, width=4.5, units="in", dpi=600)

sink()
