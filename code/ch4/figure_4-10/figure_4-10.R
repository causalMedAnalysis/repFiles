###Figure 4.10###

rm(list=ls())

packages<-c("dplyr", "tidyr", "ggplot2", "gridExtra", "metR", "foreign")
	
#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

#install_github("xiangzhou09/rwrmed")

library(rwrmed)

##office
#datadir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/data/" 
#logdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/code/ch4/_LOGS/"
#figdir <- "C:/Users/Geoffrey Wodtke/Dropbox/shared/causal_mediation_text/figures/ch4/"

##home
datadir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/data/" 
logdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/code/ch4/_LOGS/"
figdir <- "C:/Users/Geoff/Dropbox/shared/causal_mediation_text/figures/ch4/"

sink(paste(logdir, "figure_4-10_log.txt", sep=""))

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

##compute RWR estimates
Lmodel.rwr <- lm(authGovCat~plow+agricultural_suitability+tropical_climate+large_animals+rugged, data=plowData)

Mform.rwr <- ln_income~plow+agricultural_suitability+tropical_climate+large_animals+rugged

Yform.rwr <- women_politics~(plow*ln_income)+authGovCat+agricultural_suitability+tropical_climate+large_animals+rugged

rwr.fit <- rwrmed(
	treatment="plow", 
	pre_cov=c("agricultural_suitability", "tropical_climate", "large_animals", "rugged"), 
	zmodels=list(Lmodel.rwr),
	m_form=Mform.rwr,
	y_form=Yform.rwr,
	data=plowData)

rwr.decomp <- decomp(rwr.fit, a0=0, a1=1, m=7.5, bootstrap=T, rep=200)

# Create the dataframe
rwr.output <- data.frame(
	point.est = c(rwr.decomp$twocomp[,1], rwr.decomp$fourcomp[1,1]),
	ll.95ci = c(rwr.decomp$twocomp[,3], rwr.decomp$fourcomp[1,3]),
	ul.95ci = c(rwr.decomp$twocomp[,4], rwr.decomp$fourcomp[1,4])
	)

rwr.output <- data.frame(estimand = c("IDE(1,0)", "IIE(1,0)", "OE(1,0)", "CDE(1,0,1.8K)"), rwr.output)

row.names(rwr.output) <- NULL

##specify range for sensitivity parameters under D-Y confounding
sens.grid <- expand.grid(delta_UYgivCDLM=seq(-0.1, 0.1, 0.01), delta_DUgivC=seq(-0.1, 0.1, 0.01))

##compute bias-adjusted estimates
adj.grid <- cbind(sens.grid,
	ide.adj=rwr.output[1,2]-(sens.grid$delta_UYgivCDLM*sens.grid$delta_DUgivC),
	oe.adj=rwr.output[3,2]-(sens.grid$delta_UYgivCDLM*sens.grid$delta_DUgivC))

##create contour plots of bias-adjusted estimates against sensitivity parameters
ide.plot <- ggplot(adj.grid, aes(x=delta_UYgivCDLM, y=delta_DUgivC, z=ide.adj, colour=stat(level))) +  
		geom_contour(breaks=seq(round(min(adj.grid$ide.adj), 3), round(max(adj.grid$ide.adj), 3), 0.002), show.legend=FALSE) +
		scale_colour_distiller(palette="Greys", direction=1) +
		xlab(expression(delta["UY|C,D,L,M"])) +
		ylab(expression(delta["DU|C"])) +
		scale_x_continuous(breaks=seq(-0.1, 0.1, 0.02)) +
		scale_y_continuous(breaks=seq(-0.1, 0.1, 0.02)) +
		ggtitle("A. Bias-adjusted IDE(1,0) Estimates") +
		theme_bw(base_size=11) +
		theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
		geom_text_contour(
			breaks=seq(round(min(adj.grid$ide.adj), 3), round(max(adj.grid$ide.adj), 3), 0.002),
			stroke=0.3,
			size=3,
			skip=0,
			color="black")

oe.plot <- ggplot(adj.grid, aes(x=delta_UYgivCDLM,y=delta_DUgivC,z=oe.adj,colour=stat(level))) +  
		geom_contour(breaks=seq(round(min(adj.grid$oe.adj),3),round(max(adj.grid$oe.adj),3),0.002),show.legend=FALSE) +
		scale_colour_distiller(palette="Greys",direction=1) +
		xlab(expression(delta["UY|C,D,L,M"])) +
		ylab(expression(delta["DU|C"])) +
		scale_x_continuous(breaks=seq(-0.1,0.1,0.02)) +
		scale_y_continuous(breaks=seq(-0.1,0.1,0.02)) +
		ggtitle("B. Bias-adjusted OE(1,0) Estimates") +
		theme_bw(base_size=11) +
		theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
		geom_text_contour(
			breaks=seq(round(min(adj.grid$oe.adj),3),round(max(adj.grid$oe.adj),3),0.002),
			stroke=0.3,
			size=3,
			skip=0,
			color="black")

comb.plot <- grid.arrange(ide.plot, oe.plot, ncol=1)

ggsave(paste(figdir, "figure_4-10.png", sep=""), plot=comb.plot, height=8, width=4.5, units="in", dpi=600)

sink()
