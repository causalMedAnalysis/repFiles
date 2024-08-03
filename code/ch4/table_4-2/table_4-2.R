###Table 4.2###

rm(list=ls())

packages<-c("dplyr", "tidyr", "devtools")

#install.packages(packages)

for (package.i in packages) {
	suppressPackageStartupMessages(library(package.i, character.only=TRUE))
	}

#install_github("xiangzhou09/rwrmed")
library(rwrmed)

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

##RWR estimators

#OEhat^rwr, IDEhat^rwr, IIEhat^rwr, CDEhat^rwr
Lmodel <- lm(ever_unemp_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3, data=nlsy)

Mform <- faminc_adj_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

Yform <- std_cesd_age40 ~ (faminc_adj_age3539 * att22) + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + ever_unemp_age3539

rwr <- rwrmed(
	treatment="att22", 
	pre_cov=c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"), 
	zmodels=list(Lmodel),
	m_form=Mform,
	y_form=Yform,
	data=nlsy)

rwr.decomp <- decomp(rwr, a0=0, a1=1, m=50000, bootstrap=F)

rwr.est <- data.frame(est = c(rwr.decomp$twocomp[,1], rwr.decomp$fourcomp[1,1]))

rownames(rwr.est)<-c("IDE", "IIE", "OE", "CDE")

#OEhat^rwr, IDEhat^rwr, IIEhat^rwr, CDEhat^rwr w/ log(income)

Lmodel <- lm(ever_unemp_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3, data=nlsy)

Mform <- log_faminc_adj_age3539 ~ att22 + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3

Yform <- std_cesd_age40 ~ (log_faminc_adj_age3539 * att22) + female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3 + ever_unemp_age3539

rwr.ln <- rwrmed(
	treatment="att22", 
	pre_cov=c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"), 
	zmodels=list(Lmodel),
	m_form=Mform,
	y_form=Yform,
	data=nlsy)

rwr.ln.decomp <- decomp(rwr.ln, a0=0, a1=1, m=10.82, bootstrap=F)

rwr.ln.est <- data.frame(est = c(rwr.ln.decomp$twocomp[,1], rwr.ln.decomp$fourcomp[1,1]))

rownames(rwr.ln.est)<-c("IDE", "IIE", "OE", "CDE")

#OEhat^rwr+, IDEhat^rwr+, IIEhat^rwr+, CDEhat^rwr+ w/ log(income)

Lmodel <- lm(ever_unemp_age3539 ~ att22 * (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3), data=nlsy)

Mform <- log_faminc_adj_age3539 ~ att22 * (female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3)

Yform <- std_cesd_age40 ~ (log_faminc_adj_age3539 * att22) + 
	(female + black + hispan + paredu + parprof + parinc_prank + famsize + afqt3) * 
	(log_faminc_adj_age3539 + att22) + (ever_unemp_age3539 * log_faminc_adj_age3539)

rwr.ln.x <- rwrmed(
	treatment="att22", 
	pre_cov=c("female", "black", "hispan", "paredu", "parprof", "parinc_prank", "famsize", "afqt3"), 
	zmodels=list(Lmodel),
	m_form=Mform,
	y_form=Yform,
	data=nlsy)

rwr.ln.x.decomp <- decomp(rwr.ln.x, a0=0, a1=1, m=10.82, bootstrap=F)

rwr.ln.x.est <- data.frame(est = c(rwr.ln.x.decomp$twocomp[,1], rwr.ln.x.decomp$fourcomp[1,1]))

rownames(rwr.ln.x.est)<-c("IDE", "IIE", "OE", "CDE")

##print output
sink(paste(logdir, "table_4-2_log.txt", sep=""))

cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("RWR w/ DxM Interaction\n")
print(rwr.est)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("RWR w/ ln(M) and Dxln(M) Interaction\n")
print(rwr.ln.est)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("RWR w/ ln(M) and All Two-way Interactions\n")
print(rwr.ln.x.est)
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

sink()
