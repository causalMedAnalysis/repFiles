/*Table 3.3*/
capture clear all
capture log close
set more off
set maxvar 10000

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/cmed //module to perform causal mediation analysis

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch3\_LOGS\"

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_3-3.log", replace 

//load data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, female, black, ///
	hispan, paredu, parprof, parinc_prank, famsize, afqt3)
	
//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//define macros for different variables
global C female black hispan paredu parprof parinc_prank famsize afqt3 //baseline confounders
global D att22 //exposure
global M ever_unemp_age3539 //mediator
global Y std_cesd_age40 //outcome

//set seed to ensure reproducibility
set seed 3308004

//compute point estimates using simulation w/ additive models
qui cmed sim ((regress) $Y) ((logit) $M) $D = $C, nsim(2000) nointer reps(2)
mat list e(b)

qui cmed impute ((regress) $Y) $M $D = $C, m(0) nointer //use imputation for cde
mat list e(b)

//compute point estimates using simulation w/ with DxM, CxD, and CxM interactions

qui cmed sim ((regress) $Y) ((logit) $M) $D = $C, nsim(2000) cxd cxm reps(2)
mat list e(b)

qui cmed impute ((regress) $Y) $M $D = $C, m(0) cxd cxm //use imputation for cde
mat list e(b)

log close

//note that the cmed estimates differ slightly from those reported in the text,
//which are based on the R implementation. This is due only to monte carlo
//error and differences in random number seeding.
