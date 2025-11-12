/*Table 6.2*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/cmed //module to perform causal mediation analysis

//install dependencies
ssc install rforest, replace 
ssc install lassopack, replace

//note that stata does not currently support use of a superLearner
//we therefore implement DML using the LASSO only

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_6-2.log", replace 

//load data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, log_faminc_adj_age3539, ///
	female, black, hispan, paredu, parprof, parinc_prank, famsize, afqt3)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//define macros for different variables
global C female black hispan paredu parprof parinc_prank famsize afqt3 //baseline confounders
global D att22 //exposure
global M ever_unemp_age3539 //mediator
global Y std_cesd_age40 //outcome

//compute type mr1 estimates of natural effects

//parametric multiply robust estimates
qui cmed mr $Y $M $D = $C, rmpw censor(1 99) reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//DML estimates
cmed dml $Y $M $D = $C, method(lasso) rmpw censor(1 99) seed(60637)

//compute type mr2 estimates of natural effects

//parametric multiply robust estimates
qui cmed mr $Y $M $D = $C, censor(1 99) reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//DML estimates
cmed dml $Y $M $D = $C, method(lasso) censor(1 99) seed(60637)

log close

//note that the estimates differ slightly from those reported in
//the text, which are based on the R implementation. This is due to minor 
//differences in how the weights are censored and to the use of the LASSO
//exclusively rather than a super learner for the DML estimators

