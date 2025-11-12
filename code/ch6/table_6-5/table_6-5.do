/*Table 6.5*/
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
//we therefore implement DML using random forests only

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/Tatar/tatar.dta" ///
	"${datadir}Tatar\"

//open log
log using "${logdir}table_6-5.log", replace 

//load data
use "${datadir}Tatar\tatar.dta", clear

//note that -cmed mr- and -cmed dml- do not support using blocks of multiple 
//mediators in a causal sequence, so we focus on perceptions of threat only as 
//our focal mediators (i.e., fear_g1, fear_g2, fear_g3) 

//define macros for different variables
global C kulak prosoviet_pre religiosity_pre land_pre orchard_pre ///
	animals_pre carriage_pre otherprop_pre //baseline confounders
global D violence //exposure
global M1 fear_g1 //first mediator
global M2 fear_g2 //second mediator
global M3 fear_g3 //third mediator
global Y annex //outcome

//compute parametric MR (type mr2) estimates of path-specific effects
qui cmed mr $Y ($M1 $M2 $M3) $D = $C, paths censor(2 98) reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//compute DML (type mr2) estimates of path-specific effects
cmed dml $Y ($M1 $M2 $M3) $D = $C, paths method(rforest, iter(200) lsize(10)) censor(2 98) seed(60637)
	
log close

//note that the estimates differ slightly from those reported in
//the text, which are based on the R implementation.
