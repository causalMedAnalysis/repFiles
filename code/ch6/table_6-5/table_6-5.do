/*Table 6.5*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/mrpath, replace 
github install causalMedAnalysis/dmlpath, replace 

//install dependencies for dmlmed
ssc install rforest, replace 
ssc install lassopack, replace

//note that stata does not currently support use of a superLearner
//we therefore implement DML using random forests only

//specify directories 
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch6\_LOGS\"

//download data
copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/Tatar/tatar.dta" ///
	"${datadir}Tatar\"

//open log
log using "${logdir}table_6-5.log", replace 

//load data
use "${datadir}Tatar\tatar.dta", clear

//note also that the stata modules mrpath and dmlpath do not support using blocks
//of multiple mediators in a causal sequence, so we focus on perceptions of 
//threat only (fear_g1, fear_g2, fear_g3) as our focal mediators

//define macros for different variables
global C kulak prosoviet_pre religiosity_pre land_pre orchard_pre ///
	animals_pre carriage_pre otherprop_pre //baseline confounders
global D violence //exposure
global M1 fear_g1 //first mediator
global M2 fear_g2 //second mediator
global M3 fear_g3 //third mediator
global Y annex //outcome

//compute parametric MR (type mr2) estimates of path-specific effects
qui mrpath $Y $M1 $M2 $M3, dvar($D) d(1) dstar(0) cvars($C) censor(2 98) ///
	seed(60637) reps(2000) 

mat list e(b)
mat list e(ci_percentile)

//compute DML (type mr2) estimates of path-specific effects
dmlpath $Y $M1 $M2 $M3, model(rforest) dvar($D) d(1) dstar(0) cvars($C) ///
	censor(2 98) seed(60637) iter(200) lsize(10)

log close

//note that the estimates differ slightly from those reported in
//the text, which are based on the R implementation
