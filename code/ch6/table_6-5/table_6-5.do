/*Table 6.5*/
capture clear all
capture log close
set more off

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis
net install parallel, from("https://raw.github.com/gvegayon/parallel/stable/") replace //module for parallelization
mata mata mlib index

ssc install rforest, replace 
ssc install lassopack, replace

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/main/data/Tatar/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

//open log
log using "${logdir}table_6-5.log", replace 

//load data
use "${datadir}tatar.dta", clear

//note that -cmed mr- and -cmed dml- do not support using blocks of multiple 
//mediators in a causal sequence, so we focus on perceptions of threat only as 
//our focal mediators (i.e., fear_g1, fear_g2, fear_g3); see the R package 
//cmedR for support using blocks of multiple mediators

//define macros for different variables
global C kulak prosoviet_pre religiosity_pre land_pre orchard_pre ///
	animals_pre carriage_pre otherprop_pre //baseline confounders
global D violence //exposure
global M1 fear_g1 //first mediator
global M2 fear_g2 //second mediator
global M3 fear_g3 //third mediator
global Y annex //outcome

//set seed 
//note: opt parallel requires a separate seed for each child process
local myseed 3308004
qui parallel numprocessors
local ncores = max(floor(r(numprocessors) * 0.75), 1)
local parseeds 
forval i = 1/`ncores' {
    local newseed = `myseed' + `i'
    local parseeds `parseeds' " `newseed'"
}

//note that stata does not currently support use of a superLearner
//we therefore implement DML using the LASSO only

//compute parametric MR (type mr2) estimates of path-specific effects
qui cmed mr $Y ($M1 $M2 $M3) $D = $C, paths censor(1 99) ///
	reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

//compute DML (type mr2) estimates of path-specific effects
cmed dml $Y ($M1 $M2 $M3) $D = $C, paths ///
	method(lasso) censor(1 99) seed(`myseed')
	
log close

//note that the estimates differ slightly from those reported in
//the text, which are based on the R implementation. This is due to minor 
//differences in how the weights are censored and to the use of the LASSO
//exclusively rather than a super learner for the DML estimator
