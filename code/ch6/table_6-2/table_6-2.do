/*Table 6.2*/
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
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/NLSY79/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

//open log
log using "${logdir}table_6-2.log", replace 

//load data
use "${datadir}nlsy79BK_ed2.dta", clear

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

//compute type mr1 estimates of natural effects

//parametric multiply robust estimates
qui cmed mr $Y $M $D = $C, rmpw censor(1 99) reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

//DML estimates
cmed dml $Y $M $D = $C, method(lasso) rmpw censor(1 99) seed(`myseed')

//compute type mr2 estimates of natural effects

//parametric multiply robust estimates
qui cmed mr $Y $M $D = $C, censor(1 99) reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

//DML estimates
cmed dml $Y $M $D = $C, method(lasso) censor(1 99) seed(`myseed')

log close

//note that the estimates differ slightly from those reported in
//the text, which are based on the R implementation. This is due to minor 
//differences in how the weights are censored and to the use of the LASSO
//exclusively rather than a super learner for the DML estimators

