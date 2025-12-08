/*Table 3.7*/
capture clear all
capture log close
set more off
set maxvar 10000

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis
net install parallel, from("https://raw.github.com/gvegayon/parallel/stable/") replace //module for parallelization
mata mata mlib index

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/JOBSII/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch3\_LOGS\"

//open log
log using "${logdir}table_3-7.log", replace 

//load data
use "${datadir}Jobs-NoMiss-Binary.dta", clear

//define macros for different variables
global C econ_hard sex age nonwhite educ income //baseline confounders
global D treat //exposure
global M job_seek //mediator
global Y work1 //outcome

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

//compute point and interval estimates based on linear models
cmed linear $Y $M $D = $C, reps(2000) parallel seed(`parseeds')
cmed linear $Y $M $D = $C, m(4) reps(2000) parallel seed(`parseeds')

//compute point and interval estimates based on simulation/imputation
cmed sim ((logit) $Y) ((regress) $M) $D = $C, nsim(1000) reps(2000) parallel seed(`parseeds')
cmed impute ((logit) $Y) $M $D = $C, m(4) reps(2000) parallel seed(`parseeds')

//compute point and interval estimates based on inverse probability weighting
cmed ipw $Y $M $D = $C, censor(1 99) reps(2000) parallel seed(`parseeds')
cmed ipw $Y ((regress) $M) $D = $C, m(4) censor(1 99) reps(2000) parallel seed(`parseeds')

log close

//some of the estimates differ slightly from those reported in the text, 
//which are based on the R implementation. This is variously due to minor 
//differences in how the weights are censored, to Monte Carlo error, or to 
//differences in random number seeding that influence the bootstrap samples

//note: -cmed sim- with large nsim() and reps() can be time consuming to run;
//parallelization of the bootstrap is highly recommended in these instances
