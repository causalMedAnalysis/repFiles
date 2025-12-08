/*Table 5.8*/
capture clear all
capture log close
set more off

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis
net install parallel, from("https://raw.github.com/gvegayon/parallel/stable/") replace //module for parallelization
mata mata mlib index

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/Brader_et_al2008/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch5\_LOGS\"

//open log
log using "${logdir}table_5-8.log", replace 

//load data
use "${datadir}Brader_et_al2008.dta", clear

//keep complete cases
drop if missing(immigr, emo, p_harm, tone_eth, ppage, ppeducat, ppgender, ppincimp)

//standardize outcome
egen std_immigr=std(4-immigr)

//dummy code controls
qui tab ppeducat, gen(ppeducat_)
drop ppeducat_1
rename ppeducat_2 hs
rename ppeducat_3 sc
rename ppeducat_4 ba

qui tab ppgender, gen(ppgender_)
drop ppgender_1
rename ppgender_2 female

//define macros for different variables
global C ppage female hs sc ba ppincimp //baseline confounders
global D tone_eth //exposure
global M1 p_harm //first mediator
global M2 emo //second mediator
global Y std_immigr //outcome

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

//compute PSE estimates from additive linear models
qui cmed linear $Y ($M1 $M2) $D = $C, paths nointer ///
	reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

//compute PSE estimates from linear model estimates w/ Dx(M1,M2) interactions
qui cmed linear $Y ($M1 $M2) $D = $C, paths ///
	reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

//compute PSE estimates using inverse probability weighting
qui cmed ipw $Y ($M1 $M2) $D = $C, paths censor(1 99) ///
	reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

log close

//some of the estimates differ slightly from those reported in the text, 
//which are based on the R implementation. This is due to differences in
//random number seeding for the bootstrap and minor differences in how the 
//weights are censored by -cmed ipw-


