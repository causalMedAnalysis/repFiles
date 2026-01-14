/*Table 5.5*/
capture clear all
capture log close
set more off

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis
net install parallel, from("https://raw.github.com/gvegayon/parallel/stable/") replace //module for parallelization
mata mata mlib index

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/NLSY79/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch5\_LOGS\"

//open log
log using "${logdir}table_5-5.log", replace 

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
global M1 ever_unemp_age3539 //first mediator
global M2 log_faminc_adj_age3539 //second mediator
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

//compute PSE estimates from additive linear models
qui cmed linear $Y ($M2 $M1) $D = $C, paths nointer reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

//compute PSE estimates from linear model estimates w/ Dx(M1,M2) interactions
qui cmed linear $Y ($M2 $M1) $D = $C, paths reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

//compute PSE estimates using inverse probability weighting
qui cmed ipw $Y ($M2 $M1) $D = $C, paths censor(1 99) reps(2000) parallel seed(`parseeds')

mat list e(b)
mat list e(ci_percentile)

log close

//some of the estimates differ slightly from those reported in the text, 
//which are based on the R implementation. This is due to differences in
//random number seeding for the bootstrap and minor differences in how the 
//weights are censored by -cmed ipw-
