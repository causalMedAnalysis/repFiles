/*Table 4.4*/
capture clear all
capture log close
set more off

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/NLSY79/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch4\_LOGS\"

//open log
log using "${logdir}table_4-4.log", replace 

//load data
use "${datadir}nlsy79BK_ed2.dta"

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, log_faminc_adj_age3539, ///
	female, black, hispan, paredu, parprof, parinc_prank, famsize, afqt3)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//define macros for different variables
global C female black hispan paredu parprof parinc_prank famsize afqt3 //baseline confounders
global D att22 //exposure
global L ever_unemp_age3539 //exposure-induced confounder
global M log_faminc_adj_age3539 //mediator
global Y std_cesd_age40 //outcome

//compute point estimates using ipw estimates w/ additive models
qui cmed ipw $Y ((regress) $M) ((logit) $L) $D = $C, censor(1 99)
mat list e(b)

qui cmed ipw $Y ((regress) $M) ($L) $D = $C, m(10.82) censor(1 99)
mat list e(b)

//compute point estimates using ipw w/ CxD and LxD interactions
qui cmed ipw $Y ((regress) $M) ((logit) $L) $D = $C, censor(1 99) cxd lxd
mat list e(b)

qui cmed ipw $Y ((regress) $M) ($L) $D = $C, m(10.82) censor(1 99) cxd lxd
mat list e(b)

log close 

//note that the cmed estimates differ slightly from those reported in
//the text, which are based on the R implementation. This is due to 
//minor differences in how the weights are censored.
