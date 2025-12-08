/*Table 3.2*/
capture clear all
capture log close
set more off

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/NLSY79/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch3\_LOGS\"

//open log
log using "${logdir}table_3-2.log", replace 

//load data
use "${datadir}nlsy79BK_ed2.dta", clear

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

//compute point estimates using additive linear models
qui cmed linear $Y $M $D = $C, nointer
mat list e(b)

qui cmed linear $Y $M $D = $C, nointer m(0)
mat list e(b)

//compute point estimates from linear models with DxM interaction
qui cmed linear $Y $M $D = $C
mat list e(b)

qui cmed linear $Y $M $D = $C, m(0)
mat list e(b)
		
//compute point estimates from linear models with DxM, CxD, and CxM interactions
qui cmed linear $Y $M $D = $C, cxd cxm
mat list e(b)

qui cmed linear $Y $M $D = $C, m(0) cxd cxm
mat list e(b)

log close
