/*Table 3.5*/
capture clear all
capture log close
set more off

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/NLSY79/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch3\_LOGS\"

//open log
log using "${logdir}table_3-5.log", replace 

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

//set seed 
set seed 3308004

//compute interval estimates using bootstrap
qui cmed ipw $Y $M $D = $C, censor(1 99) reps(2000) ///
	saving("bootmed.dta", replace)

mat list e(ci_percentile)

qui cmed ipw $Y ((logit) $M) $D = $C, m(0) censor(1 99) reps(2000) ///
	saving("bootcde.dta", replace)

mat list e(ci_percentile)

//compute p-values for null of no effect
use "bootmed.dta", clear

foreach x in ATE NDE NIE {
	qui gen `x'_lt=0 
	qui replace `x'_lt=1 if _b_`x'<0
	qui sum `x'_lt
	scalar p_lt = r(mean)
	
	qui gen `x'_rt=0 
	qui replace `x'_rt=1 if _b_`x'>0
	qui sum `x'_rt
	scalar p_rt = r(mean)
	
	di "`x' p-value"
	di min(p_rt, p_lt)*2
	di "-----------"
	
	drop `x'_*
	}

use "bootcde.dta", clear

qui gen CDE_lt=0 
qui replace CDE_lt=1 if _b_CDE<0
qui sum CDE_lt
scalar p_lt = r(mean)

qui gen CDE_rt=0 
qui replace CDE_rt=1 if _b_CDE>0
qui sum CDE_rt
scalar p_rt = r(mean)

di "CDE p-value"
di min(p_rt, p_lt)*2

//erase bootstrap data	
erase "bootmed.dta"
erase "bootcde.dta"

log close

//note that the cmed interval estimates and pvals differ slightly 
//from those reported in the text, which are based on the R implementation. 
//This is due to minor differences in how the weights are censored and also to
//differences in random numder seeding that influence the bootstrap samples 
