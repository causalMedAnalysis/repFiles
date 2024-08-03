/*Table 3.5*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"
*global figdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\figures\ch3\" 

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"
global figdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\figures\ch3\" 

log using "${logdir}figure_3-7.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

///inverse probability weighting (ipw) estimators

//ATEhat^ipw, NDEhat^ipw, NIEhat^ipw

//define function for ipw estimator for ATE, NDE, and NIE
capture program drop ipwmed
program define ipwmed, rclass

	//fit logit model for D given C
	logit att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	est store DModel_1

	//fit logit model for D given C and M
	logit att22 ///
		ever_unemp_age3539 female black hispan paredu parprof parinc_prank famsize afqt3
	est store DModel_2

	//predict exposure probabilities
	est restore DModel_1
	predict phat_d_C, pr
	gen phat_dstar_C=1-phat_d_C

	est restore DModel_2
	predict phat_d_CM, pr
	gen phat_dstar_CM=1-phat_d_CM

	//compute stabilized weights
	logit att22
	predict phat_d, pr
	gen phat_dstar=1-phat_d
	
	gen _sw1=phat_dstar/phat_dstar_C if att22==0
	gen _sw2=phat_d/phat_d_C if att22==1
	gen _sw3=(phat_dstar_CM*phat_d)/(phat_d_CM*phat_dstar_C) if att22==1

	//censor stabilized weights at 1st and 99th percentiles
	foreach i of var _sw* {
		centile `i', c(1 99) 
		replace `i'=r(c_1) if `i'<r(c_1) & `i'!=.
		replace `i'=r(c_2) if `i'>r(c_2) & `i'!=.
		}

	//compute weighted means of Y
	reg std_cesd_age40 [pw=_sw1] if att22==0
	local Ehat_Y0M0=_b[_cons]

	reg std_cesd_age40 [pw=_sw2] if att22==1
	local Ehat_Y1M1=_b[_cons]

	reg std_cesd_age40 [pw=_sw3] if att22==1
	local Ehat_Y1M0=_b[_cons]

	//compute effect estimates
	return scalar ATE=`Ehat_Y1M1'-`Ehat_Y0M0'
	return scalar NDE=`Ehat_Y1M0'-`Ehat_Y0M0'
	return scalar NIE=`Ehat_Y1M1'-`Ehat_Y1M0'
	
	drop phat_d_C phat_dstar_C phat_d_CM phat_dstar_CM phat_d phat_dstar _sw* ///
		_est_DModel_1 _est_DModel_2
	
	//reset estimation sample for bootstrap
	reg att22
	
end

//plot histogram
quietly bootstrap ///
	ATE=r(ATE) NDE=r(NDE) NIE=r(NIE), ///
	saving("${datadir}\bootmed.dta", replace) ///
	reps(2000) seed(3308004): ipwmed

use "${datadir}\bootmed.dta", clear

centile NIE, c(2.5 97.5)

set scheme s2mono

hist NIE, ///
	width(0.00125) ///
	yscale(range(0 180)) ///
	xscale(range(-0.05 0.02)) ///
	ylabel(0(20)180) ///
	xlabel(-0.05(0.01)0.02) ///
	xline(-.0295328 .0033694, lcolor(black) lwidth(1pt) lpattern(dash)) ///
	ytitle("Frequency") freq /// 
	xtitle("`=ustrunescape("\u03B8\u0302")'{subscript:b}")

graph export "${figdir}\figure_3-7.pdf", replace

erase "${datadir}\bootmed.dta"

log close

