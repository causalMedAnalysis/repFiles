/*Table 3.5*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

log using "${logdir}table_3-5.log", replace 

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

//compute 95% bootstrap CIs
quietly bootstrap ///
	ATE=r(ATE) NDE=r(NDE) NIE=r(NIE), ///
	saving("${datadir}\bootmed.dta", replace) ///
	reps(2000) seed(3308004): ipwmed
	
mat list e(ci_percentile)

//compute p-values for null of no effect
use "${datadir}\bootmed.dta", clear

foreach x in ATE NDE NIE {
	quietly gen `x'_lt=0 
	quietly replace `x'_lt=1 if `x'<0
	quietly sum `x'_lt
	scalar p_lt = r(mean)
	
	quietly gen `x'_rt=0 
	quietly replace `x'_rt=1 if `x'>0
	quietly sum `x'_rt
	scalar p_rt = r(mean)
	
	di "`x' p-value"
	di min(p_rt, p_lt)*2
	di "-----------"
	
	drop `x'_*
	}

//reload data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

//CDEhat^ipw

//define function for ipw estimator for CDE
capture program drop ipwcde
program define ipwcde, rclass

	syntax [, mval(real 0)] 
	
	//fit logit model for D given C
	logit att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	est store DModel
	
	//fit logit model for M given D and C
	logit ever_unemp_age3539 ///
		 att22 female black hispan paredu parprof parinc_prank famsize afqt3
	est store MModel

	//predict exposure probabilities
	est restore DModel
	predict phat_d_C, pr
	gen phat_D_denom=(att22*phat_d_C)+((1-att22)*(1-phat_d_C))

	est restore MModel
	predict phat_m_CD, pr
	gen phat_M_denom=(ever_unemp_age3539*phat_m_CD)+((1-ever_unemp_age3539)*(1-phat_m_CD))

	//compute stabilized weights
	logit att22
	predict phat_d, pr
	gen phat_D_num=(att22*phat_d)+((1-att22)*(1-phat_d))
	
	logit ever_unemp_age3539 att22
	predict phat_m_D, pr
	gen phat_M_num=(ever_unemp_age3539*phat_m_D)+((1-ever_unemp_age3539)*(1-phat_m_D))
	
	gen _sw4=(phat_M_num*phat_D_num)/(phat_M_denom*phat_D_denom)

	//censor stabilized weights at 1st and 99th percentiles
	centile _sw4, c(1 99) 
	replace _sw4=r(c_1) if _sw4<r(c_1) & _sw4!=.
	replace _sw4=r(c_2) if _sw4>r(c_2) & _sw4!=.

	//fit outcome model
	reg std_cesd_age40 i.att22##i.ever_unemp_age3539 [pw=_sw4]

	//compute effect estimates
	return scalar CDE=_b[1.att22]+_b[1.att22#1.ever_unemp_age3539]*`mval'
	
	drop phat* _est_* _sw*
	
end

//compute 95% bootstrap CIs
quietly bootstrap CDE=r(CDE), ///
	saving("${datadir}\bootcde.dta", replace) ///
	reps(2000) seed(3308004): ipwcde, mval(0)
	
mat list e(ci_percentile)

//compute p-values for null of no effect
use "${datadir}\bootcde.dta", clear

quietly gen CDE_lt=0 
quietly replace CDE_lt=1 if CDE<0
quietly sum CDE_lt
scalar p_lt = r(mean)
	
gen CDE_rt=0 
quietly replace CDE_rt=1 if CDE>0
quietly sum CDE_rt
scalar p_rt = r(mean)

di "CDE p-value"
di min(p_rt, p_lt)*2
di "-----------"
	
erase "${datadir}\bootmed.dta"
erase "${datadir}\bootcde.dta"

log close

