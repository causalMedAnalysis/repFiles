/*Table 3.4*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

log using "${logdir}table_3-4.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

///inverse probability weighting (ipw) estimators

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
	
end

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

//compute estimates
quietly ipwmed

di r(ATE)
di r(NDE)
di r(NIE)

quietly ipwcde, mval(0)

di r(CDE)	

log close
