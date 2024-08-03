/*Table 5.5*/
capture clear all
capture log close
set more off

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

log using "${logdir}table_5-5.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==.  | log_faminc_adj_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

foreach x in female black hispan paredu parprof parinc_prank famsize afqt3 {
	quietly sum `x'
	replace `x' = `x'-r(mean)
	}
	
//linear model estimates w/o exposure-mediator interactions
capture program drop linpath
program define linpath, rclass

	//fit model for M1|C,D
	reg ever_unemp_age3539 att22 ///
		 female black hispan paredu parprof parinc_prank famsize afqt3
		
	local beta01 = _b[_cons]
	local beta21 = _b[att22]

	//fit model for M2|C,D 
	reg log_faminc_adj_age3539 att22 ///
		 female black hispan paredu parprof parinc_prank famsize afqt3
		
	local beta02 = _b[_cons]
	local beta22 = _b[att22]
	
	//fit model for Y|C,D,M1,M2 
	reg std_cesd_age40 att22 ever_unemp_age3539 log_faminc_adj_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3

	local gamma2 = _b[att22]
	local gamma31 = _b[ever_unemp_age3539]
	local gamma32 = _b[log_faminc_adj_age3539]

	//fit model for Y|C,D,M1
	reg std_cesd_age40 att22 ever_unemp_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	
	local lamda2 = _b[att22]
	local lamda31 = _b[ever_unemp_age3539]
	
	//compute effect estimates
	local MNDE = `gamma2'*(1-0)
	local MNIE = (`beta21'*`gamma31' + `beta22'*`gamma32')*(1-0)
	
	local NDE_M1 = `lamda2'*(1-0)
	local NIE_M1 = `beta21'*`lamda31'*(1-0)
	
	return scalar ATE = `NDE_M1' + `NIE_M1'
	return scalar PSE_DY = `MNDE'
	return scalar PSE_DM2Y = `NDE_M1' - `MNDE'
	return scalar PSE_DM1Y = `NIE_M1'
	
end

quietly linpath

di r(ATE)
di r(PSE_DY)
di r(PSE_DM2Y)
di r(PSE_DM1Y)

quietly bootstrap ///
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_AM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): linpath
		
mat list e(ci_percentile)

//linear model estimates w/ exposure-mediator interactions
capture program drop linpathx
program define linpathx, rclass

	//fit model for M1|C,D
	reg ever_unemp_age3539 att22 ///
		 female black hispan paredu parprof parinc_prank famsize afqt3
		
	local beta01 = _b[_cons]
	local beta21 = _b[att22]

	//fit model for M2|C,D 
	reg log_faminc_adj_age3539 att22 ///
		 female black hispan paredu parprof parinc_prank famsize afqt3
		
	local beta02 = _b[_cons]
	local beta22 = _b[att22]
	
	//fit model for Y|C,D,M1,M2 
	reg std_cesd_age40 ///
		c.att22##c.ever_unemp_age3539 ///
		c.att22##c.log_faminc_adj_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3

	local gamma2 = _b[c.att22]
	local gamma31 = _b[c.ever_unemp_age3539]
	local gamma41 = _b[c.att22#c.ever_unemp_age3539]
	local gamma32 = _b[c.log_faminc_adj_age3539]
	local gamma42 = _b[c.att22#c.log_faminc_adj_age3539]

	//fit model for Y|C,D,M1 
	reg std_cesd_age40 ///
		c.att22##c.ever_unemp_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3

	local lamda2 = _b[c.att22]
	local lamda31 = _b[c.ever_unemp_age3539]
	local lamda41 = _b[c.att22#c.ever_unemp_age3539]
	
	//compute effect estimates
	local MNDE = (`gamma2' + `gamma41'*(`beta01'+`beta21'*0) + `gamma42'*(`beta02'+`beta22'*0))*(1-0)
	local MNIE = (`beta21'*(`gamma31' + `gamma41'*1) + `beta22'*(`gamma32' + `gamma42'*1))*(1-0)
	
	local NDE_M1 = (`lamda2' + `lamda41'*(`beta01'+`beta21'*0))*(1-0)
	local NIE_M1 = (`beta21'*(`lamda31' + `lamda41'*1))*(1-0)
	
	return scalar ATE = `NDE_M1' + `NIE_M1'
	return scalar PSE_DY = `MNDE'
	return scalar PSE_DM2Y = `NDE_M1' - `MNDE'
	return scalar PSE_DM1Y = `NIE_M1'
	
end

quietly linpathx

di r(ATE)
di r(PSE_DY)
di r(PSE_DM2Y)
di r(PSE_DM1Y)

quietly bootstrap ///
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_DM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): linpathx
		
mat list e(ci_percentile)

//IPW estimates
capture program drop ipwpath
program define ipwpath, rclass

	//fit logit model for D|C
	logit att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	est store DModel_1

	//fit logit model for D|C,M1
	logit att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3 ///
		ever_unemp_age3539
	est store DModel_2
	
	//fit logit model for D|C,M1,M2
	logit att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3 ///
		log_faminc_adj_age3539 ever_unemp_age3539
	est store DModel_3

	//predict exposure probabilities
	est restore DModel_1
	predict phat_d_C, pr
	gen phat_dstar_C=1-phat_d_C

	est restore DModel_2
	predict phat_d_CM1, pr
	gen phat_dstar_CM1=1-phat_d_CM1

	est restore DModel_3
	predict phat_d_CM1M2, pr
	gen phat_dstar_CM1M2=1-phat_d_CM1M2
	
	//compute stabilized weights
	logit att22
	predict phat_d, pr
	gen phat_dstar=1-phat_d
	
	gen _sw1=phat_dstar/phat_dstar_C if att22==0
	gen _sw2=phat_d/phat_d_C if att22==1
	gen _sw3=(phat_dstar_CM1*phat_d)/(phat_d_CM1*phat_dstar_C) if att22==1
	gen _sw5=(phat_dstar_CM1M2*phat_d)/(phat_d_CM1M2*phat_dstar_C) if att22==1

	//censor stabilized weights at 1st and 99th percentiles
	foreach i of var _sw* {
		centile `i', c(1 99) 
		replace `i'=r(c_1) if `i'<r(c_1) & `i'!=.
		replace `i'=r(c_2) if `i'>r(c_2) & `i'!=.
		}

	//compute weighted means of Y
	reg std_cesd_age40 [pw=_sw1] if att22==0
	local Ehat_Y0 = _b[_cons]

	reg std_cesd_age40 [pw=_sw2] if att22==1
	local Ehat_Y1 = _b[_cons]

	reg std_cesd_age40 [pw=_sw3] if att22==1
	local Ehat_Y1M10 = _b[_cons]

	reg std_cesd_age40 [pw=_sw5] if att22==1
	local Ehat_Y1M10M20 = _b[_cons]
	
	//compute effect estimates
	
	local MNDE = `Ehat_Y1M10M20' - `Ehat_Y0'
	local MNIE = `Ehat_Y1' - `Ehat_Y1M10M20'

	local NDE_M1 = `Ehat_Y1M10' - `Ehat_Y0'
	local NIE_M1 = `Ehat_Y1' - `Ehat_Y1M10'
	
	return scalar ATE = `Ehat_Y1' - `Ehat_Y0'
	return scalar PSE_DY = `MNDE'
	return scalar PSE_DM2Y = `NDE_M1' - `MNDE'
	return scalar PSE_DM1Y = `NIE_M1'
	
	drop phat_d_* phat_dstar_* phat_d phat_dstar _sw* _est_DModel_* 		
	
	//reset estimation sample for bootstrap
	reg att22
	
end

quietly ipwpath

di r(ATE)
di r(PSE_DY)
di r(PSE_DM2Y)
di r(PSE_DM1Y)

quietly bootstrap ///
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_DM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): ipwpath
		
mat list e(ci_percentile)

log close
