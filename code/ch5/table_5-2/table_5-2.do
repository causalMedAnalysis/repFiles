/*Table 5.2*/
capture clear all
capture log close
set more off

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

log using "${logdir}table_5-2.log", replace 

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
	
/*linear model estimates w/o exposure-mediator interactions*/
capture program drop linmed
program define linmed, rclass

	//fit model for M1
	reg ever_unemp_age3539 att22 ///
		 female black hispan paredu parprof parinc_prank famsize afqt3
		
	scalar beta01 = _b[_cons]
	scalar beta21 = _b[att22]

	//fit model for M2 
	reg log_faminc_adj_age3539 att22 ///
		 female black hispan paredu parprof parinc_prank famsize afqt3
		
	scalar beta02 = _b[_cons]
	scalar beta22 = _b[att22]
	
	//fit model for Y 
	reg std_cesd_age40 att22 ever_unemp_age3539 log_faminc_adj_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3

	scalar gamma2 = _b[att22]
	scalar gamma31 = _b[ever_unemp_age3539]
	scalar gamma32 = _b[log_faminc_adj_age3539]
	
	//compute effect estimates
	return scalar MNDE = gamma2*(1-0)
	return scalar MNIE = (beta21*gamma31 + beta22*gamma32)*(1-0)
	return scalar ATE = (gamma2 + beta21*gamma31 + beta22*gamma32)*(1-0)
	
end

quietly linmed

di r(ATE)
di r(MNDE)
di r(MNIE)

quietly bootstrap ///
	ATE=r(ATE) MNDE=r(MNDE) MNIE=r(MNIE), ///
	reps(2000) seed(60637): linmed
		
mat list e(ci_percentile)

/*linear model estimates w/ exposure-mediator interactions*/
capture program drop linmedx
program define linmedx, rclass

	//fit model for M1
	reg ever_unemp_age3539 att22 ///
		 female black hispan paredu parprof parinc_prank famsize afqt3
		
	scalar beta01 = _b[_cons]
	scalar beta21 = _b[att22]

	//fit model for M2 
	reg log_faminc_adj_age3539 att22 ///
		 female black hispan paredu parprof parinc_prank famsize afqt3
		
	scalar beta02 = _b[_cons]
	scalar beta22 = _b[att22]
	
	//fit model for Y 
	reg std_cesd_age40 ///
		c.att22##c.ever_unemp_age3539 ///
		c.att22##c.log_faminc_adj_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3

	scalar gamma2 = _b[c.att22]
	
	scalar gamma31 = _b[c.ever_unemp_age3539]
	scalar gamma41 = _b[c.att22#c.ever_unemp_age3539]
	
	scalar gamma32 = _b[c.log_faminc_adj_age3539]
	scalar gamma42 = _b[c.att22#c.log_faminc_adj_age3539]
	
	//compute effect estimates
	return scalar MNDE = (gamma2 + gamma41*(beta01+beta21*0) + gamma42*(beta02+beta22*0))*(1-0)
	return scalar MNIE = (beta21*(gamma31 + gamma41*1) + beta22*(gamma32 + gamma42*1))*(1-0)
	return scalar ATE = (gamma2 + gamma41*(beta01+beta21*0) + gamma42*(beta02+beta22*0))*(1-0) + ///
		(beta21*(gamma31 + gamma41*1) + beta22*(gamma32 + gamma42*1))*(1-0)
	
end

quietly linmedx

di r(ATE)
di r(MNDE)
di r(MNIE)

quietly bootstrap ///
	ATE=r(ATE) MNDE=r(MNDE) MNIE=r(MNIE), ///
	reps(2000) seed(60637): linmedx
		
mat list e(ci_percentile)

// IPW estimates
capture program drop ipwmed
program define ipwmed, rclass

	//fit logit model for D given C
	logit att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	est store DModel_1

	//fit logit model for D given C, M1, and M2
	logit att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3 ///
		log_faminc_adj_age3539 ever_unemp_age3539
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
	return scalar MNDE=`Ehat_Y1M0'-`Ehat_Y0M0'
	return scalar MNIE=`Ehat_Y1M1'-`Ehat_Y1M0'
	
	drop phat_d_C phat_dstar_C phat_d_CM phat_dstar_CM phat_d phat_dstar _sw* ///
		_est_DModel_1 _est_DModel_2
	
	//reset estimation sample for bootstrap
	reg att22
	
end

quietly ipwmed

di r(ATE)
di r(MNDE)
di r(MNIE)

quietly bootstrap ///
	ATE=r(ATE) MNDE=r(MNDE) MNIE=r(MNIE), ///
	reps(2000) seed(3308004): ipwmed
	
mat list e(ci_percentile)

log close
