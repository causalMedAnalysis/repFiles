/*Table 5.7*/
capture clear all
capture log close
set more off

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

log using "${logdir}table_5-7.log", replace 

//input data
use "${datadir}Brader_et_al2008\Brader_et_al2008.dta", clear

drop if ///
	immigr==. | emo==. | p_harm==.  | tone_eth==. | ///
	ppage==. | ppeducat==. | ppgender==. | ppincimp==.
	
egen std_immigr=std(4-immigr)

tab ppeducat, gen(ppeducat_)
drop ppeducat_1
rename ppeducat_2 hs
rename ppeducat_3 sc
rename ppeducat_4 ba

recode ppgender (3=0) (4=1), gen(female)

foreach x in emo p_harm ppage female hs sc ba ppincimp {
	quietly sum `x'
	replace `x' = `x'-r(mean)
	}

// EXPOSURE D = tone_eth
// MEDIATOR M1 = p_harm
// MEDIATOR M2 = emo
// OUTCOME Y = std_immigr
// BASELINE CONFOUNDERS C = ppage female hs sc ba ppincimp	

//linear model estimates w/o exposure-mediator interactions
capture program drop linpath
program define linpath, rclass

	//fit model for M1
	reg p_harm tone_eth ppage female hs sc ba ppincimp
		
	local beta01 = _b[_cons]
	local beta21 = _b[tone_eth]

	//fit model for M2 
	reg emo tone_eth ppage female hs sc ba ppincimp
		
	local beta02 = _b[_cons]
	local beta22 = _b[tone_eth]
	
	//fit model for Y|C,D,M1,M2 
	reg std_immigr tone_eth p_harm emo ppage female hs sc ba ppincimp

	local gamma2 = _b[tone_eth]
	local gamma31 = _b[p_harm]
	local gamma32 = _b[emo]

	//fit model for Y|C,D,M1
	reg std_immigr tone_eth p_harm ppage female hs sc ba ppincimp
	
	local lamda2 = _b[tone_eth]
	local lamda31 = _b[p_harm]
	
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
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_DM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): linpath
		
mat list e(ci_percentile)

//linear model estimates w/ exposure-mediator interactions
capture program drop linpathx
program define linpathx, rclass

	//fit model for M1|C,D
	reg p_harm tone_eth ppage female hs sc ba ppincimp
		
	local beta01 = _b[_cons]
	local beta21 = _b[tone_eth]

	//fit model for M2|C,D 
	reg emo tone_eth ppage female hs sc ba ppincimp
		
	local beta02 = _b[_cons]
	local beta22 = _b[tone_eth]
		
	//fit model for Y|C,D,M1,M2 
	reg std_immigr ///
		c.tone_eth##c.p_harm ///
		c.tone_eth##c.emo ///
		ppage female hs sc ba ppincimp
	
	local gamma2 = _b[c.tone_eth]
	local gamma31 = _b[c.p_harm]
	local gamma41 = _b[c.tone_eth#c.p_harm]
	local gamma32 = _b[c.emo]
	local gamma42 = _b[c.tone_eth#c.emo]

	//fit model for Y|C,D,M1 
	reg std_immigr ///
		c.tone_eth##c.p_harm ///
		ppage female hs sc ba ppincimp

	local lamda2 = _b[c.tone_eth]
	local lamda31 = _b[c.p_harm]
	local lamda41 = _b[c.tone_eth#c.p_harm]
	
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
	logit tone_eth ppage female hs sc ba ppincimp
	est store DModel_1

	//fit logit model for D|C,M1
	logit tone_eth p_harm ppage female hs sc ba ppincimp
	est store DModel_2
	
	//fit logit model for D|C,M1,M2
	logit tone_eth p_harm emo ppage female hs sc ba ppincimp
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
	logit tone_eth
	predict phat_d, pr
	gen phat_dstar=1-phat_d
	
	gen _sw1=phat_dstar/phat_dstar_C if tone_eth==0
	gen _sw2=phat_d/phat_d_C if tone_eth==1
	gen _sw3=(phat_dstar_CM1*phat_d)/(phat_d_CM1*phat_dstar_C) if tone_eth==1
	gen _sw5=(phat_dstar_CM1M2*phat_d)/(phat_d_CM1M2*phat_dstar_C) if tone_eth==1

	//censor stabilized weights at 1st and 99th percentiles
	foreach i of var _sw* {
		centile `i', c(1 99) 
		replace `i'=r(c_1) if `i'<r(c_1) & `i'!=.
		replace `i'=r(c_2) if `i'>r(c_2) & `i'!=.
		}

	//compute weighted means of Y
	reg std_immigr [pw=_sw1] if tone_eth==0
	local Ehat_Y0 = _b[_cons]

	reg std_immigr [pw=_sw2] if tone_eth==1
	local Ehat_Y1 = _b[_cons]

	reg std_immigr [pw=_sw3] if tone_eth==1
	local Ehat_Y1M10 = _b[_cons]

	reg std_immigr [pw=_sw5] if tone_eth==1
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
	reg tone_eth
	
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
