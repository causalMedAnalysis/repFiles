/*Table 5.8*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

log using "${logdir}table_5-8.log", replace 

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

// EXPOSURE D = tone_eth
// MEDIATOR M1 = p_harm
// MEDIATOR M2 = emo
// OUTCOME Y = std_immigr
// BASELINE CONFOUNDERS C = ppage female hs sc ba ppincimp	

//pure regression imputation estimates w/o interactions
capture program drop pathimp
program define pathimp, rclass

	gen D_orig=tone_eth

	//compute regression imputations
	reg std_immigr tone_eth ppage female hs sc ba ppincimp

	replace tone_eth = 0

	predict	Y0hat, xb

	replace tone_eth = 1

	predict Y1hat, xb

	replace tone_eth = D_orig

	reg std_immigr tone_eth p_harm emo ppage female hs sc ba ppincimp

	replace tone_eth = 1

	predict Y1M1DM2Dhat, xb

	replace tone_eth = D_orig

	reg Y1M1DM2Dhat ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 0

	predict Y1M10M200hat, xb

	replace tone_eth = D_orig

	reg std_immigr tone_eth p_harm ppage female hs sc ba ppincimp

	replace tone_eth = 1

	predict Y1M1Dhat, xb

	replace tone_eth = D_orig

	reg Y1M1Dhat ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 0

	predict Y1M10M210hat, xb

	replace tone_eth = D_orig

	//average imputations over sample members
	reg Y0hat
	local Ehat_Y0=_b[_cons]
	drop Y0hat

	reg Y1hat
	local Ehat_Y1=_b[_cons]
	drop Y1hat

	reg Y1M10M200hat
	local Ehat_Y1M10M200=_b[_cons]
	drop Y1M10M200hat

	reg Y1M10M210hat
	local Ehat_Y1M10M210=_b[_cons]
	drop Y1M10M210hat

	//compute effect estimates
	return scalar ATE=`Ehat_Y1'-`Ehat_Y0'
	return scalar PSE_DY=`Ehat_Y1M10M200'-`Ehat_Y0'
	return scalar PSE_DM2Y=`Ehat_Y1M10M210'-`Ehat_Y1M10M200'
	return scalar PSE_DM1Y=`Ehat_Y1'-`Ehat_Y1M10M210'
		
	drop D_orig Y1M1DM2Dhat Y1M1Dhat
	
end

quietly pathimp

di r(ATE)
di r(PSE_DY)
di r(PSE_DM2Y)
di r(PSE_DM1Y)

quietly bootstrap ///
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_DM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): pathimp
		
mat list e(ci_percentile)

//pure regression imputation estimates w/ D x {M1,M2} interactions
capture program drop pathimpx
program define pathimpx, rclass

	gen D_orig=tone_eth

	//compute regression imputations
	reg std_immigr tone_eth ///
		ppage female hs sc ba ppincimp

	replace tone_eth = 0

	predict	Y0hat, xb

	replace tone_eth = 1

	predict Y1hat, xb

	replace tone_eth = D_orig

	reg std_immigr i.tone_eth##c.p_harm i.tone_eth##c.emo ///
		ppage female hs sc ba ppincimp

	replace tone_eth = 1

	predict Y1M1DM2Dhat, xb

	replace tone_eth = D_orig

	reg Y1M1DM2Dhat ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 0

	predict Y1M10M200hat, xb

	replace tone_eth = D_orig

	reg std_immigr i.tone_eth##c.p_harm ///
		ppage female hs sc ba ppincimp

	replace tone_eth = 1

	predict Y1M1Dhat, xb

	replace tone_eth = D_orig

	reg Y1M1Dhat ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 0

	predict Y1M10M210hat, xb

	replace tone_eth = D_orig

	//average imputations over sample members
	reg Y0hat
	local Ehat_Y0=_b[_cons]
	drop Y0hat

	reg Y1hat
	local Ehat_Y1=_b[_cons]
	drop Y1hat

	reg Y1M10M200hat
	local Ehat_Y1M10M200=_b[_cons]
	drop Y1M10M200hat

	reg Y1M10M210hat
	local Ehat_Y1M10M210=_b[_cons]
	drop Y1M10M210hat

	//compute effect estimates
	return scalar ATE=`Ehat_Y1'-`Ehat_Y0'
	return scalar PSE_DY=`Ehat_Y1M10M200'-`Ehat_Y0'
	return scalar PSE_DM2Y=`Ehat_Y1M10M210'-`Ehat_Y1M10M200'
	return scalar PSE_DM1Y=`Ehat_Y1'-`Ehat_Y1M10M210'
		
	drop D_orig Y1M1DM2Dhat Y1M1Dhat
	
end

quietly pathimpx

di r(ATE)
di r(PSE_DY)
di r(PSE_DM2Y)
di r(PSE_DM1Y)

quietly bootstrap ///
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_DM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): pathimpx
		
mat list e(ci_percentile)

//pure regression imputation estimates w/ D x {M1,M2,C} interactions
capture program drop pathimpxx
program define pathimpxx, rclass

	gen D_orig=tone_eth

	//compute regression imputations
	reg std_immigr ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 0

	predict	Y0hat, xb

	replace tone_eth = 1

	predict Y1hat, xb

	replace tone_eth = D_orig

	reg std_immigr ///
		i.tone_eth##c.p_harm ///
		i.tone_eth##c.emo ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 1

	predict Y1M1DM2Dhat, xb

	replace tone_eth = D_orig

	reg Y1M1DM2Dhat ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 0

	predict Y1M10M200hat, xb

	replace tone_eth = D_orig

	reg std_immigr ///
		i.tone_eth##c.p_harm ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 1

	predict Y1M1Dhat, xb

	replace tone_eth = D_orig

	reg Y1M1Dhat ///
		i.tone_eth##c.ppage ///
		i.tone_eth##i.female ///
		i.tone_eth##i.hs ///
		i.tone_eth##i.sc ///
		i.tone_eth##i.ba ///
		i.tone_eth##c.ppincimp

	replace tone_eth = 0

	predict Y1M10M210hat, xb

	replace tone_eth = D_orig

	//average imputations over sample members
	reg Y0hat
	local Ehat_Y0=_b[_cons]
	drop Y0hat

	reg Y1hat
	local Ehat_Y1=_b[_cons]
	drop Y1hat

	reg Y1M10M200hat
	local Ehat_Y1M10M200=_b[_cons]
	drop Y1M10M200hat

	reg Y1M10M210hat
	local Ehat_Y1M10M210=_b[_cons]
	drop Y1M10M210hat

	//compute effect estimates
	return scalar ATE=`Ehat_Y1'-`Ehat_Y0'
	return scalar PSE_DY=`Ehat_Y1M10M200'-`Ehat_Y0'
	return scalar PSE_DM2Y=`Ehat_Y1M10M210'-`Ehat_Y1M10M200'
	return scalar PSE_DM1Y=`Ehat_Y1'-`Ehat_Y1M10M210'
		
	drop D_orig Y1M1DM2Dhat Y1M1Dhat
	
end

quietly pathimpxx

di r(ATE)
di r(PSE_DY)
di r(PSE_DM2Y)
di r(PSE_DM1Y)

quietly bootstrap ///
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_DM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): pathimpxx
		
mat list e(ci_percentile)

log close
