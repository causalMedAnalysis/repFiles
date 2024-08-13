/*figure 5.6*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"
global figdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\figures\ch5\" 

log using "${logdir}figure_5-6.log", replace 

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

//compute point estimate of direct effect
quietly pathimp

scalar PSE_DY=r(PSE_DY)

//compute reference values for sensitivity parameters
qui reg std_immigr tone_eth p_harm emo ppage female hs sc ba ppincimp
local UY_female = _b[female]
local UY_colgrad = _b[ba]

qui reg female tone_eth p_harm emo ppage hs sc ba ppincimp
local DU_female = _b[tone_eth]

qui reg ba tone_eth p_harm emo ppage female hs sc ppincimp
local DU_colgrad = _b[tone_eth]

//specify range for sensitivity parameters under M-Y confounding
clear

set obs 1

gen delta_UYgivCDM1M2=.
gen delta_DUgivCM1M2=.

local counter=1 

forval i=0.0(0.1)2.1 {
	forval j=-0.5(0.025)0.5 {
		
		quietly replace delta_UYgivCDM1M2=`i' if _n==`counter'
		quietly replace delta_DUgivCM1M2=`j'  if _n==`counter'
		
		local counter=`counter'+1
		 
		qui set obs `=_N+1'
	}
}

drop if delta_UYgivCDM1M2==.

//compute bias-adjusted estimates
gen pse_DY_adj = PSE_DY - (delta_UYgivCDM1M2 * delta_DUgivCM1M2)

format pse_DY_adj %12.3f

//create contour plot of bias-adjusted estimates against sensitivity parameters
set scheme s2mono

twoway ///
	(contour pse_DY_adj delta_UYgivCDM1M2 delta_DUgivCM1M2 , ///
		levels(10) minmax heatmap ///
		yscale(range(0.0 2.0)) ///
		xscale(range(-0.5 0.5)) ///
		ylabel(0.0(0.2)2.0) ///
		xlabel(-0.5(0.1)0.5) ///
		zlabel(#10, format(%12.3f)) ///
		title("Bias-adjusted PSE{subscript:D->Y}(1,0) Estimates") ///
		xtitle("{&delta}{subscript:DU|C,M1,M2}") ///
		ytitle("{&delta}{subscript:UY|C,D,M1,M2}") ///
		ztitle("")) ///
	(scatteri `UY_female' `DU_female' (3) "female", ///
		legend(off) mcolor(black) msymbol(circle) mlabcolor(black)) ///
	(scatteri `UY_colgrad' `DU_colgrad' (9) "college graduate", ///
		legend(off) mcolor(black) msymbol(circle) mlabcolor(black))

graph export "${figdir}\figure_5-6.pdf", replace

log close
