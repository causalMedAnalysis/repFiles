/*figure 5.6*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/cmed //module to perform causal mediation analysis

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch5\_LOGS\"
global figdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\figures\ch5\" 

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/Brader_et_al2008/Brader_et_al2008.dta" ///
	"${datadir}Brader_et_al2008\"

//open log
log using "${logdir}figure_5-6.log", replace 

//load data
use "${datadir}Brader_et_al2008\Brader_et_al2008.dta", clear

//keep complete cases
drop if missing(immigr, emo, p_harm, tone_eth, ppage, ppeducat, ppgender, ppincimp)

//standardize outcome
egen std_immigr=std(4-immigr)

//dummy code controls
tab ppeducat, gen(ppeducat_)
drop ppeducat_1
rename ppeducat_2 hs
rename ppeducat_3 sc
rename ppeducat_4 ba

tab ppgender, gen(ppgender_)
drop ppgender_1
rename ppgender_2 female

//define macros for different variables
global C ppage female hs sc ba ppincimp //baseline confounders
global D tone_eth //exposure
global M1 p_harm //first mediator
global M2 emo //second mediator
global Y std_immigr //outcome

//compute pure regression imputation estimates w/o interactions
qui cmed impute ((regress) $Y) ($M1 $M2) $D = $C, paths nointer
	
scalar PSE_DY=_b[PSE_DY]

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

//create contour plot of bias-adjusted estimates
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
