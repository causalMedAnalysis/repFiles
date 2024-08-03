/*Figure 4.10*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"
*global figdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\figures\ch4\" 

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"
global figdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\figures\ch4\" 

log using "${logdir}figure_4-10.log", replace 

//input data
use "${datadir}plowUse\plowUse.dta", clear

global allvars ///
	women_politics plow ln_income ///
	agricultural_suitability tropical_climate large_animals rugged ///
	polity2_2000

foreach v in $allvars {
	drop if missing(`v')
	}

keep isocode $allvars

replace plow = round(plow)

replace women_politics = women_politics/100

recode polity2_2000 (-10/-6=1) (-5/5=2) (6/10=3), gen(authGovCat)
quietly tab authGovCat, gen(authGovCat_)

//compute RWR estimates

//compute bootstrap CIs
quietly bootstrap ///
	OE=_b[rATE] IDE=_b[rNDE] IIE=_b[rNIE] CDE=_b[CDE], ///
	reps(200) seed(60637): ///
		rwrmed women_politics authGovCat, ///
			avar(plow) mvar(ln_income) ///
			cvar(agricultural_suitability tropical_climate large_animals rugged) ///
			a(0) astar(1) m(7.5) mreg(regress)

scalar IDE=_b[IDE] 
scalar OE=_b[OE]

//specify range for sensitivity parameters under M-Y confounding
clear

set obs 1

gen delta_UYgivCDLM=.
gen delta_DUgivC=.

local counter=1 

forval i=-0.1(0.01)0.1 {
	forval j=-0.1(0.01)0.1 {
		
		quietly replace delta_UYgivCDLM=`i' if _n==`counter'
		quietly replace delta_DUgivC=`j'  if _n==`counter'
		
		local counter=`counter'+1
		 
		set obs `=_N+1'
	}
}

drop if delta_UYgivCDLM==.

//compute bias-adjusted estimates
gen ide_adj = IDE-(delta_UYgivCDLM*delta_DUgivC)
gen oe_adj = OE+(delta_UYgivCDLM*delta_DUgivC)

format ide_adj oe_adj %12.3f

//create contour plots of bias-adjusted estimates against sensitivity parameters
set scheme s2mono

twoway ///
	(contour ide_adj delta_DUgivC delta_UYgivCDLM, ///
		levels(10) minmax heatmap ///
		yscale(range(-0.1 0.1)) ///
		xscale(range(-0.1 0.1)) ///
		ylabel(-0.1(0.02)0.1) ///
		xlabel(-0.1(0.02)0.1) ///
		zlabel(#10, format(%12.3f)) ///
		title("A. Bias-adjusted IDE(1,0) Estimates") ///
		ytitle("{&delta}{subscript:DU|C}") ///
		xtitle("{&delta}{subscript:UY|C,D,L,M}") ///
		ztitle(""))

graph save "${figdir}\ide_plot.gph", replace

twoway ///
	(contour oe_adj delta_DUgivC delta_UYgivCDLM, ///
		levels(10) minmax heatmap ///
		yscale(range(-0.1 0.1)) ///
		xscale(range(-0.1 0.1)) ///
		ylabel(-0.1(0.02)0.1) ///
		xlabel(-0.1(0.02)0.1) ///
		zlabel(#10, format(%12.3f)) ///
		title("B. Bias-adjusted OE(1,0) Estimates") ///
		ytitle("{&delta}{subscript:DU|C}") ///
		xtitle("{&delta}{subscript:UY|C,D,L,M}") ///
		ztitle(""))

graph save "${figdir}\oe_plot.gph", replace

graph combine ///
	"${figdir}\ide_plot.gph" ///
	"${figdir}\oe_plot.gph", ///
	col(1) row(2) ysize(8) xsize(4.5) scheme(s2mono)

graph export "${figdir}\figure_4-10.pdf", replace

erase "${figdir}\ide_plot.gph"
erase "${figdir}\oe_plot.gph"

log close
