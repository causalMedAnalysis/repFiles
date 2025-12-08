/*Figure 4.10*/
capture clear all
capture log close
set more off

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/plowUse/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch4\_LOGS\"
global figdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\figures\ch4\" 

//open log
log using "${logdir}figure_4-10.log", replace 

//load data
use "${datadir}plowUse.dta"

//keep complete cases
drop if missing(women_politics, plow, ln_income, agricultural_suitability, ///
	tropical_climate, large_animals, rugged, polity2_2000)

//recode variables
replace plow = round(plow)

replace women_politics = women_politics/100

recode polity2_2000 (-10/-6=1) (-5/5=2) (6/10=3), gen(authGov)

//define macros for different variables
global C agricultural_suitability tropical_climate large_animals rugged //baseline confounders
global D plow //exposure
global L authGov //exposure-induced confounder
global M ln_income //mediator
global Y women_politics //outcome

//compute RWR estimates
qui cmed linear $Y $M ($L) $D = $C

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
gen ide_adj = IDE - (delta_UYgivCDLM * delta_DUgivC)
gen oe_adj = OE + (delta_UYgivCDLM * delta_DUgivC)

format ide_adj oe_adj %12.3f

//create contour plots of bias-adjusted estimates
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

graph save "${figdir}ide_plot.gph", replace

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

graph save "${figdir}oe_plot.gph", replace

graph combine ///
	"${figdir}ide_plot.gph" ///
	"${figdir}oe_plot.gph", ///
	col(1) row(2) ysize(8) xsize(4.5) scheme(s2mono)

graph export "${figdir}figure_4-10.pdf", replace

erase "${figdir}ide_plot.gph"
erase "${figdir}oe_plot.gph"

log close
