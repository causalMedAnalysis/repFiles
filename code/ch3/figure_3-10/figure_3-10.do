/*Figure 3.10*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"
*global figdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\figures\ch3\" 

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"
global figdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\figures\ch3\" 

log using "${logdir}figure_3-10.log", replace 

//input data
use "${datadir}JOBSII\Jobs-NoMiss-Binary.dta", clear

//compute point estimates using linear models w/ exposure-mediator interaction
*ssc install rwrmed 

quietly rwrmed work1, ///
	avar(treat) mvar(job_seek) cvar(econ_hard sex age nonwhite educ income) ///
	a(0) astar(1) m(4) mreg(regress)

scalar NDE=_b[NDE] 
scalar NIE=_b[NIE]

//specify range for sensitivity parameters under M-Y confounding
clear

set obs 1

gen delta_UYgivCDM=.
gen delta_DUgivCM=.

local counter=1 

forval i=-0.1(0.01)0.1 {
	forval j=-0.1(0.01)0.1 {
		
		quietly replace delta_UYgivCDM=`i' if _n==`counter'
		quietly replace delta_DUgivCM=`j'  if _n==`counter'
		
		local counter=`counter'+1
		 
		set obs `=_N+1'
	}
}

drop if delta_UYgivCDM==.

//compute bias-adjusted estimates
gen nde_adj = NDE-(delta_UYgivCDM*delta_DUgivCM)
gen nie_adj = NIE+(delta_UYgivCDM*delta_DUgivCM)

format nde_adj nie_adj %12.3f

//create contour plots of bias-adjusted estimates against sensitivity parameters
set scheme s2mono

twoway ///
	(contour nde_adj delta_DUgivCM delta_UYgivCDM, ///
		levels(10) minmax heatmap ///
		yscale(range(-0.1 0.1)) ///
		xscale(range(-0.1 0.1)) ///
		ylabel(-0.1(0.02)0.1) ///
		xlabel(-0.1(0.02)0.1) ///
		zlabel(#10, format(%12.3f)) ///
		title("A. Bias-adjusted NDE(1,0) Estimates") ///
		ytitle("{&delta}{subscript:DU|C,M}") ///
		xtitle("{&delta}{subscript:UY|C,D,M}") ///
		ztitle(""))

graph save "${figdir}\nde_plot.gph", replace

twoway ///
	(contour nie_adj delta_DUgivCM delta_UYgivCDM, ///
		levels(10) minmax heatmap ///
		yscale(range(-0.1 0.1)) ///
		xscale(range(-0.1 0.1)) ///
		ylabel(-0.1(0.02)0.1) ///
		xlabel(-0.1(0.02)0.1) ///
		zlabel(#10, format(%12.3f)) ///
		title("B. Bias-adjusted NIE(1,0) Estimates") ///
		ytitle("{&delta}{subscript:DU|C,M}") ///
		xtitle("{&delta}{subscript:UY|C,D,M}") ///
		ztitle(""))

graph save "${figdir}\nie_plot.gph", replace

graph combine ///
	"${figdir}\nde_plot.gph" ///
	"${figdir}\nie_plot.gph", ///
	col(1) row(2) ysize(8) xsize(4.5) scheme(s2mono)

graph export "${figdir}\figure_3-10.pdf", replace

erase "${figdir}\nde_plot.gph"
erase "${figdir}\nie_plot.gph"

log close
