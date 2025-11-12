/*Figure 4.9*/
capture clear all
capture log close
set more off
set maxvar 50000, perm

//install required modules
net install github, from("https://haghish.github.io/github/")

net install parallel, from("https://raw.github.com/gvegayon/parallel/stable/") replace 
mata mata mlib index

github install causalMedAnalysis/cmed

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch4\_LOGS\"
global figdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\figures\ch4\" 

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/plowUse/plowUse.dta" ///
	"${datadir}plowUse\"

//open log
log using "${logdir}figure_4-9.log", replace 

//load data
use "${datadir}plowUse\plowUse.dta"

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
qui cmed linear $Y $M ($L) $D = $C, reps(2000) seed(60637)

matrix b_matrix = e(b)'
matrix ci_matrix = e(ci_percentile)'
matrix rwrIE = b_matrix, ci_matrix

qui cmed linear $Y $M ($L) $D = $C, m(7.5) reps(2000) seed(60637)

matrix b_matrix = e(b)'
matrix ci_matrix = e(ci_percentile)'
matrix rwrCDE = b_matrix, ci_matrix

matrix rwrResults = rwrIE \ rwrCDE

//compute IPW estimates
cmed ipw $Y ((regress) $M) ((ologit) $L) $D = $C, censor(3 97) reps(2000) seed(60637)

matrix b_matrix = e(b)'
matrix ci_matrix = e(ci_percentile)'
matrix ipwIE = b_matrix, ci_matrix

cmed ipw $Y ((regress) $M) ($L) $D = $C, m(7.5) censor(3 97) reps(2000) seed(60637)

matrix b_matrix = e(b)'
matrix ci_matrix = e(ci_percentile)'
matrix ipwCDE = b_matrix, ci_matrix

matrix ipwResults = ipwIE \ ipwCDE

//compute simulation estimates w/ custom function for fitting ordinal logit 
//model for L and logistic binomial model for Y
capture program drop custsim
program define custsim, rclass

	syntax [, nsim(integer 10)] 

	ologit $L $D $C
	est store Lmodel

	reg $M $D $C
	est store Mmodel
	
	glm $Y i.$D##c.$M $L $C, family(binomial) link(logit)
	est store Ymodel

	gen D_orig=$D
	gen L_orig=$L
	gen M_orig=$M

	forval i=1/`nsim' {
		
		replace $D=0
		
		est restore Lmodel
		predict phat1_L0 phat2_L0 phat3_L0, pr
		gen runif=runiform()
		gen L0_`i'=1 if runif<=phat1_L0
		replace L0_`i'=2 if runif>phat1_L0 & runif<=(phat1_L0+phat2_L0)
		replace L0_`i'=3 if runif>(phat1_L0+phat2_L0)
		drop runif

		est restore Mmodel
		predict mhat_M0, xb
		gen M0_`i'=rnormal(mhat_M0,e(rmse))

		replace $D=1
		
		est restore Lmodel
		predict phat1_L1 phat2_L1 phat3_L1, pr
		gen runif=runiform()
		gen L1_`i'=1 if runif<=phat1_L1
		replace L1_`i'=2 if runif>phat1_L1 & runif<=(phat1_L1+phat2_L1)
		replace L1_`i'=3 if runif>(phat1_L1+phat2_L1)
		drop runif
		
		est restore Mmodel
		predict mhat_M1, xb
		gen M1_`i'=rnormal(mhat_M1,e(rmse))

		est restore Ymodel

		replace $D=0
		replace $L=L0_`i'
		replace $M=M0_`i'
		
		predict phat_Y0L0M0
		gen Y0L0M0_`i'=rbinomial(100,phat_Y0L0M0)/100 

		replace $D=1
		replace $L=L1_`i'
		replace $M=M1_`i'
		
		predict phat_Y1L1M1
		gen Y1L1M1_`i'=rbinomial(100,phat_Y1L1M1)/100

		replace $M=M0_`i'

		predict phat_Y1L1M0
		gen Y1L1M0_`i'=rbinomial(100,phat_Y1L1M0)/100

	    replace $M=7.5
	
		predict phat_Y1L1m
		gen Y1L1m_`i'=rbinomial(100,phat_Y1L1m)/100

		replace $D=0
		replace $L=L0_`i'
		
		predict phat_Y0L0m
		gen Y0L0m_`i'=rbinomial(100,phat_Y0L0m)/100 
		
		drop phat* mhat* M0_* M1_* L0_* L1_*
	}

	egen Y0L0M0_=rowmean(Y0L0M0_*)
	egen Y1L1M1_=rowmean(Y1L1M1_*)
	egen Y1L1M0_=rowmean(Y1L1M0_*)
	egen Y1L1m_=rowmean(Y1L1m_*)
	egen Y0L0m_=rowmean(Y0L0m_*)

	reg Y0L0M0_
	local Ehat_Y0L0M0=_b[_cons]
	drop Y0L0M0*

	reg Y1L1M1_
	local Ehat_Y1L1M1=_b[_cons]
	drop Y1L1M1*

	reg Y1L1M0_
	local Ehat_Y1L1M0=_b[_cons]
	drop Y1L1M0*

	reg Y1L1m_
	local Ehat_Y1L1m=_b[_cons]
	drop Y1L1m*

	reg Y0L0m_
	local Ehat_Y0L0m=_b[_cons]
	drop Y0L0m*

	return scalar oe=`Ehat_Y1L1M1'-`Ehat_Y0L0M0'
	return scalar ide=`Ehat_Y1L1M0'-`Ehat_Y0L0M0'
	return scalar iie=`Ehat_Y1L1M1'-`Ehat_Y1L1M0'
	return scalar cde=`Ehat_Y1L1m'-`Ehat_Y0L0m'

	replace $D=D_orig
	replace $L=L_orig
	replace $M=M_orig

	drop D_orig L_orig M_orig _est_Lmodel _est_Mmodel _est_Ymodel

end

//parallelize the bootstrap to reduce wall time
parallel initialize 18 //specify number of logical cores for each child process

parallel bs, expr(OE=r(oe) IDE=r(ide) IIE=r(iie) CDE=r(cde)) ///
	reps(2000): custsim, nsim(2000)

matrix b_matrix = e(b)'
matrix ci_matrix = e(ci_percentile)'
matrix simResults = b_matrix, ci_matrix

//create dot-whisker plot
//plot RWR estimates
clear
svmat rwrResults, names(col)
rename y1 pointEst

input estimand
	1
	2
	3
	4

label define estLbl ///
	1 "OE(1,0)" ///
	2 "IDE(1,0)" ///
	3 "IIE(1,0)" ///
	4 "CDE(1,0,1.8K)" 

label values estimand estLbl

set scheme s2mono

twoway (rcap ul ll estimand, horizontal) ///
       (scatter estimand pointEst, mcolor(black) msize(small)), ///
       legend(off) ytitle("") xtitle("Standard Deviations") ///
       yscale(reverse) xline(0) ///
       ylabel(1/4, valuelabel angle(horizontal) noticks nogrid) ///
       xtick(-0.1(0.02)0.06) xlabel(-0.1(0.02)0.06) ///
       title("RWR Estimates", size(medium))

graph save "${figdir}\rwr_plot.gph", replace

//plot simulation estimates
clear
svmat simResults, names(col)
rename y1 pointEst

input estimand
	1
	2
	3
	4

label define estLbl ///
	1 "OE(1,0)" ///
	2 "IDE(1,0)" ///
	3 "IIE(1,0)" ///
	4 "CDE(1,0,1.8K)" 

label values estimand estLbl

set scheme s2mono

twoway (rcap ul ll estimand, horizontal) ///
       (scatter estimand pointEst, mcolor(black) msize(small)), ///
       legend(off) ytitle("") xtitle("Standard Deviations") ///
       yscale(reverse) xline(0) ///
       ylabel(1/4, valuelabel angle(horizontal) noticks nogrid) ///
       xtick(-0.1(0.02)0.06) xlabel(-0.1(0.02)0.06) ///
       title("Simulation Estimates", size(medium))

graph save "${figdir}\sim_plot.gph", replace

//plot ipw estimates
clear
svmat ipwResults, names(col)
rename y1 pointEst

input estimand
	1
	2
	3
	4

label define estLbl ///
	1 "OE(1,0)" ///
	2 "IDE(1,0)" ///
	3 "IIE(1,0)" ///
	4 "CDE(1,0,1.8K)" 

label values estimand estLbl

set scheme s2mono

twoway (rcap ul ll estimand, horizontal) ///
       (scatter estimand pointEst, mcolor(black) msize(small)), ///
       legend(off) ytitle("") xtitle("Difference in Proportions") ///
       yscale(reverse) xline(0) ///
       ylabel(1/4, valuelabel angle(horizontal) noticks nogrid) ///
       xtick(-0.1(0.02)0.06) xlabel(-0.1(0.02)0.06) ///
       title("IPW Estimates", size(medium))

graph save "${figdir}\ipw_plot.gph", replace

//combine plots
graph combine ///
	"${figdir}\rwr_plot.gph" ///
	"${figdir}\sim_plot.gph" ///
	"${figdir}\ipw_plot.gph", ///
	col(3) row(1) ysize(4.5) xsize(10) scheme(s2mono)

graph export "${figdir}\figure_4-9.pdf", replace

//erase temporary figures
erase "${figdir}\rwr_plot.gph"
erase "${figdir}\sim_plot.gph"
erase "${figdir}\ipw_plot.gph"

log close

//some of the estimates differ slightly from those reported in the text, 
//which are based on the R implementation. This is variously due to minor 
//differences in how the weights are censored, to Monte Carlo error, or to 
//differences in random number seeding that influence the bootstrap samples

//note: -cmed sim- with large nsim() and reps() can be time consuming to run in 
//Stata because it cannot parallelize the bootstrap replications. Consider 
//switching to R implementation of -medsim- if computation time is a concern.
