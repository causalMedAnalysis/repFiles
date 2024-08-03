/*Figure 4.9*/
capture clear all
capture log close
set more off
set maxvar 40000

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"
*global figdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\figures\ch4\" 

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"
global figdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\figures\ch4\" 

set seed 60637

log using "${logdir}figure_4-9.log", replace 

*ssc install rwrmed 

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

///RWR estimates

//compute bootstrap CIs
quietly bootstrap ///
	OE=_b[rATE] IDE=_b[rNDE] IIE=_b[rNIE] CDE=_b[CDE], ///
	reps(1000) seed(60637): ///
		rwrmed women_politics authGovCat, ///
			avar(plow) mvar(ln_income) ///
			cvar(agricultural_suitability tropical_climate large_animals rugged) ///
			a(0) astar(1) m(7.5) mreg(regress)

matrix b_matrix = e(b)'
matrix ci_matrix = e(ci_percentile)'
matrix rwrResults = b_matrix, ci_matrix

///simulation estimates

//define function for simulation estimators w/ all two-way interactions
capture program drop medsim
program define medsim, rclass

	syntax [, nsim(integer 10)] 
	
	//fit ordinal logit model for L 
	ologit authGovCat plow agricultural_suitability tropical_climate large_animals rugged
		
	est store Lmodel

	//fit normal linear model for M 
	reg ln_income plow agricultural_suitability tropical_climate large_animals rugged
		
	est store Mmodel
	
	//fit log binomial model for Y 
	glm women_politics i.plow##c.ln_income authGovCat agricultural_suitability tropical_climate large_animals rugged, family(binomial) link(log)
		
	est store Ymodel

	gen D_orig=plow
	gen L_orig=authGovCat
	gen M_orig=ln_income

	//simulate potential values of exposure-induced confounder, mediator, and outcome
	forval i=1/`nsim' {
		
		replace plow=0
		
		est restore Lmodel
		predict phat1_L0 phat2_L0 phat3_L0, pr
		gen runif=runiform()
		gen L0_`i'=1 if runif<=phat1_L0 //MC draws from Multinom(phat1_L0,phat2_L0,phat3_L0)
		replace L0_`i'=2 if runif>phat1_L0 & runif<=(phat1_L0+phat2_L0)
		replace L0_`i'=3 if runif>(phat1_L0+phat2_L0)
		drop runif

		est restore Mmodel
		predict mhat_M0, xb
		gen M0_`i'=rnormal(mhat_M0,e(rmse)) //MC draws from N(mu=mhat_M0,sd=e(rmse))

		replace plow=1
		
		est restore Lmodel
		predict phat1_L1 phat2_L1 phat3_L1, pr
		gen runif=runiform()
		gen L1_`i'=1 if runif<=phat1_L1 //MC draws from Multinom(phat1_L1,phat2_L1,phat3_L1)
		replace L1_`i'=2 if runif>phat1_L1 & runif<=(phat1_L1+phat2_L1)
		replace L1_`i'=3 if runif>(phat1_L1+phat2_L1)
		drop runif
		
		est restore Mmodel
		predict mhat_M1, xb
		gen M1_`i'=rnormal(mhat_M1,e(rmse)) //MC draws from N(mu=mhat_M1,sd=e(rmse))

		est restore Ymodel

		replace plow=0
		replace authGovCat=L0_`i'
		replace ln_income=M0_`i'
		
		predict yhat_Y0L0M0, xb
		gen phat_Y0L0M0=exp(yhat_Y0L0M0)
		gen Y0L0M0_`i'=rbinomial(100,phat_Y0L0M0)/100 //MC draws from binomial(p=phat_Y0L0M0)

		replace plow=1
		replace authGovCat=L1_`i'
		replace ln_income=M1_`i'
		
		predict yhat_Y1L1M1, xb
		gen phat_Y1L1M1=exp(yhat_Y1L1M1)
		gen Y1L1M1_`i'=rbinomial(100,phat_Y1L1M1)/100 //MC draws from binomial(p=phat_Y1L1M1)

		replace ln_income=M0_`i'

		predict yhat_Y1L1M0, xb
		gen phat_Y1L1M0=exp(yhat_Y1L1M0)
		gen Y1L1M0_`i'=rbinomial(100,phat_Y1L1M0)/100 //MC draws from binomial(p=phat_Y1L1M0)

	    replace ln_income=7.5
	
		predict yhat_Y1L1m, xb
		gen phat_Y1L1m=exp(yhat_Y1L1m)
		gen Y1L1m_`i'=rbinomial(100,phat_Y1L1m)/100 //MC draws from binomial(p=phat_Y1L1m)

		replace plow=0
		replace authGovCat=L0_`i'
		
		predict yhat_Y0L0m, xb
		gen phat_Y0L0m=exp(yhat_Y0L0m)
		gen Y0L0m_`i'=rbinomial(100,phat_Y0L0m)/100 //MC draws from binomial(p=phat_Y0L0m)
		
		drop phat* yhat* mhat* M0_* M1_* L0_* L1_*
		}

	//average simulated outcomes over simulations
	egen Y0L0M0_=rowmean(Y0L0M0_*)
	egen Y1L1M1_=rowmean(Y1L1M1_*)
	egen Y1L1M0_=rowmean(Y1L1M0_*)
	egen Y1L1m_=rowmean(Y1L1m_*)
	egen Y0L0m_=rowmean(Y0L0m_*)
	
	//average simulated outcomes over sample members
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
	
	//compute effect estimates
	return scalar OE=`Ehat_Y1L1M1'-`Ehat_Y0L0M0'
	return scalar IDE=`Ehat_Y1L1M0'-`Ehat_Y0L0M0'
	return scalar IIE=`Ehat_Y1L1M1'-`Ehat_Y1L1M0'
	return scalar CDE=`Ehat_Y1L1m'-`Ehat_Y0L0m'
	
	//restore original data
	replace plow=D_orig
	replace authGovCat=L_orig
	replace ln_income=M_orig
	
	drop D_orig L_orig M_orig _est_Lmodel _est_Mmodel _est_Ymodel
	
end

//compute bootstrap CIs
quietly bootstrap ///
	OE=r(OE) IDE=r(IDE) IIE=r(IIE) CDE=r(CDE), ///
	reps(1000) seed(60637): ///
		medsim, nsim(1000)
		
matrix b_matrix = e(b)'
matrix ci_matrix = e(ci_percentile)'
matrix simResults = b_matrix, ci_matrix

///ipw estimates

//define function for ipw estimator for OE, IDE, IIE, and CDE
capture program drop ipwmedx
program define ipwmedx, rclass 
	
	gen D_orig=plow
	gen L_orig=authGovCat
	gen M_orig=ln_income
	
	//fit models
	logit plow agricultural_suitability tropical_climate large_animals rugged
	
	est store Dmodel

	ologit authGovCat plow agricultural_suitability tropical_climate large_animals rugged
	
	est store Lmodel
	
	reg ln_income plow authGovCat agricultural_suitability tropical_climate large_animals rugged

	est store Mmodel

	//predict D probabilities
	est restore Dmodel
	
	predict phatD1_C, pr
	gen phatD0_C=1-phatD1_C
	
	//predict L probabilities
	est restore Lmodel
	
	replace plow=1
	
	predict phatL1_D1C phatL2_D1C phatL3_D1C, pr
	
	replace plow=0

	predict phatL1_D0C phatL2_D0C phatL3_D0C, pr

	//predict M probabilities
	est restore Mmodel
	
	replace plow=1
	
	predict mhat_D1LC, xb
	gen phatM_D1LC=normalden(ln_income, mhat_D1LC, e(rmse))
	
	replace plow=0
	
	predict mhat_D0LC, xb
	gen phatM_D0LC=normalden(ln_income, mhat_D0LC, e(rmse))

	replace authGovCat=1
	
	predict mhat_D0L1C, xb
	gen phatM_D0L1C=normalden(ln_income, mhat_D0L1C, e(rmse))

	replace authGovCat=2
	
	predict mhat_D0L2C, xb
	gen phatM_D0L2C=normalden(ln_income, mhat_D0L2C, e(rmse))

	replace authGovCat=3
	
	predict mhat_D0L3C, xb
	gen phatM_D0L3C=normalden(ln_income, mhat_D0L3C, e(rmse))

	replace plow=1
	
	replace authGovCat=1
	
	predict mhat_D1L1C, xb
	gen phatM_D1L1C=normalden(ln_income, mhat_D1L1C, e(rmse))

	replace authGovCat=2
	
	predict mhat_D1L2C, xb
	gen phatM_D1L2C=normalden(ln_income, mhat_D1L2C, e(rmse))

	replace authGovCat=3
	
	predict mhat_D1L3C, xb
	gen phatM_D1L3C=normalden(ln_income, mhat_D1L3C, e(rmse))

	replace plow=D_orig
	replace authGovCat=L_orig
	replace ln_income=M_orig

	//compute stabilized weights
	logit plow
	predict phatD1, pr
	gen phatD0=1-phatD1
	
	reg ln_income plow 
	
	replace plow=0
	
	predict mhat_D0, xb
	gen phatM_D0=normalden(ln_income, mhat_D0, e(rmse))

	replace plow=1
	
	predict mhat_D1, xb
	gen phatM_D1=normalden(ln_income, mhat_D1, e(rmse))

	replace plow=D_orig
	replace authGovCat=L_orig
	replace ln_income=M_orig

	gen _sw1 = (phatD0/phatD0_C) * (1/phatM_D0LC) * ((phatM_D0L1C * phatL1_D0C) + (phatM_D0L2C * phatL2_D0C) + (phatM_D0L3C * phatL3_D0C)) if plow==0 
	gen _sw2 = (phatD1/phatD1_C) * (1/phatM_D1LC) * ((phatM_D1L1C * phatL1_D1C) + (phatM_D1L2C * phatL2_D1C) + (phatM_D1L3C * phatL3_D1C)) if plow==1
	gen _sw3 = (phatD1/phatD1_C) * (1/phatM_D1LC) * ((phatM_D0L1C * phatL1_D0C) + (phatM_D0L2C * phatL2_D0C) + (phatM_D0L3C * phatL3_D0C)) if plow==1
	    
	gen _sw4 = (phatD1/phatD1_C) * (phatM_D1/phatM_D1LC) if plow==1
	replace _sw4 = (phatD0/phatD0_C) * (phatM_D0/phatM_D0LC) if plow==0

	//censor stabilized weights at 1st and 99th percentiles
	foreach i of var _sw* {
		centile `i', c(3 97) 
		replace `i'=r(c_1) if `i'<r(c_1) & `i'!=.
		replace `i'=r(c_2) if `i'>r(c_2) & `i'!=.
		}

	//compute weighted means of Y
	reg women_politics [pw=_sw1] if plow==0
	local Ehat_Y0M0=_b[_cons]

	reg women_politics [pw=_sw2] if plow==1
	local Ehat_Y1M1=_b[_cons]

	reg women_politics [pw=_sw3] if plow==1
	local Ehat_Y1M0=_b[_cons]

	//compute interventional effect estimates
	return scalar OE = `Ehat_Y1M1' - `Ehat_Y0M0'
	return scalar IDE = `Ehat_Y1M0' - `Ehat_Y0M0'
	return scalar IIE = `Ehat_Y1M1' - `Ehat_Y1M0'
	
	//compute CDE estimate
	reg women_politics i.plow##c.ln_income [pw=_sw4]
	
	return scalar CDE = _b[1.plow] + 7.5*_b[1.plow#c.ln_income]
	
	drop phat* mhat* _sw* _est* D_orig L_orig M_orig
	
end

//compute bootstrap CIs
quietly bootstrap ///
	OE=r(OE) IDE=r(IDE) IIE=r(IIE) CDE=r(CDE), ///
	reps(1000) seed(60637): ///
		ipwmedx

matrix b_matrix = e(b)'
matrix ci_matrix = e(ci_percentile)'
matrix ipwResults = b_matrix, ci_matrix

///display results in dot-whisker plot

//rwr estimates
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

//simulation estimates
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

//ipw estimates
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

erase "${figdir}\rwr_plot.gph"
erase "${figdir}\sim_plot.gph"
erase "${figdir}\ipw_plot.gph"
	