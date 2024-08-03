/*Table 4.5*/
capture clear all
capture log close
set more off
set maxvar 40000

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

log using "${logdir}table_4-5.log", replace 

*ssc install rwrmed 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==. | faminc_adj_age3539==. | ///
	log_faminc_adj_age3539==.

egen std_cesd_age40=std(cesd_age40)

///RWR estimates

//compute bootstrap CIs
quietly bootstrap ///
	OE=_b[rATE] IDE=_b[rNDE] IIE=_b[rNIE] CDE=_b[CDE], ///
	saving("${datadir}\bootmed.dta", replace) ///
	reps(2000) seed(60637): ///
	rwrmed std_cesd_age40 ever_unemp_age3539, ///
		avar(att22) mvar(log_faminc_adj_age3539) cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
		a(0) astar(1) m(10.82) mreg(regress) cxa cxm lxm
		
mat list e(ci_percentile)

//compute p-values for null of no effect
use "${datadir}\bootmed.dta", clear

foreach x in OE IDE IIE CDE {
	quietly gen `x'_lt=0 
	quietly replace `x'_lt=1 if `x'<0
	quietly sum `x'_lt
	scalar p_lt = r(mean)
	
	quietly gen `x'_rt=0 
	quietly replace `x'_rt=1 if `x'>0
	quietly sum `x'_rt
	scalar p_rt = r(mean)
	
	di "`x' p-value"
	di min(p_rt, p_lt)*2
	di "-----------"
	
	drop `x'_*
	}

//reload data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==. | faminc_adj_age3539==. | ///
	log_faminc_adj_age3539==.

egen std_cesd_age40=std(cesd_age40)

///simulation estimates

//define function for simulation estimators w/ all two-way interactions
capture program drop medsimx
program define medsimx, rclass

	syntax [, nsim(integer 10)] 
	
	//fit logit model for L 
	logit ever_unemp_age3539 ///
		i.att22##i.female ///
		i.att22##i.black ///
		i.att22##i.hispan ///
		i.att22##c.paredu ///
		i.att22##i.parprof ///
		i.att22##c.parinc_prank ///
		i.att22##c.famsize ///
		i.att22##c.afqt3
		
	est store Lmodel

	//fit normal linear model for M 
	reg log_faminc_adj_age3539 ///
		i.att22##i.female ///
		i.att22##i.black ///
		i.att22##i.hispan ///
		i.att22##c.paredu ///
		i.att22##i.parprof ///
		i.att22##c.parinc_prank ///
		i.att22##c.famsize ///
		i.att22##c.afqt3
		
	est store Mmodel
	
	//fit normal linear model for Y 
	reg std_cesd_age40 ///
		i.att22##c.log_faminc_adj_age3539 ///
		i.att22##i.female ///
		i.att22##i.black ///
		i.att22##i.hispan ///
		i.att22##c.paredu ///
		i.att22##i.parprof ///
		i.att22##c.parinc_prank ///
		i.att22##c.famsize ///
		i.att22##c.afqt3 ///
		c.log_faminc_adj_age3539##i.female ///
		c.log_faminc_adj_age3539##i.black ///
		c.log_faminc_adj_age3539##i.hispan ///
		c.log_faminc_adj_age3539##c.paredu ///
		c.log_faminc_adj_age3539##i.parprof ///
		c.log_faminc_adj_age3539##c.parinc_prank ///
		c.log_faminc_adj_age3539##c.famsize ///
		c.log_faminc_adj_age3539##c.afqt3 ///
		c.log_faminc_adj_age3539##i.ever_unemp_age3539
		
	est store Ymodel
	
	gen D_orig=att22
	gen L_orig=ever_unemp_age3539
	gen M_orig=log_faminc_adj_age3539
	
	//simulate potential values of exposure-induced confounder, mediator, and outcome
	forval i=1/`nsim' {
		
		replace att22=0
		
		est restore Lmodel
		predict phat_L0, pr
		gen L0_`i'=rbinomial(1,phat_L0) //MC draws from Bern(p=phat_L0)
		
		est restore Mmodel
		predict yhat_M0, xb
		gen M0_`i'=rnormal(yhat_M0,e(rmse)) //MC draws from N(mu=yhat_M0,sd=e(rmse))
		
		replace att22=1
		
		est restore Lmodel
		predict phat_L1, pr
		gen L1_`i'=rbinomial(1,phat_L1) //MC draws from Bern(p=phat_L1)
		
		est restore Mmodel
		predict yhat_M1, xb
		gen M1_`i'=rnormal(yhat_M1,e(rmse)) //MC draws from N(mu=yhat_M1,sd=e(rmse))
		
		est restore Ymodel

		replace att22=0
		replace ever_unemp_age3539=L0_`i'
		replace log_faminc_adj_age3539=M0_`i'
		
		predict yhat_Y0L0M0, xb
		gen Y0L0M0_`i'=rnormal(yhat_Y0L0M0,e(rmse)) //MC draws from N(mu=yhat_Y0L0M0,sd=e(rmse))
		
		replace att22=1
		replace ever_unemp_age3539=L1_`i'
		replace log_faminc_adj_age3539=M1_`i'
		
		predict yhat_Y1L1M1, xb
		gen Y1L1M1_`i'=rnormal(yhat_Y1L1M1,e(rmse)) //MC draws from N(mu=yhat_Y1L1M1,sd=e(rmse))
	
		replace log_faminc_adj_age3539=M0_`i'

		predict yhat_Y1L1M0, xb
		gen Y1L1M0_`i'=rnormal(yhat_Y1L1M0,e(rmse)) //MC draws from N(mu=yhat_Y1L1M1,sd=e(rmse))
	
	    replace log_faminc_adj_age3539=ln(50000)
	
		predict yhat_Y1L1m, xb
		gen Y1L1m_`i'=rnormal(yhat_Y1L1m,e(rmse)) //MC draws from N(mu=yhat_Y1L1m,sd=e(rmse))
		
		replace att22=0
		replace ever_unemp_age3539=L0_`i'
		
		predict yhat_Y0L0m, xb
		gen Y0L0m_`i'=rnormal(yhat_Y0L0m,e(rmse)) //MC draws from N(mu=yhat_Y0L0m,sd=e(rmse))
		
		drop phat_* yhat_* M0_* M1_* L0_* L1_*
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
	replace att22=D_orig
	replace ever_unemp_age3539=L_orig
	replace log_faminc_adj_age3539=M_orig
	
	drop D_orig L_orig M_orig _est_Lmodel _est_Mmodel _est_Ymodel
	
end

//compute bootstrap CIs
quietly bootstrap ///
	OE=r(OE) IDE=r(IDE) IIE=r(IIE) CDE=r(CDE), ///
	saving("${datadir}\bootmed.dta", replace) ///
	reps(2000) seed(60637): ///
	medsimx, nsim(2000)
		
mat list e(ci_percentile)

//compute p-values for null of no effect
use "${datadir}\bootmed.dta", clear

foreach x in OE IDE IIE CDE {
	quietly gen `x'_lt=0 
	quietly replace `x'_lt=1 if `x'<0
	quietly sum `x'_lt
	scalar p_lt = r(mean)
	
	quietly gen `x'_rt=0 
	quietly replace `x'_rt=1 if `x'>0
	quietly sum `x'_rt
	scalar p_rt = r(mean)
	
	di "`x' p-value"
	di min(p_rt, p_lt)*2
	di "-----------"
	
	drop `x'_*
	}

//reload data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==. | faminc_adj_age3539==. | ///
	log_faminc_adj_age3539==.

egen std_cesd_age40=std(cesd_age40)

///ipw estimates

//define function for ipw estimator for OE, IDE, IIE, and CDE
capture program drop ipwmedx
program define ipwmedx, rclass 
	
	gen D_orig=att22
	gen L_orig=ever_unemp_age3539
	gen M_orig=log_faminc_adj_age3539
	
	//fit models
	logit att22 female black hispan paredu parprof parinc_prank famsize afqt3
	
	est store Dmodel

	logit ever_unemp_age3539 ///
	i.att22##i.female ///
	i.att22##i.black ///
	i.att22##i.hispan ///
	i.att22##c.paredu ///
	i.att22##i.parprof ///
	i.att22##c.parinc_prank ///
	i.att22##c.famsize ///
	i.att22##c.afqt3
	
	est store Lmodel
	
	reg log_faminc_adj_age3539 ///
	i.att22##i.female ///
	i.att22##i.black ///
	i.att22##i.hispan ///
	i.att22##c.paredu ///
	i.att22##i.parprof ///
	i.att22##c.parinc_prank ///
	i.att22##c.famsize ///
	i.att22##c.afqt3 ///
	i.att22##i.ever_unemp_age3539

	est store Mmodel

	//predict D probabilities
	est restore Dmodel
	
	predict phatD1_C, pr
	gen phatD0_C=1-phatD1_C
	
	//predict L probabilities
	est restore Lmodel
	
	replace att22=1
	
	predict phatL1_D1C, pr
	gen phatL0_D1C=1-phatL1_D1C
	
	replace att22=0

	predict phatL1_D0C, pr
	gen phatL0_D0C=1-phatL1_D0C

	//predict M probabilities
	est restore Mmodel
	
	replace att22=1
	
	predict mhat_D1LC, xb
	gen phatM_D1LC=normalden(log_faminc_adj_age3539, mhat_D1LC, e(rmse))
	
	replace att22=0
	
	predict mhat_D0LC, xb
	gen phatM_D0LC=normalden(log_faminc_adj_age3539, mhat_D0LC, e(rmse))
	
	replace ever_unemp_age3539=0
	
	predict mhat_D0L0C, xb
	gen phatM_D0L0C=normalden(log_faminc_adj_age3539, mhat_D0L0C, e(rmse))

	replace ever_unemp_age3539=1
	
	predict mhat_D0L1C, xb
	gen phatM_D0L1C=normalden(log_faminc_adj_age3539, mhat_D0L1C, e(rmse))
	
	replace att22=1
	
	predict mhat_D1L1C, xb
	gen phatM_D1L1C=normalden(log_faminc_adj_age3539, mhat_D1L1C, e(rmse))
	
	replace ever_unemp_age3539=0
	
	predict mhat_D1L0C, xb
	gen phatM_D1L0C=normalden(log_faminc_adj_age3539, mhat_D1L0C, e(rmse))

	replace att22=D_orig
	replace ever_unemp_age3539=L_orig
	replace log_faminc_adj_age3539=M_orig
	
	//compute stabilized weights
	logit att22
	predict phatD1, pr
	gen phatD0=1-phatD1
	
	reg log_faminc_adj_age3539 att22 
	
	replace att22=0
	
	predict mhat_D0, xb
	gen phatM_D0=normalden(log_faminc_adj_age3539, mhat_D0, e(rmse))

	replace att22=1
	
	predict mhat_D1, xb
	gen phatM_D1=normalden(log_faminc_adj_age3539, mhat_D1, e(rmse))
	
	replace att22=D_orig
	replace ever_unemp_age3539=L_orig
	replace log_faminc_adj_age3539=M_orig

	gen _sw1 = (phatD0/phatD0_C) * (1/phatM_D0LC) * ((phatM_D0L0C * phatL0_D0C) + (phatM_D0L1C * phatL1_D0C)) if att22==0 
	gen _sw2 = (phatD1/phatD1_C) * (1/phatM_D1LC) * ((phatM_D1L0C * phatL0_D1C) + (phatM_D1L1C * phatL1_D1C)) if att22==1
	gen _sw3 = (phatD1/phatD1_C) * (1/phatM_D1LC) * ((phatM_D0L0C * phatL0_D0C) + (phatM_D0L1C * phatL1_D0C)) if att22==1
    
	gen _sw4 = (phatD1/phatD1_C) * (phatM_D1/phatM_D1LC) if att22==1
	replace _sw4 = (phatD0/phatD0_C) * (phatM_D0/phatM_D0LC) if att22==0
	
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

	//compute interventional effect estimates
	return scalar OE = `Ehat_Y1M1' - `Ehat_Y0M0'
	return scalar IDE = `Ehat_Y1M0' - `Ehat_Y0M0'
	return scalar IIE = `Ehat_Y1M1' - `Ehat_Y1M0'
	
	//compute CDE estimate
	reg std_cesd_age40 i.att22##c.log_faminc_adj_age3539 [pw=_sw4]
	
	return scalar CDE = _b[1.att22] + ln(50000)*_b[1.att22#c.log_faminc_adj_age3539]
	
	drop phat* mhat* _sw* _est* D_orig L_orig M_orig
	
end

//compute bootstrap CIs
quietly bootstrap ///
	OE=r(OE) IDE=r(IDE) IIE=r(IIE) CDE=r(CDE), ///
	saving("${datadir}\bootmed.dta", replace) ///
	reps(2000) seed(60637): ///
	ipwmedx
		
mat list e(ci_percentile)

//compute p-values for null of no effect
use "${datadir}\bootmed.dta", clear

foreach x in OE IDE IIE CDE {
	quietly gen `x'_lt=0 
	quietly replace `x'_lt=1 if `x'<0
	quietly sum `x'_lt
	scalar p_lt = r(mean)
	
	quietly gen `x'_rt=0 
	quietly replace `x'_rt=1 if `x'>0
	quietly sum `x'_rt
	scalar p_rt = r(mean)
	
	di "`x' p-value"
	di min(p_rt, p_lt)*2
	di "-----------"
	
	drop `x'_*
	}

capture clear
erase "${datadir}\bootmed.dta"

log close
