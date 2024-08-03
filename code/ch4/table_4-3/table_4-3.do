/*Table 4.3*/
capture clear all
capture log close
set more off
set maxvar 20000

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

log using "${logdir}table_4-3.log", replace 

set seed 3308004

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==. | faminc_adj_age3539==. | ///
	log_faminc_adj_age3539==.

egen std_cesd_age40=std(cesd_age40)

///simulation estimators

//OEhat^sim, IDEhat^sim, IIEhat^sim, CDEhat^sim

//define function for simulation estimators w/ exposure-mediator interaction
capture program drop medsim
program define medsim, rclass

	syntax [, nsim(integer 10)] 
	
	//fit logit model for L 
	logit ever_unemp_age3539 ///
		att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	est store Lmodel

	//fit normal linear model for M 
	reg log_faminc_adj_age3539 ///
		att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	est store Mmodel
	
	//fit normal linear model for Y 
	reg std_cesd_age40 ///
		c.log_faminc_adj_age3539##i.att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3 ever_unemp_age3539
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

//compute effect estimates
quietly medsim, nsim(2000)

di r(OE)
di r(IDE)
di r(IIE)
di r(CDE)

quietly medsimx, nsim(2000)

di r(OE)
di r(IDE)
di r(IIE)
di r(CDE)

log close


