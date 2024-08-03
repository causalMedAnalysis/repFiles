/*Table 3.3*/
capture clear all
capture log close
set more off
set maxvar 20000

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

log using "${logdir}table_3-3.log", replace 

set seed 3308004

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

///simulation and imputation estimators

//ATEhat^sim, NDEhat^sim, NIEhat^sim, CDEhat^imp

//define function for simulation/imputation estimators w/o covariate interactions
capture program drop medsim
program define medsim, rclass

	syntax [, nsim(integer 10)] 
	
	//fit logit model for M 
	logit ever_unemp_age3539 ///
		att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	est store Mmodel
	
	//fit normal linear model for Y 
	reg std_cesd_age40 ///
		i.ever_unemp_age3539 i.att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	est store Ymodel
	
	gen D_orig=att22
	gen M_orig=ever_unemp_age3539
	
	//simulate potential values of mediator and outcome for NDE, NIE, and ATE
	forval i=1/`nsim' {
		
		est restore Mmodel

		replace att22=0
		predict phat_M0, pr
		gen M0_`i'=rbinomial(1,phat_M0) //MC draws from Bern(p=phat_M0)
		
		replace att22=1
		predict phat_M1, pr
		gen M1_`i'=rbinomial(1,phat_M1) //MC draws from Bern(p=phat_M1)
		
		est restore Ymodel

		replace att22=0
		replace ever_unemp_age3539=M0_`i'
		predict yhat_Y0M0, xb
		gen Y0M0_`i'=rnormal(yhat_Y0M0,e(rmse)) //MC draws from N(mu=yhat_Y0M0,sd=e(rmse))
		
		replace ever_unemp_age3539=M1_`i'
		predict yhat_Y0M1, xb
		gen Y0M1_`i'=rnormal(yhat_Y0M1,e(rmse)) //MC draws from N(mu=yhat_Y0M1,sd=e(rmse))
	
		replace att22=1
		predict yhat_Y1M1, xb
		gen Y1M1_`i'=rnormal(yhat_Y1M1,e(rmse)) //MC draws from N(mu=yhat_Y1M1,sd=e(rmse))

		replace ever_unemp_age3539=M0_`i'
		predict yhat_Y1M0, xb
		gen Y1M0_`i'=rnormal(yhat_Y1M0,e(rmse)) //MC draws from N(mu=yhat_Y1M0,sd=e(rmse))
		
		drop phat_* yhat_* M0_* M1_*
		}

	//average simulated outcomes over simulations
	egen Y0M0_=rowmean(Y0M0_*)
	egen Y0M1_=rowmean(Y0M1_*)
	egen Y1M1_=rowmean(Y1M1_*)
	egen Y1M0_=rowmean(Y1M0_*)
	
	//average simulated outcomes over sample members
	reg Y0M0_
	local Ehat_Y0M0=_b[_cons]
	drop Y0M0*

	reg Y0M1_
	local Ehat_Y0M1=_b[_cons]
	drop Y0M1*

	reg Y1M1_
	local Ehat_Y1M1=_b[_cons]
	drop Y1M1*

	reg Y1M0_
	local Ehat_Y1M0=_b[_cons]
	drop Y1M0*

	//compute effect estimates
	return scalar ATE=`Ehat_Y1M1'-`Ehat_Y0M0'
	return scalar NDE=`Ehat_Y1M0'-`Ehat_Y0M0'
	return scalar NIE=`Ehat_Y1M1'-`Ehat_Y1M0'
	
	//impute potential outcomes for CDE
	est restore Ymodel
	
	replace att22=1
	replace ever_unemp_age3539=0
	predict Yhat10, xb

	replace att22=0
	predict Yhat00, xb
	
	//compute effect estimates
	gen _diff=Yhat10-Yhat00
	
	reg _diff
	return scalar CDE=_b[_cons]
	
	drop Yhat10 Yhat00 _diff

	//restore original data
	replace att22=D_orig
	replace ever_unemp_age3539=M_orig
	
	drop D_orig M_orig _est_Mmodel _est_Ymodel
	
end

quietly medsim, nsim(2000)

di r(ATE)
di r(NDE)
di r(NIE)
di r(CDE)

//define function for simulation/imputation estimators w/ covariate interactions
capture program drop medsim
program define medsim, rclass

	syntax [, nsim(integer 10)] 
	
	//fit logit model for M 
	logit ever_unemp_age3539 ///
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
		i.att22##i.female ///
		i.att22##i.black ///
		i.att22##i.hispan ///
		i.att22##c.paredu ///
		i.att22##i.parprof ///
		i.att22##c.parinc_prank ///
		i.att22##c.famsize ///
		i.att22##c.afqt3 ///
		i.ever_unemp_age3539##i.female ///
		i.ever_unemp_age3539##i.black ///
		i.ever_unemp_age3539##i.hispan ///
		i.ever_unemp_age3539##c.paredu ///
		i.ever_unemp_age3539##i.parprof ///
		i.ever_unemp_age3539##c.parinc_prank ///
		i.ever_unemp_age3539##c.famsize ///
		i.ever_unemp_age3539##c.afqt3 ///
		i.ever_unemp_age3539##i.att22
		
	est store Ymodel
	
	gen D_orig=att22
	gen M_orig=ever_unemp_age3539
	
	//simulate potential values of mediator and outcome for NDE, NIE, and ATE
	forval i=1/`nsim' {
		
		est restore Mmodel

		replace att22=0
		predict phat_M0, pr
		gen M0_`i'=rbinomial(1,phat_M0) //MC draws from Bern(p=phat_M0)
		
		replace att22=1
		predict phat_M1, pr
		gen M1_`i'=rbinomial(1,phat_M1) //MC draws from Bern(p=phat_M1)
		
		est restore Ymodel

		replace att22=0
		replace ever_unemp_age3539=M0_`i'
		predict yhat_Y0M0, xb
		gen Y0M0_`i'=rnormal(yhat_Y0M0,e(rmse)) //MC draws from N(mu=yhat_Y0M0,sd=e(rmse))
		
		replace ever_unemp_age3539=M1_`i'
		predict yhat_Y0M1, xb
		gen Y0M1_`i'=rnormal(yhat_Y0M1,e(rmse)) //MC draws from N(mu=yhat_Y0M1,sd=e(rmse))
	
		replace att22=1
		predict yhat_Y1M1, xb
		gen Y1M1_`i'=rnormal(yhat_Y1M1,e(rmse)) //MC draws from N(mu=yhat_Y1M1,sd=e(rmse))

		replace ever_unemp_age3539=M0_`i'
		predict yhat_Y1M0, xb
		gen Y1M0_`i'=rnormal(yhat_Y1M0,e(rmse)) //MC draws from N(mu=yhat_Y1M0,sd=e(rmse))
		
		drop phat_* yhat_* M0_* M1_*
		}

	//average simulated outcomes over simulations
	egen Y0M0_=rowmean(Y0M0_*)
	egen Y0M1_=rowmean(Y0M1_*)
	egen Y1M1_=rowmean(Y1M1_*)
	egen Y1M0_=rowmean(Y1M0_*)
	
	//average simulated outcomes over sample members
	reg Y0M0_
	local Ehat_Y0M0=_b[_cons]
	drop Y0M0*

	reg Y0M1_
	local Ehat_Y0M1=_b[_cons]
	drop Y0M1*

	reg Y1M1_
	local Ehat_Y1M1=_b[_cons]
	drop Y1M1*

	reg Y1M0_
	local Ehat_Y1M0=_b[_cons]
	drop Y1M0*

	//compute effect estimates
	return scalar ATE=`Ehat_Y1M1'-`Ehat_Y0M0'
	return scalar NDE=`Ehat_Y1M0'-`Ehat_Y0M0'
	return scalar NIE=`Ehat_Y1M1'-`Ehat_Y1M0'
	
	//impute potential outcomes for CDE
	est restore Ymodel
	
	replace att22=1
	replace ever_unemp_age3539=0
	predict Yhat10, xb

	replace att22=0
	predict Yhat00, xb
	
	//compute effect estimates
	gen _diff=Yhat10-Yhat00
	
	reg _diff
	return scalar CDE=_b[_cons]
	
	drop Yhat10 Yhat00 _diff

	//restore original data
	replace att22=D_orig
	replace ever_unemp_age3539=M_orig
	
	drop D_orig M_orig _est_Mmodel _est_Ymodel
	
end

quietly medsim, nsim(2000)

di r(ATE)
di r(NDE)
di r(NIE)
di r(CDE)

log close


