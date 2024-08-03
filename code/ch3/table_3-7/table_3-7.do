/*Table 3.7*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

log using "${logdir}table_3-7.log", replace 

//input data
use "${datadir}JOBSII\Jobs-NoMiss-Binary.dta", clear

///linear model estimators
*ssc install rwrmed 

//ATEhat^lmi, NDEhat^lmi, NIEhat^lmi, CDE4hat^lmi

quietly bootstrap ATE=_b[ATE] NDE=_b[NDE] NIE=_b[NIE] CDE=_b[CDE], reps(2000) seed(60637): ///
	rwrmed work1, ///
		avar(treat) mvar(job_seek) cvar(econ_hard sex age nonwhite educ income) ///
		a(0) astar(1) m(4) mreg(regress)
		
mat list e(b)
mat list e(ci_percentile)

///simulation estimators

//define function for simulation estimator
capture program drop medsim
program define medsim, rclass

	syntax [, nsim(integer 10)] 
	
	//fit normal linear model for M 
	reg job_seek ///
		treat ///
		econ_hard sex age nonwhite educ income
	est store Mmodel
	
	//fit logit model for Y 
	logit work1 ///
		c.job_seek##i.treat ///
		econ_hard sex age nonwhite educ income
	est store Ymodel
	
	gen D_orig=treat
	gen M_orig=job_seek
	
	//simulate potential values of mediator and outcome for NDE, NIE, and ATE
	forval i=1/`nsim' {
		
		est restore Mmodel

		replace treat=0
		predict Ehat_M0, xb
		gen M0_`i'=rnormal(Ehat_M0,e(rmse)) //MC draws from N(mu=yhat_M0,sd=e(rmse))
		
		replace treat=1
		predict Ehat_M1, xb
		gen M1_`i'=rnormal(Ehat_M1,e(rmse)) //MC draws from N(mu=yhat_M1,sd=e(rmse))
		
		est restore Ymodel

		replace treat=0
		replace job_seek=M0_`i'
		predict phat_Y0M0, pr
		gen Y0M0_`i'=rbinomial(1,phat_Y0M0) //MC draws from Bern(p=phat_Y0M0)
		
		replace job_seek=M1_`i'
		predict phat_Y0M1, pr
		gen Y0M1_`i'=rbinomial(1,phat_Y0M1) //MC draws from Bern(p=phat_Y0M1)
	
		replace treat=1
		predict phat_Y1M1, pr
		gen Y1M1_`i'=rbinomial(1,phat_Y1M1) //MC draws from Bern(p=phat_Y1M1)
		
		replace job_seek=M0_`i'
		predict phat_Y1M0, pr
		gen Y1M0_`i'=rbinomial(1,phat_Y1M0) //MC draws from Bern(p=phat_Y1M0)
		
		drop phat_* Ehat_* M0_* M1_*
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
	
	replace treat=1
	replace job_seek=4
	predict phat_Y1m, pr

	replace treat=0
	predict phat_Y0m, pr
	
	//compute effect estimates
	gen _diff=phat_Y1m-phat_Y0m
	
	reg _diff
	return scalar CDE=_b[_cons]
	
	drop phat_Y1m phat_Y0m _diff

	//restore original data
	replace treat=D_orig
	replace job_seek=M_orig

	drop D_orig M_orig _est_Mmodel _est_Ymodel

end

quietly bootstrap ATE=r(ATE) NDE=r(NDE) NIE=r(NIE) CDE=r(CDE), reps(2000) seed(60637): ///
	medsim, nsim(1000)
	
mat list e(b)
mat list e(ci_percentile)

///inverse probability weighting (ipw) estimators

//ATEhat^ipw, NDEhat^ipw, NIEhat^ipw

//define function for ipw estimator for ATE, NDE, and NIE
capture program drop ipwmed
program define ipwmed, rclass

	//fit logit model for D given C
	logit treat ///
		econ_hard sex age nonwhite educ income
	est store DModel_1

	//fit logit model for D given C and M
	logit treat ///
		job_seek econ_hard sex age nonwhite educ income
	est store DModel_2

	//predict exposure probabilities
	est restore DModel_1
	predict phat_d_C, pr
	gen phat_dstar_C=1-phat_d_C

	est restore DModel_2
	predict phat_d_CM, pr
	gen phat_dstar_CM=1-phat_d_CM

	//compute stabilized weights
	logit treat
	predict phat_d, pr
	gen phat_dstar=1-phat_d
	
	gen _sw1=phat_dstar/phat_dstar_C if treat==0
	gen _sw2=phat_d/phat_d_C if treat==1
	gen _sw3=(phat_dstar_CM*phat_d)/(phat_d_CM*phat_dstar_C) if treat==1

	//compute weighted means of Y
	reg work1 [pw=_sw1] if treat==0
	local Ehat_Y0M0=_b[_cons]

	reg work1 [pw=_sw2] if treat==1
	local Ehat_Y1M1=_b[_cons]

	reg work1 [pw=_sw3] if treat==1
	local Ehat_Y1M0=_b[_cons]

	//compute effect estimates
	return scalar ATE=`Ehat_Y1M1'-`Ehat_Y0M0'
	return scalar NDE=`Ehat_Y1M0'-`Ehat_Y0M0'
	return scalar NIE=`Ehat_Y1M1'-`Ehat_Y1M0'
	
	drop phat_d_C phat_dstar_C phat_d_CM phat_dstar_CM phat_d phat_dstar _sw* ///
		_est_DModel_1 _est_DModel_2
	
	//reset estimation sample for bootstrap
	reg treat
	
end

quietly bootstrap ATE=r(ATE) NDE=r(NDE) NIE=r(NIE), reps(2000) seed(60637): ipwmed
	
mat list e(b)
mat list e(ci_percentile)

//CDEhat^ipw

//define function for ipw estimator for CDE
capture program drop ipwcde
program define ipwcde, rclass

	syntax [, mval(real 0)] 
	
	//fit logit model for D given C
	logit treat ///
		econ_hard sex age nonwhite educ income
	est store DModel
	
	//fit normal linear model for M given D and C
	reg job_seek ///
		 treat econ_hard sex age nonwhite educ income
	est store MModel

	//predict exposure probabilities
	est restore DModel
	predict phat_d_C, pr
	gen phat_D_denom=(treat*phat_d_C)+((1-treat)*(1-phat_d_C))

	est restore MModel
	predict Ehat_M_CD, xb
	gen phat_M_denom=normalden(job_seek,Ehat_M_CD,e(rmse))

	//compute stabilized weights
	logit treat
	predict phat_d, xb
	gen phat_D_num=(treat*phat_d)+((1-treat)*(1-phat_d))
	
	reg job_seek treat
	predict Ehat_M_D, xb
	gen phat_M_num=normalden(job_seek,Ehat_M_D,e(rmse))
	
	gen _sw4=(phat_M_num*phat_D_num)/(phat_M_denom*phat_D_denom)

	//censor stabilized weights at 1st and 99th percentiles
	centile _sw4, c(1 99) 
	replace _sw4=r(c_1) if _sw4<r(c_1) & _sw4!=.
	replace _sw4=r(c_2) if _sw4>r(c_2) & _sw4!=.

	//fit outcome model
	reg work1 i.treat##c.job_seek [pw=_sw4]

	//compute effect estimates
	return scalar CDE=_b[1.treat]+_b[1.treat#c.job_seek]*`mval'
	
	drop phat* Ehat_* _est_* _sw*
	
end

quietly bootstrap CDE=r(CDE), reps(2000) seed(60637): ipwcde, mval(4)
	
mat list e(b)
mat list e(ci_percentile)

log close
