/*Table 3.7*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/linmed, replace 
github install causalMedAnalysis/lincde, replace
github install causalMedAnalysis/medsim, replace 
github install causalMedAnalysis/impcde, replace 
github install causalMedAnalysis/ipwmed, replace 

//specify directories 
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

//download data
copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/JOBSII/Jobs-NoMiss-Binary.dta" ///
	"${datadir}JOBSII\"

//open log
log using "${logdir}table_3-7.log", replace 

//load data
use "${datadir}JOBSII\Jobs-NoMiss-Binary.dta", clear

//define macros for different variables
global C econ_hard sex age nonwhite educ income //baseline confounders
global D treat //exposure
global M job_seek //mediator
global Y work1 //outcome

//set seed to ensure reproducibility
set seed 3308004

//compute point and interval estimates based on linear models
linmed $Y $M, dvar($D) d(1) dstar(0) cvars($C) reps(2000) seed(60637)
lincde $Y, dvar($D) mvar($M) d(1) dstar(0) m(4) cvars($C) reps(2000) seed(60637)

//compute point and interval estimates based on simulation/imputation
medsim $Y, dvar($D) mvar($M) d(1) dstar(0) yreg(logit) mreg(regress) ///
	cvars($C) nsim(1000) reps(2000) seed(60637)

impcde $Y, dvar($D) mvar($M) d(1) dstar(0) m(4) yreg(logit) cvars($C) ///
	reps(2000) seed(60637) 

//compute point and interval estimates based on inverse probability weighting
ipwmed $Y $M, dvar($D) d(1) dstar(0) cvars($C) censor reps(2000) seed(60637)

//define function to estimate CDE using IPW with a continuous mediator
capture program drop ipwcde_Mcon
program define ipwcde_Mcon, rclass

	//fit logit model for D given C
	logit treat econ_hard sex age nonwhite educ income
	est store DModel
	
	//fit normal linear model for M given D and C
	reg job_seek treat econ_hard sex age nonwhite educ income
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
	return scalar cde=_b[1.treat]+_b[1.treat#c.job_seek]*4
	
	drop phat* Ehat_* _est_* _sw*
	
end

bootstrap CDE=r(cde), reps(2000) seed(60637) noheader notable: ipwcde_Mcon
estat bootstrap, p noheader

log close

//some of the estimates differ slightly from those reported in the text, 
//which are based on the R implementation. This is variously due to slight 
//differences in how the weights are censored, to Monte Carlo error, or to 
//differences in random number seeding that influence the bootstrap samples

//note: medsim with large nsim() and reps() can be time consuming to run in 
//Stata because it cannot parallelize the bootstrap replications. Consider 
//switching to R implementation of medsim if computation time is a concern.
