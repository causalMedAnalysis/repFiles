/*Table 6.1*/
capture clear all
capture log close
set more off

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch6\_LOGS\"

//download data
copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_6-1.log", replace 

//load data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, log_faminc_adj_age3539, ///
	female, black, hispan, paredu, parprof, parinc_prank, famsize, afqt3)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//define macros for different variables
global C female black hispan paredu parprof parinc_prank famsize afqt3 //baseline confounders
global D att22 //exposure
global Y std_cesd_age40 //outcome

//define regression imputation estimator for ATE
capture program drop regimp
program define regimp, rclass
	
	tempvar dvar_orig
	qui gen `dvar_orig' = $D
	
	//fit outcome regression
	reg $Y $D $C
	
	//impute potential outcomes
	replace $D = 0
	tempvar yhat0
	predict `yhat0'
	
	replace $D = 1
	tempvar yhat1
	predict `yhat1'
	
	reg `yhat0'
	return scalar psi0 = _b[_cons]

	reg `yhat1'
	return scalar psi1 = _b[_cons]
	
	replace $D = `dvar_orig'
	
end regimp

//compute regression imputation estimates for ATE
qui bootstrap ///
	psi0=r(psi0) ///
	psi1=r(psi1) ///
	ATE=(r(psi1)-r(psi0)), ///
		reps(2000) seed(60637): regimp

estat bootstrap, p noheader

//define inverse probability weighting estimator for ATE
capture program drop ipwate
program define ipwate, rclass
	
	//fit exposure model
	logit $D $C
	
	//compute IPWs
	tempvar phat_D1_C
	predict `phat_D1_C', pr
	
	tempvar ipwt
	gen `ipwt' = ($D/`phat_D1_C') + ((1-$D)/(1-`phat_D1_C'))
	
	//censor IPWs at 1st and 99th percentiles
	forval d=0/1 {
		qui centile `ipwt' if `ipwt'!=. & $D==`d', c(1 99) 
		qui replace `ipwt'=r(c_1) if `ipwt'<r(c_1) & `ipwt'!=. & $D ==`d'
		qui replace `ipwt'=r(c_2) if `ipwt'>r(c_2) & `ipwt'!=. & $D ==`d'
		}
	
	//compute weighted outcome means
	reg $Y $D [pw=`ipwt']
	return scalar psi0 = _b[_cons]
	return scalar psi1 = _b[_cons] + _b[$D]
	
end ipwate

//compute IPW estimates for ATE	
qui bootstrap ///
	psi0=r(psi0) ///
	psi1=r(psi1) ///
	ATE=(r(psi1)-r(psi0)), ///
		reps(2000) seed(60637): ipwate

estat bootstrap, p noheader

//define parametric doubly-robust estimator for ATE
capture program drop drate
program define drate, rclass
	
	//fit exposure model
	logit $D $C
	
	//compute IPWs
	tempvar phat_D1_C
	predict `phat_D1_C', pr
	
	tempvar ipwt
	gen `ipwt' = ($D/`phat_D1_C') + ((1-$D)/(1-`phat_D1_C'))
	
	//censor IPWs at 1st and 99th percentiles
	forval d=0/1 {
		qui centile `ipwt' if `ipwt'!=. & $D==`d', c(1 99) 
		qui replace `ipwt'=r(c_1) if `ipwt'<r(c_1) & `ipwt'!=. & $D ==`d'
		qui replace `ipwt'=r(c_2) if `ipwt'>r(c_2) & `ipwt'!=. & $D ==`d'
	}
	
	tempvar dvar_orig
	qui gen `dvar_orig' = $D
	
	//fit outcome regression
	reg $Y $D $C
	
	//compute residuals
	tempvar yres
	predict `yres', resid
	
	//impute outcomes
	replace $D = 0
	tempvar yhat0
	predict `yhat0'
	
	replace $D = 1
	tempvar yhat1
	predict `yhat1'
	
	replace $D = `dvar_orig'
	
	//compute eif summand
	tempvar dr0_summand
	qui gen `dr0_summand' = `yhat0' + `ipwt'*`yres'*(1-$D)
	
	//compute POMs
	reg `dr0_summand'
	return scalar psi0 = _b[_cons]
	
	tempvar dr1_summand
	qui gen `dr1_summand' = `yhat1' + `ipwt'*`yres'*($D)
	
	reg `dr1_summand'
	return scalar psi1 = _b[_cons]
	
end drate

//compute doubly robust estimates for ATE
qui bootstrap ///
	psi0=r(psi0) ///
	psi1=r(psi1) ///
	ATE=(r(psi1)-r(psi0)), ///
		reps(2000) seed(60637): drate

estat bootstrap, p noheader

//define DML estimator for ATE

//note that stata does not currently support use of a superLearner
//we therefore implement the DML estimator using only random forests

ssc install rforest

capture program drop dmlate
program define dmlate, rclass

	tempvar dvar_orig
	qui gen `dvar_orig' = $D
	
	//randomly partition the sample for cross-fitting
	tempvar kpart
	set seed 60637
	qui gen u = uniform()
	qui sort u
	qui gen `kpart' = ceil(_n/(_N/5))
	drop u

	tempvar yres yhat1 yhat0 phat_D1_C
	qui gen `yres' = .
	qui gen `yhat1' = .
	qui gen `yhat0' = .
	qui gen `phat_D1_C' = . 

	//initiate cross-fitting
	forval k=1/5 {

		//fit random forest for the exposure
		qui rforest $D $C if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60657)

		tempvar	xxphatD0_C xxphatD1_C
		qui predict `xxphatD0_C' `xxphatD1_C' if `kpart'==`k', pr
		qui replace `phat_D1_C' = `xxphatD1_C' if `kpart'==`k'

		//fit random forest for the outcome
		qui rforest $Y $D $C if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60657)

		//compute residuals	
		tempvar	xxyhat
		qui predict `xxyhat' if `kpart'==`k'
		qui replace `yres' = std_cesd_age40 - `xxyhat' if `kpart'==`k'

		//impute outcomes
		qui replace $D = 1
		tempvar	xxyhat1
		qui predict `xxyhat1' if `kpart'==`k'
		qui replace `yhat1' = `xxyhat1' if `kpart'==`k'

		qui replace $D = 0 
		tempvar xxyhat0
		qui predict `xxyhat0' if `kpart'==`k'
		qui replace `yhat0' = `xxyhat0' if `kpart'==`k'

		qui replace $D = `dvar_orig'
	}

	//compute IPWs
	tempvar ipwt
	qui gen `ipwt' = ($D/`phat_D1_C') + ((1-$D)/(1-`phat_D1_C'))
	
	//censor the IPWs at 1st and 99th percentiles
	forval d=0/1 {
		qui centile `ipwt' if `ipwt'!=. & $D == `d', c(1 99) 
		qui replace `ipwt'=r(c_1) if `ipwt'<r(c_1) & `ipwt'!=. & $D == `d'
		qui replace `ipwt'=r(c_2) if `ipwt'>r(c_2) & `ipwt'!=. & $D == `d'
	}

	//compute eif summands
	tempvar dml0_summand
	qui gen `dml0_summand' = `yhat0' + `ipwt'*`yres'*(1-$D)

	tempvar dml1_summand
	qui gen `dml1_summand' = `yhat1' + `ipwt'*`yres'*($D)
	
	tempvar dml1v0_summand
	qui gen `dml1v0_summand' = `dml1_summand' - `dml0_summand'
	
	//compute POMs and ATE
	qui mean `dml0_summand' `dml1_summand' `dml1v0_summand'
	
	scalar psi0 = _b[`dml0_summand']
	scalar psi1 = _b[`dml1_summand']
	scalar ate = _b[`dml1v0_summand']
	
	//compute normal-based CIs
	scalar ll95_psi0 = _b[`dml0_summand']-1.96*_se[`dml0_summand']
	scalar ul95_psi0 = _b[`dml0_summand']+1.96*_se[`dml0_summand']

	scalar ll95_psi1 = _b[`dml1_summand']-1.96*_se[`dml1_summand']
	scalar ul95_psi1 = _b[`dml1_summand']+1.96*_se[`dml1_summand']	
	
	scalar ll95_ate = _b[`dml1v0_summand']-1.96*_se[`dml1v0_summand']
	scalar ul95_ate = _b[`dml1v0_summand']+1.96*_se[`dml1v0_summand']	
	
	//print results
	matrix results = ///
		(psi0, ll95_psi0, ul95_psi0 \ ///
		psi1, ll95_psi1, ul95_psi1 \ ///
		ate, ll95_ate, ul95_ate)
	
	matrix colnames results = "Est." "[95% Conf." "Interval]"
	
	matrix rownames results = "psi0" "psi1" "ATE"

	matrix list results
	
end dmlate

//compute DML estimates for ATE
dmlate

log close

//note the estimates differ slightly from those reported in
//the text, which are based on the R implementation. This is due to minor 
//differences in how the weights are censored and to the use of random 
//forests exclusively rather than a super learner in the DML estimator.
