/*Table 6.1*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch6\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

log using "${logdir}table_6-1.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if missing(cesd_age40, att22, ever_unemp_age3539, female, black, ///
	hispan, paredu, parprof, parinc_prank, famsize, afqt3, log_faminc_adj_age3539)

egen std_cesd_age40=std(cesd_age40)

///regression imputation estimator
capture program drop regimp
program define regimp, rclass
	
	tempvar dvar_orig
	qui gen `dvar_orig' = att22
	
	reg std_cesd_age40 att22 female black hispan paredu parprof parinc_prank famsize afqt3
	
	replace att22 = 0
	tempvar yhat0
	predict `yhat0'
	
	replace att22 = 1
	tempvar yhat1
	predict `yhat1'
	
	reg `yhat0'
	return scalar psi0 = _b[_cons]

	reg `yhat1'
	return scalar psi1 = _b[_cons]
	
	replace att22 = `dvar_orig'
	
end regimp
	
qui bootstrap psi0=r(psi0) psi1=r(psi1) ATE=(r(psi1)-r(psi0)), ///
	reps(2000) seed(60637): regimp

estat bootstrap, p noheader

///inverse probability weighting estimator
capture program drop ipwate
program define ipwate, rclass
	
	logit att22 female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar phat_D1_C
	predict `phat_D1_C', pr
	
	tempvar ipwt
	gen `ipwt' = (att22/`phat_D1_C') + ((1-att22)/(1-`phat_D1_C'))
	
	forval d=0/1 {
		qui centile `ipwt' if `ipwt'!=. & att22==`d', c(1 99) 
		qui replace `ipwt'=r(c_1) if `ipwt'<r(c_1) & `ipwt'!=. & att22==`d'
		qui replace `ipwt'=r(c_2) if `ipwt'>r(c_2) & `ipwt'!=. & att22==`d'
		}
	
	reg std_cesd_age40 att22 [pw=`ipwt']
	return scalar psi0 = _b[_cons]
	return scalar psi1 = _b[_cons] + _b[att22]
	
end ipwate
	
qui bootstrap psi0=r(psi0) psi1=r(psi1) ATE=(r(psi1)-r(psi0)), ///
	reps(2000) seed(60637): ipwate

estat bootstrap, p noheader

///parametric DR estimator
capture program drop drate
program define drate, rclass
	
	logit att22 female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar phat_D1_C
	predict `phat_D1_C', pr
	
	tempvar ipwt
	gen `ipwt' = (att22/`phat_D1_C') + ((1-att22)/(1-`phat_D1_C'))
	
	forval d=0/1 {
		qui centile `ipwt' if `ipwt'!=. & att22==`d', c(1 99) 
		qui replace `ipwt'=r(c_1) if `ipwt'<r(c_1) & `ipwt'!=. & att22==`d'
		qui replace `ipwt'=r(c_2) if `ipwt'>r(c_2) & `ipwt'!=. & att22==`d'
		}
	
	tempvar dvar_orig
	qui gen `dvar_orig' = att22
	
	reg std_cesd_age40 att22 female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar yres
	predict `yres', resid
	
	replace att22 = 0
	tempvar yhat0
	predict `yhat0'
	
	replace att22 = 1
	tempvar yhat1
	predict `yhat1'
	
	replace att22 = `dvar_orig'
	
	tempvar dr0_summand
	qui gen `dr0_summand' = `yhat0' + `ipwt'*`yres'*(1-att2)
	
	reg `dr0_summand'
	return scalar psi0 = _b[_cons]
	
	tempvar dr1_summand
	qui gen `dr1_summand' = `yhat1' + `ipwt'*`yres'*(att2)
	
	reg `dr1_summand'
	return scalar psi1 = _b[_cons]
	
end drate

qui bootstrap psi0=r(psi0) psi1=r(psi1) ATE=(r(psi1)-r(psi0)), ///
	reps(2000) seed(60637): drate

estat bootstrap, p noheader

///DML estimator

//note that stata does not currently support use of a superLearner
//we therefore implement the DML estimator using only random forests

*ssc install rforest

capture program drop dmlate
program define dmlate, rclass

	tempvar dvar_orig
	qui gen `dvar_orig' = att22
	
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

	forval k=1/5 {

		qui rforest att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60657)

		tempvar	xxphatD0_C xxphatD1_C
		qui predict `xxphatD0_C' `xxphatD1_C' if `kpart'==`k', pr
		qui replace `phat_D1_C' = `xxphatD1_C' if `kpart'==`k'

		qui rforest std_cesd_age40 att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60657)

		tempvar	xxyhat
		qui predict `xxyhat' if `kpart'==`k'
		qui replace `yres' = std_cesd_age40 - `xxyhat' if `kpart'==`k'

		qui replace att22 = 1
		tempvar	xxyhat1
		qui predict `xxyhat1' if `kpart'==`k'
		qui replace `yhat1' = `xxyhat1' if `kpart'==`k'

		qui replace att22 = 0 
		tempvar xxyhat0
		qui predict `xxyhat0' if `kpart'==`k'
		qui replace `yhat0' = `xxyhat0' if `kpart'==`k'

		qui replace att22 = `dvar_orig'
		}

	tempvar ipwt
	qui gen `ipwt' = (att22/`phat_D1_C') + ((1-att22)/(1-`phat_D1_C'))
	
	forval d=0/1 {
		qui centile `ipwt' if `ipwt'!=. & att22==`d', c(1 99) 
		qui replace `ipwt'=r(c_1) if `ipwt'<r(c_1) & `ipwt'!=. & att22==`d'
		qui replace `ipwt'=r(c_2) if `ipwt'>r(c_2) & `ipwt'!=. & att22==`d'
		}
	
	tempvar dml0_summand
	qui gen `dml0_summand' = `yhat0' + `ipwt'*`yres'*(1-att2)
	
	qui reg `dml0_summand'
	return scalar psi0 = _b[_cons]
	return scalar se_psi0 = sqrt((e(rmse)^2)/_N)
	
	tempvar dml1_summand
	qui gen `dml1_summand' = `yhat1' + `ipwt'*`yres'*(att2)
	
	qui reg `dml1_summand'
	return scalar psi1 = _b[_cons]
	return scalar se_psi1 = sqrt((e(rmse)^2)/_N)
	
	tempvar dml1v0_summand
	qui gen `dml1v0_summand' = `dml1_summand' - `dml0_summand'
	
	reg `dml1v0_summand'
	return scalar ate = _b[_cons]
	return scalar se_ate = sqrt((e(rmse)^2)/_N)
	
end dmlate

dmlate

di r(psi0) // point estimate
di r(psi0)-1.96*r(se_psi0) // lower limit 95% CI
di r(psi0)+1.96*r(se_psi0) // upper limit 95% CI

di r(psi1)
di r(psi1)-1.96*r(se_psi1)
di r(psi1)+1.96*r(se_psi1)

di r(ate)
di r(ate)-1.96*r(se_ate)
di r(ate)+1.96*r(se_ate)

log close
