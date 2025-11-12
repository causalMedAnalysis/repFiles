/*Table 6.3*/
capture clear all
capture log close
set more off

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_6-3.log", replace 

//load data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, log_faminc_adj_age3539, ///
	female, black, hispan, paredu, parprof, parinc_prank, famsize, afqt3)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//define multiply robust estimator for interventional effects
capture program drop mrvent
program define mrvent, rclass

	tempvar dvar_orig lvar_orig
	qui gen `dvar_orig' = att22
	qui gen `lvar_orig' = ever_unemp_age3539
	
	logit att22 female black hispan paredu parprof parinc_prank famsize afqt3
		
	tempvar phat_D0_C phat_D1_C phat_D_C
	predict `phat_D1_C', pr
	gen double `phat_D0_C' = 1 - `phat_D1_C'
	gen double `phat_D_C' = att22*`phat_D1_C' + (1-att22)*`phat_D0_C'
		
	logit att22 log_faminc_adj_age3539 female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar phat_D0_CM phat_D1_CM phat_D_CM
	predict `phat_D1_CM', pr
	gen double `phat_D0_CM' = 1 - `phat_D1_CM'
	gen double `phat_D_CM' = att22*`phat_D1_CM' + (1-att22)*`phat_D0_CM'

	logit ever_unemp_age3539 att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar phat_L0_CD phat_L1_CD phat_L_CD
	predict `phat_L1_CD', pr
	gen double `phat_L0_CD' = 1 - `phat_L1_CD'
	gen double `phat_L_CD' = ever_unemp_age3539*`phat_L1_CD' + (1-ever_unemp_age3539)*`phat_L0_CD'
	
	replace att22 = 1
	tempvar phat_L0_CD1 phat_L1_CD1 phat_L_CD1
	predict `phat_L1_CD1', pr
	gen double `phat_L0_CD1' = 1 - `phat_L1_CD1'	
	gen double `phat_L_CD1' = ever_unemp_age3539*`phat_L1_CD1' + (1-ever_unemp_age3539)*`phat_L0_CD1'

	replace att22 = 0
	tempvar phat_L0_CD0 phat_L1_CD0 phat_L_CD0
	predict `phat_L1_CD0', pr
	gen double `phat_L0_CD0' = 1 - `phat_L1_CD0'
	gen double `phat_L_CD0' = ever_unemp_age3539*`phat_L1_CD0' + (1-ever_unemp_age3539)*`phat_L0_CD0'
	
	replace att22 = `dvar_orig'
	
	logit ever_unemp_age3539 log_faminc_adj_age3539 att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar phat_L0_CDM phat_L1_CDM phat_L_CDM
	predict `phat_L1_CDM', pr
	gen double `phat_L0_CDM' = 1 - `phat_L1_CDM'
	gen double `phat_L_CDM' = ever_unemp_age3539*`phat_L1_CDM' + (1-ever_unemp_age3539)*`phat_L0_CDM'
	
	replace att22 = 1
	tempvar phat_L0_CD1M phat_L1_CD1M phat_L_CD1M
	predict `phat_L1_CD1M', pr
	gen double `phat_L0_CD1M' = 1 - `phat_L1_CD1M'	
	gen double `phat_L_CD1M' = ever_unemp_age3539*`phat_L1_CD1M' + (1-ever_unemp_age3539)*`phat_L0_CD1M'

	replace att22 = 0
	tempvar phat_L0_CD0M phat_L1_CD0M phat_L_CD0M
	predict `phat_L1_CD0M', pr
	gen double `phat_L0_CD0M' = 1 - `phat_L1_CD0M'
	gen double `phat_L_CD0M' = ever_unemp_age3539*`phat_L1_CD0M' + (1-ever_unemp_age3539)*`phat_L0_CD0M'
	
	replace att22 = `dvar_orig'	
			
	reg std_cesd_age40 i.att22##c.log_faminc_adj_age3539 ever_unemp_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar mu_CDLM
	predict `mu_CDLM'
	
	replace att22 = 1
	tempvar mu1_CML
	predict `mu1_CML'
	
	replace att22 = 0
	tempvar mu0_CML
	predict `mu0_CML'
	
	replace att22 = 0
	replace ever_unemp_age3539 = 0
	tempvar mu0_CML0
	predict `mu0_CML0'
	
	replace att22 = 0
	replace ever_unemp_age3539 = 1
	tempvar mu0_CML1
	predict `mu0_CML1'

	replace att22 = 1
	replace ever_unemp_age3539 = 0
	tempvar mu1_CML0
	predict `mu1_CML0'
	
	replace att22 = 1
	replace ever_unemp_age3539 = 1
	tempvar mu1_CML1
	predict `mu1_CML1'
	
	replace att22 = `dvar_orig'
	replace ever_unemp_age3539 = `lvar_orig'
	
	tempvar w_00 w_11 w_01 
	gen double `w_11' = `phat_L_CD1' / `phat_L_CD1M'
	gen double `w_00' = `phat_L_CD0' / `phat_L_CD0M'
	gen double `w_01' = (`phat_D1_C'*`phat_L_CD1'*`phat_D0_CM') / (`phat_D0_C'*`phat_L_CD1M'*`phat_D1_CM')
	
	foreach w in `w_11' `w_00' `w_01' {
		qui centile `w' if `w'!=., c(1 99) 
		qui replace `w'=r(c_1) if `w'<r(c_1) & `w'!=.
		qui replace `w'=r(c_2) if `w'>r(c_2) & `w'!=.
	}
	
	tempvar u_00_regressand u_11_regressand u_01_regressand
	
	gen double `u_00_regressand' = `mu0_CML'*`w_00'
	gen double `u_11_regressand' = `mu1_CML'*`w_11'
	gen double `u_01_regressand' = `mu1_CML'*`w_01'
	
	reg `u_00_regressand' i.att22 i.ever_unemp_age3539 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3
		
	replace att22 = 0
	tempvar u_00_CL
	predict `u_00_CL'
	
	replace att22 = 0
	replace ever_unemp_age3539 = 0
	tempvar u_00_CL0
	predict `u_00_CL0'
	
	replace att22 = 0
	replace ever_unemp_age3539 = 1
	tempvar u_00_CL1
	predict `u_00_CL1'
	
	replace att22 = `dvar_orig'
	replace ever_unemp_age3539 = `lvar_orig'
	
	reg `u_11_regressand' i.att22 i.ever_unemp_age3539 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3
		
	replace att22 = 1
	tempvar u_11_CL
	predict `u_11_CL'
	
	replace att22 = 1
	replace ever_unemp_age3539 = 0
	tempvar u_11_CL0
	predict `u_11_CL0'
	
	replace att22 = 1
	replace ever_unemp_age3539 = 1
	tempvar u_11_CL1
	predict `u_11_CL1'
		
	replace att22 = `dvar_orig'
	replace ever_unemp_age3539 = `lvar_orig'

	reg `u_01_regressand' i.att22 i.ever_unemp_age3539 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3 
		
	replace att22 = 1
	tempvar u_01_CL
	predict `u_01_CL'
	
	replace att22 = 1
	replace ever_unemp_age3539 = 0
	tempvar u_01_CL0
	predict `u_01_CL0'
	
	replace att22 = 1
	replace ever_unemp_age3539 = 1
	tempvar u_01_CL1
	predict `u_01_CL1'
		
	replace att22 = `dvar_orig'
	replace ever_unemp_age3539 = `lvar_orig'
	
	tempvar v_CD1_regressand v_CD0_regressand
	gen double `v_CD1_regressand' = `mu1_CML0'*`phat_L0_CD1' + `mu1_CML1'*`phat_L1_CD1'
	gen double `v_CD0_regressand' = `mu0_CML0'*`phat_L0_CD0' + `mu0_CML1'*`phat_L1_CD0'
	
	reg `v_CD1_regressand' i.att22 i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3
	
	tempvar v_CD11 v_CD01
	
	replace att22 = 1
	predict `v_CD11'
	
	replace att22 = 0
	predict `v_CD01'
	
	replace att22 = `dvar_orig'

	reg `v_CD0_regressand' i.att22 i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3
	
	tempvar v_CD00
	
	replace att22 = 0
	predict `v_CD00'
	
	replace att22 = `dvar_orig'
	
	tempvar ipw_D1 ipw_D0 
	gen double `ipw_D1' = att22/`phat_D1_C'
	gen double `ipw_D0' = (1-att22)/`phat_D0_C'
	
	tempvar ipw_X_w_11 ipw_X_w_00 ipw_X_w_01
	gen double `ipw_X_w_11' = `ipw_D1'*`w_11'
	gen double `ipw_X_w_00' = `ipw_D0'*`w_00'
	gen double `ipw_X_w_01' = `ipw_D1'*`w_01'
	
	centile `ipw_D0' if `ipw_D0'!=. & att22==0, c(1 99) 
	replace `ipw_D0'=r(c_1) if `ipw_D0'<r(c_1) & `ipw_D0'!=. & att22==0
	replace `ipw_D0'=r(c_2) if `ipw_D0'>r(c_2) & `ipw_D0'!=. & att22==0
	
	centile `ipw_D1' if `ipw_D1'!=. & att22==1, c(1 99) 
	replace `ipw_D1'=r(c_1) if `ipw_D1'<r(c_1) & `ipw_D1'!=. & att22==1
	replace `ipw_D1'=r(c_2) if `ipw_D1'>r(c_2) & `ipw_D1'!=. & att22==1
	
	centile `ipw_X_w_00' if `ipw_X_w_00'!=. & att22==0, c(1 99) 
	replace `ipw_X_w_00'=r(c_1) if `ipw_X_w_00'<r(c_1) & `ipw_X_w_00'!=. & att22==0
	replace `ipw_X_w_00'=r(c_2) if `ipw_X_w_00'>r(c_2) & `ipw_X_w_00'!=. & att22==0
		
	centile `ipw_X_w_11' if `ipw_X_w_11'!=. & att22==1, c(1 99) 
	replace `ipw_X_w_11'=r(c_1) if `ipw_X_w_11'<r(c_1) & `ipw_X_w_11'!=. & att22==1
	replace `ipw_X_w_11'=r(c_2) if `ipw_X_w_11'>r(c_2) & `ipw_X_w_11'!=. & att22==1
		
	centile `ipw_X_w_01' if `ipw_X_w_01'!=. & att22==1, c(1 99) 
	replace `ipw_X_w_01'=r(c_1) if `ipw_X_w_01'<r(c_1) & `ipw_X_w_01'!=. & att22==1
	replace `ipw_X_w_01'=r(c_2) if `ipw_X_w_01'>r(c_2) & `ipw_X_w_01'!=. & att22==1
	
	tempvar psi11_summand
	gen double `psi11_summand' = `ipw_X_w_11'*(std_cesd_age40 - `mu1_CML') ///
		+ `ipw_D1'*(`u_11_CL' - (`u_11_CL0'*`phat_L0_CD1' + `u_11_CL1'*`phat_L1_CD1')) ///
		+ `ipw_D1'*((`mu1_CML0'*`phat_L0_CD1' + `mu1_CML1'*`phat_L1_CD1') - `v_CD11') ///
		+ `v_CD11'
	
	tempvar psi00_summand
	gen double `psi00_summand' = `ipw_X_w_00'*(std_cesd_age40 - `mu0_CML') ///
		+ `ipw_D0'*(`u_00_CL' - (`u_00_CL0'*`phat_L0_CD0' + `u_00_CL1'*`phat_L1_CD0')) ///
		+ `ipw_D0'*((`mu0_CML0'*`phat_L0_CD0' + `mu0_CML1'*`phat_L1_CD0') - `v_CD00') ///
		+ `v_CD00'

	tempvar psi01_summand
	gen double `psi01_summand' = `ipw_X_w_01'*(std_cesd_age40 - `mu1_CML') ///
		+ `ipw_D1'*(`u_01_CL' - (`u_01_CL0'*`phat_L0_CD1' + `u_01_CL1'*`phat_L1_CD1')) ///
		+ `ipw_D0'*((`mu1_CML0'*`phat_L0_CD1' + `mu1_CML1'*`phat_L1_CD1') - `v_CD01') ///
		+ `v_CD01'
			
	reg `psi11_summand'
	return scalar psi11 = _b[_cons]
	
	reg `psi00_summand'
	return scalar psi00 = _b[_cons]

	reg `psi01_summand'
	return scalar psi01 = _b[_cons]
		
end mrvent

//parametric multiply robust estimates of interventional effects
qui bootstrap ///
	OE=(r(psi11)-r(psi00)) ///
	IDE=(r(psi01)-r(psi00)) ///
	IIE=(r(psi11)-r(psi01)), ///
	reps(2000) seed(60637): mrvent

estat bootstrap, p noheader

//define DML estimator for interventional effects

//note that stata does not currently support use of a superLearner
//we therefore implement the DML estimator using only random forests

ssc install rforest

capture program drop dmlvent
program define dmlvent, rclass

	tempvar dvar_orig lvar_orig
	qui gen `dvar_orig' = att22
	qui gen `lvar_orig' = ever_unemp_age3539
	
	tempvar kpart
	qui set seed 60637
	qui gen u = uniform()
	qui sort u
	qui gen `kpart' = ceil(_n/(_N/5))
	qui drop u

	local tvars	///
		phat_D0_C phat_D1_C phat_D_C ///
		phat_D0_CM phat_D1_CM phat_D_CM ///
		phat_L0_CD phat_L1_CD phat_L_CD ///
		phat_L0_CD1 phat_L1_CD1 phat_L_CD1 ///
		phat_L0_CD0 phat_L1_CD0 phat_L_CD0 ///
		phat_L0_CDM phat_L1_CDM phat_L_CDM ///
		phat_L0_CD1M phat_L1_CD1M phat_L_CD1M ///
		phat_L0_CD0M phat_L1_CD0M phat_L_CD0M ///
		mu_CDLM mu1_CML mu0_CML mu0_CML0 mu0_CML1 mu1_CML0 mu1_CML1 ///
		w_00 w_11 w_01 ///
		u_00_CL u_00_CL0 u_00_CL1 ///
		u_11_CL u_11_CL0 u_11_CL1 ///
		u_01_CL u_01_CL0 u_01_CL1 ///
		v_CD11 v_CD01 ///
		v_CD00
		
	tempvar `tvars'
	
	foreach v in `tvars' {
		qui gen double ``v'' = .
		}

	forval k=1/5 {

		qui rforest att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

		tempvar	xxphat_D0_C xxphat_D1_C xxphat_D_C
		qui predict `xxphat_D0_C' `xxphat_D1_C', pr
		qui gen `xxphat_D_C' = att22*`xxphat_D1_C' + (1-att22)*`xxphat_D0_C'
		qui replace `phat_D0_C' = `xxphat_D0_C' if `kpart'==`k'
		qui replace `phat_D1_C' = `xxphat_D1_C' if `kpart'==`k'
		qui replace `phat_D_C' = `xxphat_D_C' if `kpart'==`k'
		
		qui rforest att22 log_faminc_adj_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(class) iter(200) lsize(40) seed(60637)
				
		tempvar xxphat_D0_CM xxphat_D1_CM xxphat_D_CM
		qui predict `xxphat_D0_CM' `xxphat_D1_CM', pr
		qui gen `xxphat_D_CM' = att22*`xxphat_D1_CM' + (1-att22)*`xxphat_D0_CM'
		qui replace `phat_D0_CM' = `xxphat_D0_CM' if `kpart'==`k'
		qui replace `phat_D1_CM' = `xxphat_D1_CM' if `kpart'==`k'
		qui replace `phat_D_CM' = `xxphat_D_CM' if `kpart'==`k'

		qui rforest ever_unemp_age3539 att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)		

		tempvar xxphat_L0_CD xxphat_L1_CD xxphat_L_CD
		qui predict `xxphat_L0_CD' `xxphat_L1_CD', pr
		qui gen `xxphat_L_CD' = ever_unemp_age3539*`xxphat_L1_CD' + (1-ever_unemp_age3539)*`xxphat_L0_CD'
		qui replace `phat_L0_CD' = `xxphat_L0_CD' if `kpart'==`k'
		qui replace `phat_L1_CD' = `xxphat_L1_CD' if `kpart'==`k'
		qui replace `phat_L_CD' = `xxphat_L_CD' if `kpart'==`k'
	
		qui replace att22 = 1
		tempvar xxphat_L0_CD1 xxphat_L1_CD1 xxphat_L_CD1 
		qui predict `xxphat_L0_CD1' `xxphat_L1_CD1', pr
		qui gen `xxphat_L_CD1' = ever_unemp_age3539*`xxphat_L1_CD1' + (1-ever_unemp_age3539)*`xxphat_L0_CD1'
		qui replace `phat_L0_CD1' = `xxphat_L0_CD1' if `kpart'==`k'
		qui replace `phat_L1_CD1' = `xxphat_L1_CD1' if `kpart'==`k'
		qui replace `phat_L_CD1' = `xxphat_L_CD1' if `kpart'==`k'

		qui replace att22 = 0
		tempvar xxphat_L0_CD0 xxphat_L1_CD0 xxphat_L_CD0 
		qui predict `xxphat_L0_CD0' `xxphat_L1_CD0', pr
		qui gen `xxphat_L_CD0' = ever_unemp_age3539*`xxphat_L1_CD0' + (1-ever_unemp_age3539)*`xxphat_L0_CD0'
		qui replace `phat_L0_CD0' = `xxphat_L0_CD0' if `kpart'==`k'
		qui replace `phat_L1_CD0' = `xxphat_L1_CD0' if `kpart'==`k'
		qui replace `phat_L_CD0' = `xxphat_L_CD0' if `kpart'==`k'

		qui replace att22 = `dvar_orig'
		
		qui rforest ever_unemp_age3539 log_faminc_adj_age3539 att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(class) iter(200) lsize(40) seed(60637)		
	
		tempvar xxphat_L0_CDM xxphat_L1_CDM xxphat_L_CDM 
		qui predict `xxphat_L0_CDM' `xxphat_L1_CDM', pr
		qui gen `xxphat_L_CDM' = ever_unemp_age3539*`xxphat_L1_CDM' + (1-ever_unemp_age3539)*`xxphat_L0_CDM'
		qui replace `phat_L0_CDM' = `xxphat_L0_CDM' if `kpart'==`k'
		qui replace `phat_L1_CDM' = `xxphat_L1_CDM' if `kpart'==`k'
		qui replace `phat_L_CDM' = `xxphat_L_CDM' if `kpart'==`k'
			
		qui replace att22 = 1
		tempvar xxphat_L0_CD1M xxphat_L1_CD1M xxphat_L_CD1M
		qui predict `xxphat_L0_CD1M' `xxphat_L1_CD1M', pr
		qui gen `xxphat_L_CD1M' = ever_unemp_age3539*`xxphat_L1_CD1M' + (1-ever_unemp_age3539)*`xxphat_L0_CD1M'
		qui replace `phat_L0_CD1M' = `xxphat_L0_CD1M' if `kpart'==`k'
		qui replace `phat_L1_CD1M' = `xxphat_L1_CD1M' if `kpart'==`k'
		qui replace `phat_L_CD1M' = `xxphat_L_CD1M' if `kpart'==`k'

		qui replace att22 = 0
		tempvar xxphat_L0_CD0M xxphat_L1_CD0M xxphat_L_CD0M
		qui predict `xxphat_L0_CD0M' `xxphat_L1_CD0M', pr
		qui gen `xxphat_L_CD0M' = ever_unemp_age3539*`xxphat_L1_CD0M' + (1-ever_unemp_age3539)*`xxphat_L0_CD0M'
		qui replace `phat_L0_CD0M' = `xxphat_L0_CD0M' if `kpart'==`k'
		qui replace `phat_L1_CD0M' = `xxphat_L1_CD0M' if `kpart'==`k'
		qui replace `phat_L_CD0M' = `xxphat_L_CD0M' if `kpart'==`k'

		qui replace att22 = `dvar_orig'				
		
		qui rforest std_cesd_age40 ever_unemp_age3539 log_faminc_adj_age3539 att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(10) seed(60637)		
				
		tempvar xxmu_CDLM
		qui predict `xxmu_CDLM'
		qui replace `mu_CDLM' = `xxmu_CDLM' if `kpart'==`k'
			
		qui replace att22 = 1
		tempvar xxmu1_CML
		qui predict `xxmu1_CML'
		qui replace `mu1_CML' = `xxmu1_CML' if `kpart'==`k'
	
		qui replace att22 = 0
		tempvar xxmu0_CML
		qui predict `xxmu0_CML'
		qui replace `mu0_CML' = `xxmu0_CML' if `kpart'==`k'
	
		qui replace att22 = 0
		qui replace ever_unemp_age3539 = 0
		tempvar xxmu0_CML0
		qui predict `xxmu0_CML0'
		qui replace `mu0_CML0' = `xxmu0_CML0' if `kpart'==`k'
	
		qui replace att22 = 0
		qui replace ever_unemp_age3539 = 1
		tempvar xxmu0_CML1
		qui predict `xxmu0_CML1'
		qui replace `mu0_CML1' = `xxmu0_CML1' if `kpart'==`k'

		qui replace att22 = 1
		qui replace ever_unemp_age3539 = 0
		tempvar xxmu1_CML0
		qui predict `xxmu1_CML0'
		qui replace `mu1_CML0' = `xxmu1_CML0' if `kpart'==`k'
	
		qui replace att22 = 1
		qui replace ever_unemp_age3539 = 1
		tempvar xxmu1_CML1
		qui predict `xxmu1_CML1'
		qui replace `mu1_CML1' = `xxmu1_CML1' if `kpart'==`k'
	
		qui replace att22 = `dvar_orig'
		qui replace ever_unemp_age3539 = `lvar_orig'
			
		tempvar xxw_00 xxw_11 xxw_01 
		qui gen double `xxw_11' = (`xxphat_L_CD1') / (`xxphat_L_CD1M')
		qui gen double `xxw_00' = (`xxphat_L_CD0') / (`xxphat_L_CD0M')
		qui gen double `xxw_01' = (`xxphat_D1_C'*`xxphat_L_CD1'*`xxphat_D0_CM') / (`xxphat_D0_C'*`xxphat_L_CD1M'*`xxphat_D1_CM')
		
		qui replace `w_11' = `xxw_11' if `kpart'==`k'
		qui replace `w_00' = `xxw_00' if `kpart'==`k'
		qui replace `w_01' = `xxw_01' if `kpart'==`k'
		
		foreach w in `w_11' `w_00' `w_01' {
			qui centile `w' if `w'!=., c(1 99) 
			qui replace `w'=r(c_1) if `w'<r(c_1) & `w'!=.
			qui replace `w'=r(c_2) if `w'>r(c_2) & `w'!=.
		}
	
		tempvar u_00_regressand u_11_regressand u_01_regressand
		qui gen double `u_00_regressand' = `xxmu0_CML'*`xxw_00' if `kpart'!=`k'
		qui gen double `u_11_regressand' = `xxmu1_CML'*`xxw_11' if `kpart'!=`k'
		qui gen double `u_01_regressand' = `xxmu1_CML'*`xxw_01' if `kpart'!=`k'
		
		qui rforest `u_00_regressand' att22 ever_unemp_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(10) seed(60637)		
				
		qui replace att22 = 0
		tempvar xxu_00_CL
		qui predict `xxu_00_CL'
		qui replace `u_00_CL' = `xxu_00_CL' if `kpart'==`k'
	
		qui replace att22 = 0
		qui replace ever_unemp_age3539 = 0
		tempvar xxu_00_CL0
		qui predict `xxu_00_CL0'
		qui replace `u_00_CL0' = `xxu_00_CL0' if `kpart'==`k'
			
		qui replace att22 = 0
		qui replace ever_unemp_age3539 = 1
		tempvar xxu_00_CL1
		qui predict `xxu_00_CL1'
		qui replace `u_00_CL1' = `xxu_00_CL1' if `kpart'==`k'
	
		qui replace att22 = `dvar_orig'
		qui replace ever_unemp_age3539 = `lvar_orig'

		qui rforest `u_11_regressand' att22 ever_unemp_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(10) seed(60637)		
				
		qui replace att22 = 1
		tempvar xxu_11_CL
		qui predict `xxu_11_CL'
		qui replace `u_11_CL' = `xxu_11_CL' if `kpart'==`k'
	
		qui replace att22 = 1
		qui replace ever_unemp_age3539 = 0
		tempvar xxu_11_CL0
		qui predict `xxu_11_CL0'
		qui replace `u_11_CL0' = `xxu_11_CL0' if `kpart'==`k'
	
		qui replace att22 = 1
		qui replace ever_unemp_age3539 = 1
		tempvar xxu_11_CL1
		qui predict `xxu_11_CL1'
		qui replace `u_11_CL1' = `xxu_11_CL1' if `kpart'==`k'
		
		qui replace att22 = `dvar_orig'
		qui replace ever_unemp_age3539 = `lvar_orig'

		qui rforest `u_01_regressand' att22 ever_unemp_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(10) seed(60637)		
				
		qui replace att22 = 1
		tempvar xxu_01_CL
		qui predict `xxu_01_CL'
		qui replace `u_01_CL' = `xxu_01_CL' if `kpart'==`k'
	
		qui replace att22 = 1
		qui replace ever_unemp_age3539 = 0
		tempvar xxu_01_CL0
		qui predict `xxu_01_CL0'
		qui replace `u_01_CL0' = `xxu_01_CL0' if `kpart'==`k'
	
		qui replace att22 = 1
		qui replace ever_unemp_age3539 = 1
		tempvar xxu_01_CL1
		qui predict `xxu_01_CL1'
		qui replace `u_01_CL1' = `xxu_01_CL1' if `kpart'==`k'
		
		qui replace att22 = `dvar_orig'
		qui replace ever_unemp_age3539 = `lvar_orig'
	
		tempvar v_CD1_regressand v_CD0_regressand
		qui gen double `v_CD1_regressand' = `xxmu1_CML0'*`xxphat_L0_CD1' + `xxmu1_CML1'*`xxphat_L1_CD1' if `kpart'!=`k'
		qui gen double `v_CD0_regressand' = `xxmu0_CML0'*`xxphat_L0_CD0' + `xxmu0_CML1'*`xxphat_L1_CD0' if `kpart'!=`k'

		qui rforest `v_CD1_regressand' att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(10) seed(60637)		
				
		tempvar xxv_CD11 xxv_CD01

		qui replace att22 = 1
		qui predict `xxv_CD11'
		qui replace `v_CD11' = `xxv_CD11' if `kpart'==`k'
	
		qui replace att22 = 0
		qui predict `xxv_CD01'
		qui replace `v_CD01' = `xxv_CD01' if `kpart'==`k'
	
		qui replace att22 = `dvar_orig'
	
		qui rforest `v_CD0_regressand' att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(10) seed(60637)		

		tempvar xxv_CD00
	
		qui replace att22 = 0
		qui predict `xxv_CD00'
		qui replace `v_CD00' = `xxv_CD00' if `kpart'==`k'
	
		qui replace att22 = `dvar_orig'
		
		capture drop ///
			`xxphat_D0_C' `xxphat_D1_C' `xxphat_D_C' ///
			`xxphat_D0_CM' `xxphat_D1_CM' `xxphat_D_CM' ///
			`xxphat_L0_CD' `xxphat_L1_CD' `xxphat_L_CD' ///
			`xxphat_L0_CD1' `xxphat_L1_CD1' `xxphat_L_CD1' ///
			`xxphat_L0_CD0' `xxphat_L1_CD0' `xxphat_L_CD0' ///
			`xxphat_L0_CDM' `xxphat_L1_CDM' `xxphat_L_CDM' ///
			`xxphat_L0_CD1M' `xxphat_L1_CD1M' `xxphat_L_CD1M' ///
			`xxphat_L0_CD0M' `xxphat_L1_CD0M' `xxphat_L_CD0M' ///
			`xxmu_CDLM' ///
			`xxmu1_CML' `xxmu0_CML' ///
			`xxmu0_CML0' `xxmu0_CML1' ///
			`xxmu1_CML0' `xxmu1_CML1' ///
			`xxw_00' `xxw_11' `xxw_01' ///
			`u_00_regressand' `u_11_regressand' `u_01_regressand' ///
			`xxu_00_CL' `xxu_00_CL0' `xxu_00_CL1' ///
			`xxu_11_CL' `xxu_11_CL0' `xxu_11_CL1' ///
			`xxu_01_CL' `xxu_01_CL0' `xxu_01_CL1' ///
			`xxv_CD11' `xxv_CD01' `xxv_CD00'
			
		}
	
	tempvar ipw_D1 ipw_D0 
	qui gen double `ipw_D1' = att22/`phat_D1_C'
	qui gen double `ipw_D0' = (1-att22)/`phat_D0_C'
	
	tempvar ipw_X_w_11 ipw_X_w_00 ipw_X_w_01
	qui gen double `ipw_X_w_11' = `ipw_D1'*`w_11'
	qui gen double `ipw_X_w_00' = `ipw_D0'*`w_00'
	qui gen double `ipw_X_w_01' = `ipw_D1'*`w_01'
	
	qui centile `ipw_D0' if `ipw_D0'!=. & att22==0, c(1 99) 
	qui replace `ipw_D0'=r(c_1) if `ipw_D0'<r(c_1) & `ipw_D0'!=. & att22==0
	qui replace `ipw_D0'=r(c_2) if `ipw_D0'>r(c_2) & `ipw_D0'!=. & att22==0
	
	qui centile `ipw_D1' if `ipw_D1'!=. & att22==1, c(1 99) 
	qui replace `ipw_D1'=r(c_1) if `ipw_D1'<r(c_1) & `ipw_D1'!=. & att22==1
	qui replace `ipw_D1'=r(c_2) if `ipw_D1'>r(c_2) & `ipw_D1'!=. & att22==1
	
	qui centile `ipw_X_w_00' if `ipw_X_w_00'!=. & att22==0, c(1 99) 
	qui replace `ipw_X_w_00'=r(c_1) if `ipw_X_w_00'<r(c_1) & `ipw_X_w_00'!=. & att22==0
	qui replace `ipw_X_w_00'=r(c_2) if `ipw_X_w_00'>r(c_2) & `ipw_X_w_00'!=. & att22==0
		
	qui centile `ipw_X_w_11' if `ipw_X_w_11'!=. & att22==1, c(1 99) 
	qui replace `ipw_X_w_11'=r(c_1) if `ipw_X_w_11'<r(c_1) & `ipw_X_w_11'!=. & att22==1
	qui replace `ipw_X_w_11'=r(c_2) if `ipw_X_w_11'>r(c_2) & `ipw_X_w_11'!=. & att22==1
		
	qui centile `ipw_X_w_01' if `ipw_X_w_01'!=. & att22==1, c(1 99) 
	qui replace `ipw_X_w_01'=r(c_1) if `ipw_X_w_01'<r(c_1) & `ipw_X_w_01'!=. & att22==1
	qui replace `ipw_X_w_01'=r(c_2) if `ipw_X_w_01'>r(c_2) & `ipw_X_w_01'!=. & att22==1
		
	tempvar psi11_summand
	qui gen double `psi11_summand' = `ipw_X_w_11'*(std_cesd_age40 - `mu1_CML') ///
		+ `ipw_D1'*(`u_11_CL' - (`u_11_CL0'*`phat_L0_CD1' + `u_11_CL1'*`phat_L1_CD1')) ///
		+ `ipw_D1'*((`mu1_CML0'*`phat_L0_CD1' + `mu1_CML1'*`phat_L1_CD1') - `v_CD11') ///
		+ `v_CD11'
	
	tempvar psi00_summand
	qui gen double `psi00_summand' = `ipw_X_w_00'*(std_cesd_age40 - `mu0_CML') ///
		+ `ipw_D0'*(`u_00_CL' - (`u_00_CL0'*`phat_L0_CD0' + `u_00_CL1'*`phat_L1_CD0')) ///
		+ `ipw_D0'*((`mu0_CML0'*`phat_L0_CD0' + `mu0_CML1'*`phat_L1_CD0') - `v_CD00') ///
		+ `v_CD00'

	tempvar psi01_summand
	qui gen double `psi01_summand' = `ipw_X_w_01'*(std_cesd_age40 - `mu1_CML') ///
		+ `ipw_D1'*(`u_01_CL' - (`u_01_CL0'*`phat_L0_CD1' + `u_01_CL1'*`phat_L1_CD1')) ///
		+ `ipw_D0'*((`mu1_CML0'*`phat_L0_CD1' + `mu1_CML1'*`phat_L1_CD1') - `v_CD01') ///
		+ `v_CD01'
	
	qui gen OE = `psi11_summand' - `psi00_summand'
	qui gen IDE = `psi01_summand' - `psi00_summand'
	qui gen IIE = `psi11_summand' - `psi01_summand'
		
	mean OE IDE IIE, noheader
		
	qui drop OE IDE IIE
		
end dmlvent

//compute DML estimates of interventional effects 
dmlvent

log close

//note that the estimates differ from those reported in the text, which are 
//based on the R implementation. This is due to differences in how the 
//weights are censored and to the use of random forests exclusively rather 
//than a super learner in the DML estimators.
