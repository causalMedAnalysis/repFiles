/*Table 6.5*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch6\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

log using "${logdir}table_6-5.log", replace 

//input data
use "${datadir}Tatar\tatar.dta", clear

/****************************************************************
define parametric MR estimator for path-specific effects (type 2)
*****************************************************************/
capture program drop mrpath
program define mrpath, rclass

	local cvars ///
		kulak prosoviet_pre religiosity_pre land_pre orchard_pre ///
		animals_pre carriage_pre otherprop_pre
	local dvar violence
	local yvar annex
	local mvar1 trust_g1 victim_g1 fear_g1
	local mvar2 trust_g2 victim_g2 fear_g2
	local mvar3 trust_g3 victim_g3 fear_g3

	tempvar dvar_orig
	qui gen `dvar_orig' = `dvar'
	
	logit `dvar' `cvars'
	
	tempvar pi1_C
	predict `pi1_C', pr
		
	logit `dvar' `mvar1' `cvars' 
	
	tempvar pi1_CM1
	predict `pi1_CM1', pr

	logit `dvar' `mvar1' `mvar2' `cvars'
	
	tempvar pi1_CM1M2
	predict `pi1_CM1M2', pr

	logit `dvar' `mvar1' `mvar2' `mvar3' `cvars'
	
	tempvar pi1_CM1M2M3
	predict `pi1_CM1M2M3', pr

	logit `yvar' `dvar' `mvar1' `mvar2' `mvar3' `cvars'

	replace `dvar' = 1
	tempvar mu1_CM1M2M3
	predict `mu1_CM1M2M3', pr
	
	replace `dvar' = 0
	tempvar mu0_CM1M2M3
	predict `mu0_CM1M2M3', pr
	
	replace `dvar' = `dvar_orig'
	
	glm `mu1_CM1M2M3' `dvar' `mvar1' `mvar2' `cvars', family(binomial) link(logit)
	
	replace `dvar' = 1
	tempvar nu11_CM1M2
	predict `nu11_CM1M2'
	
	replace `dvar' = 0
	tempvar nu01_CM1M2
	predict `nu01_CM1M2'

	replace `dvar' = `dvar_orig'

	glm `mu0_CM1M2M3' `dvar' `mvar1' `mvar2' `cvars', family(binomial) link(logit)
		
	replace `dvar' = 0
	tempvar nu00_CM1M2
	predict `nu00_CM1M2'

	replace `dvar' = `dvar_orig'
	
	glm `nu00_CM1M2' `dvar' `mvar1' `cvars', family(binomial) link(logit)
		
	replace `dvar' = 0
	tempvar xi000_CM1
	predict `xi000_CM1'
	
	replace `dvar' = `dvar_orig'
	
	glm `nu01_CM1M2' `dvar' `mvar1' `cvars', family(binomial) link(logit)
		
	replace `dvar' = 0
	tempvar xi001_CM1
	predict `xi001_CM1'
	
	replace `dvar' = `dvar_orig'
	
	glm `nu11_CM1M2' `dvar' `mvar1' `cvars', family(binomial) link(logit)
		
	replace `dvar' = 0
	tempvar xi011_CM1
	predict `xi011_CM1'

	replace `dvar' = 1
	tempvar xi111_CM1
	predict `xi111_CM1'
	
	replace `dvar' = `dvar_orig'
	
	glm `xi000_CM1' `dvar' `cvars', family(binomial) link(logit)
	
	replace `dvar' = 0
	tempvar zi0000_C
	predict `zi0000_C'
	
	replace `dvar' = `dvar_orig'
	
	glm `xi001_CM1' `dvar' `cvars', family(binomial) link(logit)
	
	replace `dvar' = 0
	tempvar zi0001_C
	predict `zi0001_C'
	
	replace `dvar' = `dvar_orig'
	
	glm `xi011_CM1' `dvar' `cvars', family(binomial) link(logit)
	
	replace `dvar' = 0
	tempvar zi0011_C
	predict `zi0011_C'
	
	replace `dvar' = `dvar_orig'
	
	glm `xi111_CM1' `dvar' `cvars', family(binomial) link(logit)
	
	replace `dvar' = 0
	tempvar zi0111_C
	predict `zi0111_C'

	replace `dvar' = 1
	tempvar zi1111_C
	predict `zi1111_C'
	
	replace `dvar' = `dvar_orig'
	
	tempvar ipw_term1_0000 ipw_term1_0001 ipw_term1_0011 ipw_term1_0111 ipw_term1_1111
	gen `ipw_term1_0000' = (1-`dvar')/(1 - `pi1_C')
	gen `ipw_term1_0001' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1M2M3')/`pi1_CM1M2M3')
	gen `ipw_term1_0011' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1M2')/`pi1_CM1M2')
	gen `ipw_term1_0111' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1')/`pi1_CM1')
	gen `ipw_term1_1111' = (`dvar'/`pi1_C')
	
	tempvar ipw_term2_000 ipw_term2_001 ipw_term2_011 ipw_term2_111
	gen `ipw_term2_000' = ((1-`dvar')/(1 - `pi1_C'))
	gen `ipw_term2_001' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1M2')/`pi1_CM1M2')
	gen `ipw_term2_011' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1')/`pi1_CM1')
	gen `ipw_term2_111' = (`dvar'/`pi1_C')

	tempvar ipw_term3_00 ipw_term3_01 ipw_term3_11
	gen `ipw_term3_00' = (1-`dvar')/(1 - `pi1_C')
	gen `ipw_term3_01' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1') / `pi1_CM1')
	gen `ipw_term3_11' = (`dvar'/`pi1_C')

	tempvar ipw_term4_0 ipw_term4_1
	gen `ipw_term4_0' = (1-`dvar')/(1 - `pi1_C')
	gen `ipw_term4_1' = (`dvar'/`pi1_C')
  
	foreach v in `ipw_term1_0000' `ipw_term2_000' `ipw_term3_00' `ipw_term4_0' {
		centile `v' if `v'!=. & `dvar'==0, c(1 99) 
		replace `v'=r(c_1) if `v'<r(c_1) & `v'!=. & `dvar'==0
		replace `v'=r(c_2) if `v'>r(c_2) & `v'!=. & `dvar'==0
		sum `v' if `dvar'==0
		}

	foreach v in ///
		`ipw_term1_0001' `ipw_term1_0011' `ipw_term1_0111' `ipw_term1_1111' ///
		`ipw_term2_001' `ipw_term2_011' `ipw_term2_111' ///
		`ipw_term3_01' `ipw_term3_11' ///
		`ipw_term4_1' {
			centile `v' if `v'!=. & `dvar'==1, c(1 99) 
			replace `v'=r(c_1) if `v'<r(c_1) & `v'!=. & `dvar'==1
			replace `v'=r(c_2) if `v'>r(c_2) & `v'!=. & `dvar'==1
			sum `v' if `dvar'==1
			}
		
	tempvar psi0000_summand
	gen `psi0000_summand' = `ipw_term1_0000'*(`yvar' - `mu0_CM1M2M3') ///
		+ `ipw_term2_000'*(`mu0_CM1M2M3' - `nu00_CM1M2') ///
		+ `ipw_term3_00'*(`nu00_CM1M2' - `xi000_CM1') ///
		+ `ipw_term4_0'*(`xi000_CM1' - `zi0000_C') ///
		+ `zi0000_C'

	tempvar psi0001_summand
	gen `psi0001_summand' = `ipw_term1_0001'*(`yvar' - `mu1_CM1M2M3') ///
		+ `ipw_term2_000'*(`mu1_CM1M2M3' - `nu01_CM1M2') ///
		+ `ipw_term3_00'*(`nu01_CM1M2' - `xi001_CM1') ///
		+ `ipw_term4_0'*(`xi001_CM1' - `zi0001_C') ///
		+ `zi0001_C'
		
	tempvar psi0011_summand
	gen `psi0011_summand' = `ipw_term1_0011'*(`yvar' - `mu1_CM1M2M3') ///
		+ `ipw_term2_001'*(`mu1_CM1M2M3' - `nu11_CM1M2') ///
		+ `ipw_term3_00'*(`nu11_CM1M2' - `xi011_CM1') ///
		+ `ipw_term4_0'*(`xi011_CM1' - `zi0011_C') ///
		+ `zi0011_C'

	tempvar psi0111_summand
	gen `psi0111_summand' = `ipw_term1_0111'*(`yvar' - `mu1_CM1M2M3') ///
		+ `ipw_term2_011'*(`mu1_CM1M2M3' - `nu11_CM1M2') ///
		+ `ipw_term3_01'*(`nu11_CM1M2' - `xi111_CM1') ///
		+ `ipw_term4_0'*(`xi111_CM1' - `zi0111_C') ///
		+ `zi0111_C'

	tempvar psi1111_summand
	gen `psi1111_summand' = `ipw_term1_1111'*(`yvar' - `mu1_CM1M2M3') ///
		+ `ipw_term2_111'*(`mu1_CM1M2M3' - `nu11_CM1M2') ///
		+ `ipw_term3_11'*(`nu11_CM1M2' - `xi111_CM1') ///
		+ `ipw_term4_1'*(`xi111_CM1' - `zi1111_C') ///
		+ `zi1111_C'
		
	reg `psi0000_summand'
	return scalar psi0000 = _b[_cons]

	reg `psi0001_summand'
	return scalar psi0001 = _b[_cons]
	
	reg `psi0011_summand'
	return scalar psi0011 = _b[_cons]
	
	reg `psi0111_summand'
	return scalar psi0111 = _b[_cons]

	reg `psi1111_summand'
	return scalar psi1111 = _b[_cons]
	
end mrpath

//parametric MR estimates of path-specific effects (type 2)
qui bootstrap ///
	ATE=(r(psi1111)-r(psi0000)) ///
	PSE_DY=(r(psi0001)-r(psi0000)) ///
	PSE_DM3Y=(r(psi0011)-r(psi0001)) ///
	PSE_DM2Y=(r(psi0111)-r(psi0011)) ///
	PSE_DM1Y=(r(psi1111)-r(psi0111)), ///
	reps(2000) seed(60637): mrpath

estat bootstrap, p noheader

/******************************************************
define DML estimator for path-specific effects (type 2)
*******************************************************/
//note that stata does not currently support use of a superLearner
//we therefore implement the DML estimator using only random forests

*ssc install rforest //install rforest module if not already installed

capture program drop dmlpath
program define dmlpath, rclass

	local cvars ///
		kulak prosoviet_pre religiosity_pre land_pre orchard_pre ///
		animals_pre carriage_pre otherprop_pre
	local dvar violence
	local yvar annex
	local mvar1 trust_g1 victim_g1 fear_g1
	local mvar2 trust_g2 victim_g2 fear_g2
	local mvar3 trust_g3 victim_g3 fear_g3

	tempvar dvar_orig
	qui gen `dvar_orig' = `dvar'
	
	tempvar u kpart
	set seed 60637
	qui gen `u' = uniform()
	qui sort `u'
	qui gen `kpart' = ceil(_n/(_N/5))

	local tvars	///
		pi1_C pi1_CM1 ///
		pi1_CM1M2 pi1_CM1M2M3 ///
		mu1_CM1M2M3 mu0_CM1M2M3 ///
		nu11_CM1M2 nu01_CM1M2 nu00_CM1M2 ///
		xi000_CM1 xi001_CM1 xi011_CM1 xi111_CM1 ///
		zi0000_C zi0001_C zi0011_C zi0111_C zi1111_C
		
	tempvar `tvars'
	
	foreach v in `tvars' {
		qui gen ``v'' = .
		}

	forval k=1/5 {
	
		qui rforest `dvar' `cvars' ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

		tempvar	xxpi0_C xxpi1_C
		qui predict `xxpi0_C' `xxpi1_C' if `kpart'==`k', pr
		qui replace `pi1_C' = `xxpi1_C' if `kpart'==`k'
		
		qui rforest `dvar' `mvar1' `cvars' ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

		tempvar	xxpi0_CM1 xxpi1_CM1
		qui predict `xxpi0_CM1' `xxpi1_CM1' if `kpart'==`k', pr
		qui replace `pi1_CM1' = `xxpi1_CM1' if `kpart'==`k'
		
		qui rforest `dvar' `mvar1' `mvar2' `cvars' ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

		tempvar	xxpi0_CM1M2 xxpi1_CM1M2
		qui predict `xxpi0_CM1M2' `xxpi1_CM1M2' if `kpart'==`k', pr
		qui replace `pi1_CM1M2' = `xxpi1_CM1M2' if `kpart'==`k'
		
		qui rforest `dvar' `mvar1' `mvar2' `mvar3' `cvars' ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

		tempvar	xxpi0_CM1M2M3 xxpi1_CM1M2M3
		qui predict `xxpi0_CM1M2M3' `xxpi1_CM1M2M3' if `kpart'==`k', pr
		qui replace `pi1_CM1M2M3' = `xxpi1_CM1M2M3' if `kpart'==`k'
		
		qui rforest `yvar' `dvar' `mvar1' `mvar2' `mvar3' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)

		qui replace `dvar' = 1
		tempvar xxmu1_CM1M2M3
		qui predict `xxmu1_CM1M2M3'
		qui replace `mu1_CM1M2M3' = `xxmu1_CM1M2M3' if `kpart'==`k'
	
		qui replace `dvar' = 0
		tempvar xxmu0_CM1M2M3_ xxmu0_CM1M2M3
		qui predict `xxmu0_CM1M2M3'
		qui replace `mu0_CM1M2M3' = `xxmu0_CM1M2M3' if `kpart'==`k'
			
		qui replace `dvar' = `dvar_orig'
	
		qui rforest `xxmu1_CM1M2M3' `dvar' `mvar1' `mvar2' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
		qui replace `dvar' = 1
		tempvar xxnu11_CM1M2
		qui predict `xxnu11_CM1M2'
		qui replace `nu11_CM1M2' = `xxnu11_CM1M2' if `kpart'==`k'
	
		qui replace `dvar' = 0
		tempvar xxnu01_CM1M2
		qui predict `xxnu01_CM1M2'
		qui replace `nu01_CM1M2' = `xxnu01_CM1M2' if `kpart'==`k'

		qui replace `dvar' = `dvar_orig'

		qui rforest `xxmu0_CM1M2M3' `dvar' `mvar1' `mvar2' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace `dvar' = 0
		tempvar xxnu00_CM1M2
		qui predict `xxnu00_CM1M2'
		qui replace `nu00_CM1M2' = `xxnu00_CM1M2' if `kpart'==`k'

		qui replace `dvar' = `dvar_orig'
	
		qui rforest `xxnu00_CM1M2' `dvar' `mvar1' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace `dvar' = 0
		tempvar xxxi000_CM1
		qui predict `xxxi000_CM1'
		qui replace `xi000_CM1' = `xxxi000_CM1' if `kpart'==`k'
	
		qui replace `dvar' = `dvar_orig'
	
		qui rforest `xxnu01_CM1M2' `dvar' `mvar1' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace `dvar' = 0
		tempvar xxxi001_CM1
		qui predict `xxxi001_CM1'
		qui replace `xi001_CM1' = `xxxi001_CM1' if `kpart'==`k'
	
		qui replace `dvar' = `dvar_orig'
	
		qui rforest `xxnu11_CM1M2' `dvar' `mvar1' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace `dvar' = 0
		tempvar xxxi011_CM1
		qui predict `xxxi011_CM1'
		qui replace `xi011_CM1' = `xxxi011_CM1' if `kpart'==`k'

		qui replace `dvar' = 1
		tempvar xxxi111_CM1
		qui predict `xxxi111_CM1'
		qui replace `xi111_CM1' = `xxxi111_CM1' if `kpart'==`k'
	
		qui replace `dvar' = `dvar_orig'
	
		qui rforest `xxxi000_CM1' `dvar' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
		qui replace `dvar' = 0
		tempvar xxzi0000_C
		qui predict `xxzi0000_C'
		qui replace `zi0000_C' = `xxzi0000_C' if `kpart'==`k'
	
		qui replace `dvar' = `dvar_orig'
	
		qui rforest `xxxi001_CM1' `dvar' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
		qui replace `dvar' = 0
		tempvar xxzi0001_C
		qui predict `xxzi0001_C'
		qui replace `zi0001_C' = `xxzi0001_C' if `kpart'==`k'
	
		qui replace `dvar' = `dvar_orig'
	
		qui rforest `xxxi011_CM1' `dvar' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
		qui replace `dvar' = 0
		tempvar xxzi0011_C
		qui predict `xxzi0011_C'
		qui replace `zi0011_C' = `xxzi0011_C' if `kpart'==`k'
	
		qui replace `dvar' = `dvar_orig'
	
		qui rforest `xxxi111_CM1' `dvar' `cvars' ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
		qui replace `dvar' = 0
		tempvar xxzi0111_C
		qui predict `xxzi0111_C'
		qui replace `zi0111_C' = `xxzi0111_C' if `kpart'==`k'

		qui replace `dvar' = 1
		tempvar xxzi1111_C
		qui predict `xxzi1111_C'
		qui replace `zi1111_C' = `xxzi1111_C' if `kpart'==`k'
	
		qui replace `dvar' = `dvar_orig'
		
		capture drop ///
			`xxpi0_C' `xxpi1_C' ///
			`xxpi0_CM1' `xxpi1_CM1' ///
			`xxpi0_CM1M2' `xxpi1_CM1M2' ///
			`xxpi0_CM1M2M3' `xxpi1_CM1M2M3' ///
			`xxmu1_CM1M2M3_' `xxmu1_CM1M2M3' `xxmu0_CM1M2M3_' `xxmu0_CM1M2M3' ///
			`xxnu11_CM1M2' `xxnu01_CM1M2' `xxnu00_CM1M2' ///
			`xxxi000_CM1' `xxxi001_CM1' `xxxi011_CM1' `xxxi111_CM1' ///
			`xxzi0000_C' `xxzi0001_C' `xxzi0011_C' `xxzi0111_C' `xxzi1111_C'
			
		}
	
	tempvar ipw_term1_0000 ipw_term1_0001 ipw_term1_0011 ipw_term1_0111 ipw_term1_1111
	qui gen `ipw_term1_0000' = (1-`dvar')/(1 - `pi1_C')
	qui gen `ipw_term1_0001' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1M2M3')/`pi1_CM1M2M3')
	qui gen `ipw_term1_0011' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1M2')/`pi1_CM1M2')
	qui gen `ipw_term1_0111' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1')/`pi1_CM1')
	qui gen `ipw_term1_1111' = (`dvar'/`pi1_C')
	
	tempvar ipw_term2_000 ipw_term2_001 ipw_term2_011 ipw_term2_111
	qui gen `ipw_term2_000' = ((1-`dvar')/(1 - `pi1_C'))
	qui gen `ipw_term2_001' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1M2')/`pi1_CM1M2')
	qui gen `ipw_term2_011' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1')/`pi1_CM1')
	qui gen `ipw_term2_111' = (`dvar'/`pi1_C')

	tempvar ipw_term3_00 ipw_term3_01 ipw_term3_11
	qui gen `ipw_term3_00' = (1-`dvar')/(1 - `pi1_C')
	qui gen `ipw_term3_01' = (`dvar'/(1 - `pi1_C')) * ((1 - `pi1_CM1') / `pi1_CM1')
	qui gen `ipw_term3_11' = (`dvar'/`pi1_C')

	tempvar ipw_term4_0 ipw_term4_1
	qui gen `ipw_term4_0' = (1-`dvar')/(1 - `pi1_C')
	qui gen `ipw_term4_1' = (`dvar'/`pi1_C')
  
	foreach v in `ipw_term1_0000' `ipw_term2_000' `ipw_term3_00' `ipw_term4_0' {
		qui centile `v' if `v'!=. & `dvar'==0, c(1 99) 
		qui replace `v'=r(c_1) if `v'<r(c_1) & `v'!=. & `dvar'==0
		qui replace `v'=r(c_2) if `v'>r(c_2) & `v'!=. & `dvar'==0
		}

	foreach v in ///
		`ipw_term1_0001' `ipw_term1_0011' `ipw_term1_0111' `ipw_term1_1111' ///
		`ipw_term2_001' `ipw_term2_011' `ipw_term2_111' ///
		`ipw_term3_01' `ipw_term3_11' ///
		`ipw_term4_1' {
			qui centile `v' if `v'!=. & `dvar'==1, c(1 99) 
			qui replace `v'=r(c_1) if `v'<r(c_1) & `v'!=. & `dvar'==1
			qui replace `v'=r(c_2) if `v'>r(c_2) & `v'!=. & `dvar'==1
			}
		
	tempvar psi0000_summand
	qui gen `psi0000_summand' = `ipw_term1_0000'*(`yvar' - `mu0_CM1M2M3') ///
		+ `ipw_term2_000'*(`mu0_CM1M2M3' - `nu00_CM1M2') ///
		+ `ipw_term3_00'*(`nu00_CM1M2' - `xi000_CM1') ///
		+ `ipw_term4_0'*(`xi000_CM1' - `zi0000_C') ///
		+ `zi0000_C'

	tempvar psi0001_summand
	qui gen `psi0001_summand' = `ipw_term1_0001'*(`yvar' - `mu1_CM1M2M3') ///
		+ `ipw_term2_000'*(`mu1_CM1M2M3' - `nu01_CM1M2') ///
		+ `ipw_term3_00'*(`nu01_CM1M2' - `xi001_CM1') ///
		+ `ipw_term4_0'*(`xi001_CM1' - `zi0001_C') ///
		+ `zi0001_C'
		
	tempvar psi0011_summand
	qui gen `psi0011_summand' = `ipw_term1_0011'*(`yvar' - `mu1_CM1M2M3') ///
		+ `ipw_term2_001'*(`mu1_CM1M2M3' - `nu11_CM1M2') ///
		+ `ipw_term3_00'*(`nu11_CM1M2' - `xi011_CM1') ///
		+ `ipw_term4_0'*(`xi011_CM1' - `zi0011_C') ///
		+ `zi0011_C'

	tempvar psi0111_summand
	qui gen `psi0111_summand' = `ipw_term1_0111'*(`yvar' - `mu1_CM1M2M3') ///
		+ `ipw_term2_011'*(`mu1_CM1M2M3' - `nu11_CM1M2') ///
		+ `ipw_term3_01'*(`nu11_CM1M2' - `xi111_CM1') ///
		+ `ipw_term4_0'*(`xi111_CM1' - `zi0111_C') ///
		+ `zi0111_C'

	tempvar psi1111_summand
	qui gen `psi1111_summand' = `ipw_term1_1111'*(`yvar' - `mu1_CM1M2M3') ///
		+ `ipw_term2_111'*(`mu1_CM1M2M3' - `nu11_CM1M2') ///
		+ `ipw_term3_11'*(`nu11_CM1M2' - `xi111_CM1') ///
		+ `ipw_term4_1'*(`xi111_CM1' - `zi1111_C') ///
		+ `zi1111_C'

	gen ATE = `psi1111_summand' - `psi0000_summand'
	gen PSE_DY = `psi0001_summand' - `psi0000_summand'
	gen PSE_DM3Y = `psi0011_summand' - `psi0001_summand'
	gen PSE_DM2Y = `psi0111_summand' - `psi0011_summand'
	gen PSE_DM1Y = `psi1111_summand' - `psi0111_summand'
		
	mean ATE PSE_DY PSE_DM3Y PSE_DM2Y PSE_DM1Y, noheader
		
	qui drop ATE PSE_DY PSE_DM3Y PSE_DM2Y PSE_DM1Y
		
end dmlpath

//DML estimates of path-specific effects (type 2)
dmlpath

log close
