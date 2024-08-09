/*Table 6.4*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch6\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

log using "${logdir}table_6-4.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if missing(cesd_age40, att22, ever_unemp_age3539, female, black, ///
	hispan, paredu, parprof, parinc_prank, famsize, afqt3, log_faminc_adj_age3539)

egen std_cesd_age40=std(cesd_age40)

/****************************************************************
define parametric MR estimator for path-specific effects (type 2)
*****************************************************************/
capture program drop mrpath
program define mrpath, rclass

	tempvar dvar_orig
	qui gen `dvar_orig' = att22
	
	logit att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar pi0_C pi1_C
	predict `pi1_C', pr
	qui gen `pi0_C' = 1 - `pi1_C'
		
	logit att22 ever_unemp_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar pi1_CM1 pi0_CM1
	predict `pi1_CM1', pr
	qui gen `pi0_CM1' = 1 - `pi1_CM1'

	logit att22 ever_unemp_age3539 log_faminc_adj_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	
	tempvar pi1_CM1M2 pi0_CM1M2
	predict `pi1_CM1M2', pr
	qui gen `pi0_CM1M2' = 1 - `pi1_CM1M2'
		
	reg std_cesd_age40 att22 ever_unemp_age3539 log_faminc_adj_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3

	replace att22 = 1
	tempvar mu1_CM1M2
	predict `mu1_CM1M2'
	
	replace att22 = 0
	tempvar mu0_CM1M2
	predict `mu0_CM1M2'
	
	replace att22 = `dvar_orig'
	
	reg `mu1_CM1M2' att22 ever_unemp_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
	
	replace att22 = 1
	tempvar nu11_CM1
	predict `nu11_CM1'
	
	replace att22 = 0
	tempvar nu01_CM1
	predict `nu01_CM1'

	replace att22 = `dvar_orig'

	reg `mu0_CM1M2' att22 ever_unemp_age3539 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
		
	replace att22 = 0
	tempvar nu00_CM1
	predict `nu00_CM1'

	replace att22 = `dvar_orig'
	
	reg `nu00_CM1' att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
		
	replace att22 = 0
	tempvar xi000_C
	predict `xi000_C'
	
	replace att22 = `dvar_orig'
	
	reg `nu01_CM1' att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
		
	replace att22 = 0
	tempvar xi001_C
	predict `xi001_C'
	
	replace att22 = `dvar_orig'
	
	reg `nu11_CM1' att22 ///
		female black hispan paredu parprof parinc_prank famsize afqt3
		
	replace att22 = 0
	tempvar xi011_C
	predict `xi011_C'

	replace att22 = 1
	tempvar xi111_C
	predict `xi111_C'
	
	replace att22 = `dvar_orig'
		
	tempvar ipw_term1_000 ipw_term1_001 ipw_term1_011 ipw_term1_111
	gen `ipw_term1_000' = (1-att22) / `pi0_C'
	gen `ipw_term1_001' = (att22 / `pi0_C') * ((`pi0_CM1M2'*`pi0_CM1') / (`pi1_CM1M2'*`pi0_CM1'))
	gen `ipw_term1_011' = (att22 / `pi0_C') * ((`pi1_CM1M2'*`pi0_CM1') / (`pi1_CM1M2'*`pi1_CM1'))
	gen `ipw_term1_111' = (att22 / `pi1_C')
	
	tempvar ipw_term2_00 ipw_term2_01 ipw_term2_11
	gen `ipw_term2_00' = (1-att22) / `pi0_C'
	gen `ipw_term2_01' = (att22 / `pi0_C') * (`pi0_CM1' / `pi1_CM1')
	gen `ipw_term2_11' = (att22 / `pi1_C')

	tempvar ipw_term3_0 ipw_term3_1
	gen `ipw_term3_0' = (1-att22) / `pi0_C'
	gen `ipw_term3_1' = (att22 / `pi1_C')
	
	foreach v in `ipw_term1_000' `ipw_term2_00' `ipw_term3_0' {
		centile `v' if `v'!=. & att22==0, c(1 99) 
		replace `v'=r(c_1) if `v'<r(c_1) & `v'!=. & att22==0
		replace `v'=r(c_2) if `v'>r(c_2) & `v'!=. & att22==0
		sum `v' if att22==0
		}

	foreach v in ///
		`ipw_term1_001' `ipw_term1_011' `ipw_term1_111' ///
		`ipw_term2_01' `ipw_term2_11' ///
		`ipw_term3_1' {
			centile `v' if `v'!=. & att22==1, c(1 99) 
			replace `v'=r(c_1) if `v'<r(c_1) & `v'!=. & att22==1
			replace `v'=r(c_2) if `v'>r(c_2) & `v'!=. & att22==1
			sum `v' if att22==1
			}
		
	tempvar psi000_summand
	gen `psi000_summand' = `ipw_term1_000'*(std_cesd_age40 - `mu0_CM1M2') ///
		+ `ipw_term2_00'*(`mu0_CM1M2' - `nu00_CM1') ///
		+ `ipw_term3_0'*(`nu00_CM1' - `xi000_C') ///
		+ `xi000_C'

	tempvar psi001_summand
	gen `psi001_summand' = `ipw_term1_001'*(std_cesd_age40 - `mu1_CM1M2') ///
		+ `ipw_term2_00'*(`mu1_CM1M2' - `nu01_CM1') ///
		+ `ipw_term3_0'*(`nu01_CM1' - `xi001_C') ///
		+ `xi001_C'
		
	tempvar psi011_summand
	gen `psi011_summand' = `ipw_term1_011'*(std_cesd_age40 - `mu1_CM1M2') ///
		+ `ipw_term2_01'*(`mu1_CM1M2' - `nu11_CM1') ///
		+ `ipw_term3_0'*(`nu11_CM1' - `xi011_C') ///
		+ `xi011_C'

	tempvar psi111_summand
	gen `psi111_summand' = `ipw_term1_111'*(std_cesd_age40 - `mu1_CM1M2') ///
		+ `ipw_term2_11'*(`mu1_CM1M2' - `nu11_CM1') ///
		+ `ipw_term3_1'*(`nu11_CM1' - `xi111_C') ///
		+ `xi111_C'
		
	reg `psi000_summand'
	return scalar psi000 = _b[_cons]

	reg `psi001_summand'
	return scalar psi001 = _b[_cons]
	
	reg `psi011_summand'
	return scalar psi011 = _b[_cons]
	
	reg `psi111_summand'
	return scalar psi111 = _b[_cons]
	
end mrpath

//parametric MR estimates of path-specific effects (type 2)
qui bootstrap ///
	ATE=(r(psi111)-r(psi000)) ///
	PSE_DY=(r(psi001)-r(psi000)) ///
	PSE_DM2Y=(r(psi011)-r(psi001)) ///
	PSE_DM1Y=(r(psi111)-r(psi011)), ///
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

	tempvar dvar_orig
	qui gen `dvar_orig' = att22

	tempvar u kpart
	set seed 60637
	qui gen `u' = uniform()
	qui sort `u'
	qui gen `kpart' = ceil(_n/(_N/5))

	local tvars	///
		pi0_C pi1_C ///
		pi1_CM1 pi0_CM1 ///
		pi1_CM1M2 pi0_CM1M2 ///
		mu1_CM1M2 mu0_CM1M2 ///
		nu11_CM1 nu01_CM1 nu00_CM1 ///
		xi000_C xi001_C xi011_C xi111_C
		
	tempvar `tvars'
	
	foreach v in `tvars' {
		qui gen ``v'' = .
		}

	forval k=1/5 {

		qui rforest att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

		tempvar	xxpi0_C xxpi1_C
		qui predict `xxpi0_C' `xxpi1_C' if `kpart'==`k', pr
		qui replace `pi0_C' = `xxpi0_C' if `kpart'==`k'
		qui replace `pi1_C' = `xxpi1_C' if `kpart'==`k'
			
		qui rforest att22 ever_unemp_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

		tempvar	xxpi0_CM1 xxpi1_CM1
		qui predict `xxpi0_CM1' `xxpi1_CM1' if `kpart'==`k', pr
		qui replace `pi0_CM1' = `xxpi0_CM1' if `kpart'==`k'
		qui replace `pi1_CM1' = `xxpi1_CM1' if `kpart'==`k'

		qui rforest att22 ever_unemp_age3539 log_faminc_adj_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

		tempvar	xxpi0_CM1M2 xxpi1_CM1M2
		qui predict `xxpi0_CM1M2' `xxpi1_CM1M2' if `kpart'==`k', pr
		qui replace `pi0_CM1M2' = `xxpi0_CM1M2' if `kpart'==`k'
		qui replace `pi1_CM1M2' = `xxpi1_CM1M2' if `kpart'==`k'

		qui rforest std_cesd_age40 att22 ever_unemp_age3539 log_faminc_adj_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace att22 = 1
		tempvar xxmu1_CM1M2
		qui predict `xxmu1_CM1M2'
		qui replace `mu1_CM1M2' = `xxmu1_CM1M2' if `kpart'==`k'
	
		qui replace att22 = 0
		tempvar xxmu0_CM1M2
		qui predict `xxmu0_CM1M2'
		qui replace `mu0_CM1M2' = `xxmu0_CM1M2' if `kpart'==`k'
	
		qui replace att22 = `dvar_orig'
	
		qui rforest `xxmu1_CM1M2' att22 ever_unemp_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
		qui replace att22 = 1
		tempvar xxnu11_CM1
		qui predict `xxnu11_CM1'
		qui replace `nu11_CM1' = `xxnu11_CM1' if `kpart'==`k'
	
		qui replace att22 = 0
		tempvar xxnu01_CM1
		qui predict `xxnu01_CM1'
		qui replace `nu01_CM1' = `xxnu01_CM1' if `kpart'==`k'

		qui replace att22 = `dvar_orig'

		qui rforest `xxmu0_CM1M2' att22 ever_unemp_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace att22 = 0
		tempvar xxnu00_CM1
		qui predict `xxnu00_CM1'
		qui replace `nu00_CM1' = `xxnu00_CM1' if `kpart'==`k'

		qui replace att22 = `dvar_orig'	

		qui rforest `xxnu00_CM1' att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace att22 = 0
		tempvar xxxi000_C
		qui predict `xxxi000_C'
		qui replace `xi000_C' = `xxxi000_C' if `kpart'==`k'
	
		qui replace att22 = `dvar_orig'

		qui rforest `xxnu01_CM1' att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace att22 = 0
		tempvar xxxi001_C
		qui predict `xxxi001_C'
		qui replace `xi001_C' = `xxxi001_C' if `kpart'==`k'
	
		qui replace att22 = `dvar_orig'
		
		qui rforest `xxnu11_CM1' att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3 ///
			if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
		
		qui replace att22 = 0
		tempvar xxxi011_C
		qui predict `xxxi011_C'
		qui replace `xi011_C' = `xxxi011_C' if `kpart'==`k'

		qui replace att22 = 1
		tempvar xxxi111_C
		qui predict `xxxi111_C'
		qui replace `xi111_C' = `xxxi111_C' if `kpart'==`k'
	
		qui replace att22 = `dvar_orig'
		
		qui drop ///
			`xxpi0_C' `xxpi1_C' ///
			`xxpi1_CM1' `xxpi0_CM1' ///
			`xxpi1_CM1M2' `xxpi0_CM1M2' ///
			`xxmu1_CM1M2' `xxmu0_CM1M2' ///
			`xxnu11_CM1' `xxnu01_CM1' `xxnu00_CM1' ///
			`xxxi000_C' `xxxi001_C' `xxxi011_C' `xxxi111_C'
		
		}
		
	tempvar ipw_term1_000 ipw_term1_001 ipw_term1_011 ipw_term1_111
	qui gen `ipw_term1_000' = (1-att22) / `pi0_C'
	qui gen `ipw_term1_001' = (att22 / `pi0_C') * ((`pi0_CM1M2'*`pi0_CM1') / (`pi1_CM1M2'*`pi0_CM1'))
	qui gen `ipw_term1_011' = (att22 / `pi0_C') * ((`pi1_CM1M2'*`pi0_CM1') / (`pi1_CM1M2'*`pi1_CM1'))
	qui gen `ipw_term1_111' = (att22 / `pi1_C')
	
	tempvar ipw_term2_00 ipw_term2_01 ipw_term2_11
	qui gen `ipw_term2_00' = (1-att22) / `pi0_C'
	qui gen `ipw_term2_01' = (att22 / `pi0_C') * (`pi0_CM1' / `pi1_CM1')
	qui gen `ipw_term2_11' = (att22 / `pi1_C')

	tempvar ipw_term3_0 ipw_term3_1
	qui gen `ipw_term3_0' = (1-att22) / `pi0_C'
	qui gen `ipw_term3_1' = (att22 / `pi1_C')
	
	foreach v in `ipw_term1_000' `ipw_term2_00' `ipw_term3_0' {
		qui centile `v' if `v'!=. & att22==0, c(1 99) 
		qui replace `v'=r(c_1) if `v'<r(c_1) & `v'!=. & att22==0
		qui replace `v'=r(c_2) if `v'>r(c_2) & `v'!=. & att22==0
		}

	foreach v in ///
		`ipw_term1_001' `ipw_term1_011' `ipw_term1_111' ///
		`ipw_term2_01' `ipw_term2_11' ///
		`ipw_term3_1' {
			qui centile `v' if `v'!=. & att22==1, c(1 99) 
			qui replace `v'=r(c_1) if `v'<r(c_1) & `v'!=. & att22==1
			qui replace `v'=r(c_2) if `v'>r(c_2) & `v'!=. & att22==1
			}
		
	tempvar psi000_summand
	qui gen `psi000_summand' = `ipw_term1_000'*(std_cesd_age40 - `mu0_CM1M2') ///
		+ `ipw_term2_00'*(`mu0_CM1M2' - `nu00_CM1') ///
		+ `ipw_term3_0'*(`nu00_CM1' - `xi000_C') ///
		+ `xi000_C'

	tempvar psi001_summand
	qui gen `psi001_summand' = `ipw_term1_001'*(std_cesd_age40 - `mu1_CM1M2') ///
		+ `ipw_term2_00'*(`mu1_CM1M2' - `nu01_CM1') ///
		+ `ipw_term3_0'*(`nu01_CM1' - `xi001_C') ///
		+ `xi001_C'
		
	tempvar psi011_summand
	qui gen `psi011_summand' = `ipw_term1_011'*(std_cesd_age40 - `mu1_CM1M2') ///
		+ `ipw_term2_01'*(`mu1_CM1M2' - `nu11_CM1') ///
		+ `ipw_term3_0'*(`nu11_CM1' - `xi011_C') ///
		+ `xi011_C'

	tempvar psi111_summand
	qui gen `psi111_summand' = `ipw_term1_111'*(std_cesd_age40 - `mu1_CM1M2') ///
		+ `ipw_term2_11'*(`mu1_CM1M2' - `nu11_CM1') ///
		+ `ipw_term3_1'*(`nu11_CM1' - `xi111_C') ///
		+ `xi111_C'
		
	qui gen ATE = `psi111_summand' - `psi000_summand' 
	qui gen PSE_DY = `psi001_summand' - `psi000_summand' 
	qui gen PSE_DM2Y = `psi011_summand' - `psi001_summand' 
	qui gen PSE_DM1Y = `psi111_summand' - `psi011_summand' 
	
	mean ATE PSE_DY PSE_DM2Y PSE_DM1Y, noheader
		
	qui drop ATE PSE_DY PSE_DM2Y PSE_DM1Y
		
end dmlpath

//DML estimates of path-specific effects (type 2)
dmlpath

log close
