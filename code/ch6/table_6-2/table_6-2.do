/*Table 6.2*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch6\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\D\projects\causal_mediation_text\code\ch6\_LOGS\"

log using "${logdir}table_6-2.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if missing(cesd_age40, att22, ever_unemp_age3539, female, black, ///
	hispan, paredu, parprof, parinc_prank, famsize, afqt3, log_faminc_adj_age3539)

egen std_cesd_age40=std(cesd_age40)

/*****************************************************************
define parametric MR estimator for natural effects (types 1 and 2)
******************************************************************/
capture program drop mrmed
program define mrmed, rclass

	syntax , type(string)
	
	if ("`type'"=="type1") {

		tempvar dvar_orig mvar_orig
		qui gen `dvar_orig' = att22
		qui gen `mvar_orig' = ever_unemp_age3539
	
		logit att22 female black hispan paredu parprof parinc_prank famsize afqt3
	
		tempvar phat_D1_C
		predict `phat_D1_C', pr
	
		logit ever_unemp_age3539 att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3
	
		replace att22 = 0
		tempvar f_M1_CD0 f_M0_CD0 f_M_CD0
		predict `f_M1_CD0', pr
		qui gen `f_M0_CD0' = 1 - `f_M1_CD0'
		qui gen `f_M_CD0' = `f_M1_CD0'*ever_unemp_age3539 + `f_M0_CD0'*(1-ever_unemp_age3539)

		replace att22 = 1
		tempvar f_M1_CD1 f_M0_CD1 f_M_CD1
		predict `f_M1_CD1', pr
		qui gen `f_M0_CD1' = 1 - `f_M1_CD1'	
		qui gen `f_M_CD1' = `f_M1_CD1'*ever_unemp_age3539 + `f_M0_CD1'*(1-ever_unemp_age3539)
	
		replace att22 = `dvar_orig'

		reg std_cesd_age40 i.ever_unemp_age3539##i.att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3
	
		replace att22 = 0
		tempvar mu0_CM
		predict `mu0_CM'
	
		replace att22 = 1
		tempvar mu1_CM
		predict `mu1_CM'
	
		replace att22 = 0
		replace ever_unemp_age3539 = 0
		tempvar mu0_CM0
		predict `mu0_CM0'
	
		replace att22 = 0
		replace ever_unemp_age3539 = 1
		tempvar mu0_CM1
		predict `mu0_CM1'

		replace att22 = 1
		replace ever_unemp_age3539 = 0
		tempvar mu1_CM0
		predict `mu1_CM0'
	
		replace att22 = 1
		replace ever_unemp_age3539 = 1
		tempvar mu1_CM1
		predict `mu1_CM1'
	
		replace att22 = `dvar_orig'
		replace ever_unemp_age3539 = `mvar_orig'
	
		tempvar ipw0 ipw1 rmpw
		qui gen `ipw0' = (1-att22)/(1-`phat_D1_C')
		qui gen `ipw1' = att22/`phat_D1_C'
		qui gen `rmpw' = (`f_M_CD0'/`f_M_CD1')

		qui centile `ipw0' if `ipw0'!=. & att22==0, c(1 99) 
		qui replace `ipw0'=r(c_1) if `ipw0'<r(c_1) & `ipw0'!=. & att22==0
		qui replace `ipw0'=r(c_2) if `ipw0'>r(c_2) & `ipw0'!=. & att22==0
	
		qui centile `ipw1' if `ipw1'!=. & att22==1, c(1 99) 
		qui replace `ipw1'=r(c_1) if `ipw1'<r(c_1) & `ipw1'!=. & att22==1
		qui replace `ipw1'=r(c_2) if `ipw1'>r(c_2) & `ipw1'!=. & att22==1
	
		qui centile `rmpw' if `rmpw'!=., c(1 99) 
		qui replace `rmpw'=r(c_1) if `rmpw'<r(c_1) & `rmpw'!=. 
		qui replace `rmpw'=r(c_2) if `rmpw'>r(c_2) & `rmpw'!=. 
	
		tempvar dr11_summand
		qui gen `dr11_summand' = `ipw1'*(std_cesd_age40 - `mu1_CM') ///
			+ `ipw1'*(`mu1_CM' - (`mu1_CM0'*`f_M0_CD1' + `mu1_CM1'*`f_M1_CD1')) ///
			+ (`mu1_CM0'*`f_M0_CD1' + `mu1_CM1'*`f_M1_CD1')
	
		tempvar dr00_summand
		qui gen `dr00_summand' = `ipw0'*(std_cesd_age40 - `mu0_CM') ///
			+ `ipw0'*(`mu0_CM' - (`mu0_CM0'*`f_M0_CD0' + `mu0_CM1'*`f_M1_CD0')) ///
			+ (`mu0_CM0'*`f_M0_CD0' + `mu0_CM1'*`f_M1_CD0')

		tempvar dr01_summand
		qui gen `dr01_summand' = `ipw1'*`rmpw'*(std_cesd_age40 - `mu1_CM') ///
			+ `ipw0'*(`mu1_CM' - (`mu1_CM0'*`f_M0_CD0' + `mu1_CM1'*`f_M1_CD0')) ///
			+ (`mu1_CM0'*`f_M0_CD0' + `mu1_CM1'*`f_M1_CD0')
			
		reg `dr11_summand'
		return scalar psi11 = _b[_cons]
	
		reg `dr00_summand'
		return scalar psi00 = _b[_cons]

		reg `dr01_summand'
		return scalar psi01 = _b[_cons]
		}
	
	if ("`type'"=="type2") {
		
		tempvar dvar_orig mvar_orig
		qui gen `dvar_orig' = att22
	
		logit att22 female black hispan paredu parprof parinc_prank famsize afqt3
	
		tempvar pi0_C pi1_C
		predict `pi1_C', pr
		qui gen `pi0_C' = 1 - `pi1_C'
		
		logit att22 ever_unemp_age3539 ///
			female black hispan paredu parprof parinc_prank famsize afqt3
	
		tempvar pi1_CM pi0_CM
		predict `pi1_CM', pr
		qui gen `pi0_CM' = 1 - `pi1_CM'
		
		reg std_cesd_age40 i.ever_unemp_age3539##i.att22 ///
			female black hispan paredu parprof parinc_prank famsize afqt3
	
		replace att22 = 0
		tempvar mu0_CM
		predict `mu0_CM'
	
		replace att22 = 1
		tempvar mu1_CM
		predict `mu1_CM'
	
		replace att22 = `dvar_orig'

		reg `mu1_CM' i.att22##(i.female i.black i.hispan ///
			c.paredu i.parprof c.parinc_prank c.famsize c.afqt3)
	
		replace att22 = 1
		tempvar nu1Ofmu1_C
		predict `nu1Ofmu1_C'

		replace att22 = 0
		tempvar nu0Ofmu1_C
		predict `nu0Ofmu1_C'

		replace att22 = `dvar_orig'

		reg `mu0_CM' i.att22##(i.female i.black i.hispan ///
			c.paredu i.parprof c.parinc_prank c.famsize c.afqt3)
	
		replace att22 = 0
		tempvar nu0Ofmu0_C
		predict `nu0Ofmu0_C'

		replace att22 = `dvar_orig'
	
		tempvar ipw0C ipw1C ipw01CM
		qui gen `ipw0C' = (1-att22)/`pi0_C'
		qui gen `ipw1C' = att22/`pi1_C'
		qui gen `ipw01CM' = (att22/`pi0_C')*(`pi0_CM'/`pi1_CM')
		
		qui centile `ipw0C' if `ipw0C'!=. & att22==0, c(1 99) 
		qui replace `ipw0C'=r(c_1) if `ipw0C'<r(c_1) & `ipw0C'!=. & att22==0
		qui replace `ipw0C'=r(c_2) if `ipw0C'>r(c_2) & `ipw0C'!=. & att22==0
	
		qui centile `ipw1C' if `ipw1C'!=. & att22==1, c(1 99) 
		qui replace `ipw1C'=r(c_1) if `ipw1C'<r(c_1) & `ipw1C'!=. & att22==1
		qui replace `ipw1C'=r(c_2) if `ipw1C'>r(c_2) & `ipw1C'!=. & att22==1
	
		qui centile `ipw01CM' if `ipw01CM'!=. & att22==1, c(1 99) 
		qui replace `ipw01CM'=r(c_1) if `ipw01CM'<r(c_1) & `ipw01CM'!=. & att22==1
		qui replace `ipw01CM'=r(c_2) if `ipw01CM'>r(c_2) & `ipw01CM'!=. & att22==1
		
		tempvar dr11_summand
		qui gen `dr11_summand' = `ipw1C'*(std_cesd_age40 - `mu1_CM') ///
			+ `ipw1C'*(`mu1_CM' - `nu1Ofmu1_C') ///
			+ `nu1Ofmu1_C'
		
		tempvar dr00_summand
		qui gen `dr00_summand' = `ipw0C'*(std_cesd_age40 - `mu0_CM') ///
			+ `ipw0C'*(`mu0_CM' - `nu0Ofmu0_C') ///
			+ `nu0Ofmu0_C'

		tempvar dr01_summand
		qui gen `dr01_summand' = `ipw01CM'*(std_cesd_age40 - `mu1_CM') ///
			+ `ipw0C'*(`mu1_CM' - `nu0Ofmu1_C') ///
			+ `nu0Ofmu1_C'
			
		reg `dr11_summand'
		return scalar psi11 = _b[_cons]
	
		reg `dr00_summand'
		return scalar psi00 = _b[_cons]

		reg `dr01_summand'
		return scalar psi01 = _b[_cons]
		}
	
end mrmed

//parametric MR estimates of natural effects (type 1)
qui bootstrap ///
	ATE=(r(psi11)-r(psi00)) ///
	NDE=(r(psi01)-r(psi00)) ///
	NIE=(r(psi11)-r(psi01)), ///
	reps(2000) seed(60637): mrmed, type(type1)

estat bootstrap, p noheader

//parametric MR estimates of natural effects (type 2)

qui bootstrap ///
	ATE=(r(psi11)-r(psi00)) ///
	NDE=(r(psi01)-r(psi00)) ///
	NIE=(r(psi11)-r(psi01)), ///
	reps(2000) seed(60637): mrmed, type(type2)

estat bootstrap, p noheader

/*******************************************************
define DML estimator for natural effects (types 1 and 2)
********************************************************/
//note that stata does not currently support use of a superLearner
//we therefore implement the DML estimator using only random forests

*ssc install rforest //install rforest module if not already installed

capture program drop dmlmed
program define dmlmed, rclass

	syntax , type(string)
	
	if ("`type'"=="type1") {
		tempvar dvar_orig mvar_orig
		qui gen `dvar_orig' = att22
		qui gen `mvar_orig' = ever_unemp_age3539
	
		tempvar kpart
		set seed 60637
		qui gen u = uniform()
		qui sort u
		qui gen `kpart' = ceil(_n/(_N/5))
		drop u

		local tvars	phat_D0_C phat_D1_C ///
			f_M1_CD0 f_M0_CD0 f_M_CD0 f_M1_CD1 f_M0_CD1 f_M_CD1 ///
			mu0_CM mu1_CM mu0_CM0 mu0_CM1 mu1_CM0 mu1_CM1
		
		tempvar `tvars'
	
		foreach v in `tvars' {
			qui gen ``v'' = .
			}

		forval k=1/5 {

			qui rforest att22 ///
				female black hispan paredu parprof parinc_prank famsize afqt3 ///
				if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

			tempvar	xxphat_D0_C xxphat_D1_C
			qui predict `xxphat_D0_C' `xxphat_D1_C' if `kpart'==`k', pr
			qui replace `phat_D0_C' = `xxphat_D0_C' if `kpart'==`k'
			qui replace `phat_D1_C' = `xxphat_D1_C' if `kpart'==`k'
		
			qui rforest ever_unemp_age3539 att22 ///
				female black hispan paredu parprof parinc_prank famsize afqt3 ///
				if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)		

			qui replace att22 = 0
			tempvar	xxf_M0_CD0 xxf_M1_CD0
			qui predict `xxf_M0_CD0' `xxf_M1_CD0' if `kpart'==`k', pr
			qui replace `f_M0_CD0' = `xxf_M0_CD0' if `kpart'==`k'
			qui replace `f_M1_CD0' = `xxf_M1_CD0' if `kpart'==`k'
			qui replace `f_M_CD0' = `f_M1_CD0'*ever_unemp_age3539 + `f_M0_CD0'*(1-ever_unemp_age3539) if `kpart'==`k'

			qui replace att22 = 1
			tempvar	xxf_M0_CD1 xxf_M1_CD1
			qui predict `xxf_M0_CD1' `xxf_M1_CD1' if `kpart'==`k', pr
			qui replace `f_M0_CD1' = `xxf_M0_CD1' if `kpart'==`k'
			qui replace `f_M1_CD1' = `xxf_M1_CD1' if `kpart'==`k'
			qui replace `f_M_CD1' = `f_M1_CD1'*ever_unemp_age3539 + `f_M0_CD1'*(1-ever_unemp_age3539) if `kpart'==`k'
		
			qui replace att22 = `dvar_orig'

			qui rforest std_cesd_age40 att22 ever_unemp_age3539 ///
				female black hispan paredu parprof parinc_prank famsize afqt3 ///
				if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
			qui replace att22 = 0
			tempvar xxmu0_CM
			qui predict `xxmu0_CM' if `kpart'==`k'
			qui replace `mu0_CM' = `xxmu0_CM' if `kpart'==`k'
	
			qui replace att22 = 1
			tempvar xxmu1_CM
			qui predict `xxmu1_CM' if `kpart'==`k'
			qui replace `mu1_CM' = `xxmu1_CM' if `kpart'==`k'

			qui replace att22 = 0
			qui replace ever_unemp_age3539 = 0
			tempvar xxmu0_CM0
			qui predict `xxmu0_CM0' if `kpart'==`k'
			qui replace `mu0_CM0' = `xxmu0_CM0' if `kpart'==`k'
	
			qui replace att22 = 0
			qui replace ever_unemp_age3539 = 1
			tempvar xxmu0_CM1
			qui predict `xxmu0_CM1' if `kpart'==`k'
			qui replace `mu0_CM1' = `xxmu0_CM1' if `kpart'==`k'

			qui replace att22 = 1
			qui replace ever_unemp_age3539 = 0
			tempvar xxmu1_CM0
			qui predict `xxmu1_CM0' if `kpart'==`k'
			qui replace `mu1_CM0' = `xxmu1_CM0' if `kpart'==`k'
		
			qui replace att22 = 1
			qui replace ever_unemp_age3539 = 1
			tempvar xxmu1_CM1
			qui predict `xxmu1_CM1' if `kpart'==`k'
			qui replace `mu1_CM1' = `xxmu1_CM1' if `kpart'==`k'
	
			qui replace att22 = `dvar_orig'
			qui replace ever_unemp_age3539 = `mvar_orig'
			}
	
		tempvar ipw0 ipw1 rmpw
		qui gen `ipw0' = (1-att22)/(1-`phat_D1_C')
		qui gen `ipw1' = att22/`phat_D1_C'
		qui gen `rmpw' = att22*(`f_M_CD0'/`f_M_CD1')

		qui centile `ipw0' if `ipw0'!=. & att22==0, c(1 99) 
		qui replace `ipw0'=r(c_1) if `ipw0'<r(c_1) & `ipw0'!=. & att22==0
		qui replace `ipw0'=r(c_2) if `ipw0'>r(c_2) & `ipw0'!=. & att22==0
	
		qui centile `ipw1' if `ipw1'!=. & att22==1, c(1 99) 
		qui replace `ipw1'=r(c_1) if `ipw1'<r(c_1) & `ipw1'!=. & att22==1
		qui replace `ipw1'=r(c_2) if `ipw1'>r(c_2) & `ipw1'!=. & att22==1
	
		qui centile `rmpw' if `rmpw'!=. & att22==1, c(1 99) 
		qui replace `rmpw'=r(c_1) if `rmpw'<r(c_1) & `rmpw'!=. & att22==1
		qui replace `rmpw'=r(c_2) if `rmpw'>r(c_2) & `rmpw'!=. & att22==1
	
		tempvar dr11_summand
		qui gen `dr11_summand' = `ipw1'*(std_cesd_age40 - `mu1_CM') ///
			+ `ipw1'*(`mu1_CM' - (`mu1_CM0'*`f_M0_CD1' + `mu1_CM1'*`f_M1_CD1')) ///
			+ (`mu1_CM0'*`f_M0_CD1' + `mu1_CM1'*`f_M1_CD1')
	
		tempvar dr00_summand
		qui gen `dr00_summand' = `ipw0'*(std_cesd_age40 - `mu0_CM') ///
			+ `ipw0'*(`mu0_CM' - (`mu0_CM0'*`f_M0_CD0' + `mu0_CM1'*`f_M1_CD0')) ///
			+ (`mu0_CM0'*`f_M0_CD0' + `mu0_CM1'*`f_M1_CD0')

		tempvar dr01_summand
		qui gen `dr01_summand' = `ipw1'*`rmpw'*(std_cesd_age40 - `mu1_CM') ///
			+ `ipw0'*(`mu1_CM' - (`mu1_CM0'*`f_M0_CD0' + `mu1_CM1'*`f_M1_CD0')) ///
			+ (`mu1_CM0'*`f_M0_CD0' + `mu1_CM1'*`f_M1_CD0')
	
		qui gen ATE = `dr11_summand' - `dr00_summand'
		qui gen NDE = `dr01_summand' - `dr00_summand'
		qui gen NIE = `dr11_summand' - `dr01_summand'
		
		mean ATE NDE NIE, noheader
		
		drop ATE NDE NIE
		}
	
	if ("`type'"=="type2") {
		tempvar dvar_orig mvar_orig
		qui gen `dvar_orig' = att22

		tempvar kpart
		set seed 60637
		qui gen u = uniform()
		qui sort u
		qui gen `kpart' = ceil(_n/(_N/5))
		drop u

		local tvars	pi0_C pi1_C pi1_CM pi0_CM ///
			mu0_CM mu1_CM nu1Ofmu1_C nu0Ofmu1_C nu0Ofmu0_C
		
		tempvar `tvars'
	
		foreach v in `tvars' {
			qui gen ``v'' = .
			}

		forval k=1/5 {

			qui rforest att22 ///
				female black hispan paredu parprof parinc_prank famsize afqt3 ///
				if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

			tempvar	xxphat_D0_C xxphat_D1_C
			qui predict `xxphat_D0_C' `xxphat_D1_C' if `kpart'==`k', pr
			qui replace `pi0_C' = `xxphat_D0_C' if `kpart'==`k'
			qui replace `pi1_C' = `xxphat_D1_C' if `kpart'==`k'

			qui rforest att22 ever_unemp_age3539 ///
				female black hispan paredu parprof parinc_prank famsize afqt3 ///
				if `kpart'!=`k', type(class) iter(200) lsize(20) seed(60637)

			tempvar	xxphat_D0_CM xxphat_D1_CM
			qui predict `xxphat_D0_CM' `xxphat_D1_CM' if `kpart'==`k', pr
			qui replace `pi0_CM' = `xxphat_D0_CM' if `kpart'==`k'
			qui replace `pi1_CM' = `xxphat_D1_CM' if `kpart'==`k'
		
			qui rforest std_cesd_age40 att22 ever_unemp_age3539 ///
				female black hispan paredu parprof parinc_prank famsize afqt3 ///
				if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
			qui replace att22 = 0
			tempvar xxmu0_CM
			qui predict `xxmu0_CM'
			qui replace `mu0_CM' = `xxmu0_CM' if `kpart'==`k'
	
			qui replace att22 = 1
			tempvar xxmu1_CM
			qui predict `xxmu1_CM'
			qui replace `mu1_CM' = `xxmu1_CM' if `kpart'==`k'
		
			qui replace att22 = `dvar_orig'

			qui rforest `xxmu1_CM' att22 ///
				female black hispan paredu parprof parinc_prank famsize afqt3 ///
				if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
			qui replace att22 = 1
			tempvar xxnu1Ofmu1_C
			qui predict `xxnu1Ofmu1_C' if `kpart'==`k'
			qui replace `nu1Ofmu1_C' = `xxnu1Ofmu1_C' if `kpart'==`k'

			qui replace att22 = 0
			tempvar xxnu0Ofmu1_C
			qui predict `xxnu0Ofmu1_C' if `kpart'==`k'
			qui replace `nu0Ofmu1_C' = `xxnu0Ofmu1_C' if `kpart'==`k'
		
			qui replace att22 = `dvar_orig'

			qui rforest `xxmu0_CM' att22 ///
				female black hispan paredu parprof parinc_prank famsize afqt3 ///
				if `kpart'!=`k', type(reg) iter(200) lsize(20) seed(60637)
	
			qui replace att22 = 0
			tempvar xxnu0Ofmu0_C
			qui predict `xxnu0Ofmu0_C' if `kpart'==`k'
			qui replace `nu0Ofmu0_C' = `xxnu0Ofmu0_C' if `kpart'==`k'
		
			qui replace att22 = `dvar_orig'
			}

		tempvar ipw0C ipw1C ipw01CM
		qui gen `ipw0C' = (1-att22)/`pi0_C'
		qui gen `ipw1C' = att22/`pi1_C'
		qui gen `ipw01CM' = (att22/`pi0_C')*(`pi0_CM'/`pi1_CM')
		
		qui centile `ipw0C' if `ipw0C'!=. & att22==0, c(1 99) 
		qui replace `ipw0C'=r(c_1) if `ipw0C'<r(c_1) & `ipw0C'!=. & att22==0
		qui replace `ipw0C'=r(c_2) if `ipw0C'>r(c_2) & `ipw0C'!=. & att22==0
	
		qui centile `ipw1C' if `ipw1C'!=. & att22==1, c(1 99) 
		qui replace `ipw1C'=r(c_1) if `ipw1C'<r(c_1) & `ipw1C'!=. & att22==1
		qui replace `ipw1C'=r(c_2) if `ipw1C'>r(c_2) & `ipw1C'!=. & att22==1
	
		qui centile `ipw01CM' if `ipw01CM'!=. & att22==1, c(1 99) 
		qui replace `ipw01CM'=r(c_1) if `ipw01CM'<r(c_1) & `ipw1C'!=. & att22==1
		qui replace `ipw01CM'=r(c_2) if `ipw01CM'>r(c_2) & `ipw1C'!=. & att22==1
		
		tempvar dr11_summand
		qui gen `dr11_summand' = `ipw1C'*(std_cesd_age40 - `mu1_CM') ///
			+ `ipw1C'*(`mu1_CM' - `nu1Ofmu1_C') ///
			+ `nu1Ofmu1_C'
		
		tempvar dr00_summand
		qui gen `dr00_summand' = `ipw0C'*(std_cesd_age40 - `mu0_CM') ///
			+ `ipw0C'*(`mu0_CM' - `nu0Ofmu0_C') ///
			+ `nu0Ofmu0_C'

		tempvar dr01_summand
		qui gen `dr01_summand' = `ipw01CM'*(std_cesd_age40 - `mu1_CM') ///
			+ `ipw0C'*(`mu1_CM' - `nu0Ofmu1_C') ///
			+ `nu0Ofmu1_C'

		qui gen ATE = `dr11_summand' - `dr00_summand'
		qui gen NDE = `dr01_summand' - `dr00_summand'
		qui gen NIE = `dr11_summand' - `dr01_summand'
		
		mean ATE NDE NIE, noheader
		
		drop ATE NDE NIE
		
		}
	
end dmlmed

//DML estimates of natural effects (type 1)
dmlmed, type(type1)

//DML estimates of natural effects (type 2)
dmlmed, type(type2)

log close
