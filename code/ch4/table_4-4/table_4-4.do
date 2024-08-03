/*Table 4.4*/
capture clear all
capture log close
set more off

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

log using "${logdir}table_4-4.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==. | faminc_adj_age3539==. | ///
	log_faminc_adj_age3539==.

egen std_cesd_age40=std(cesd_age40)

///inverse probability weighting (ipw) estimators

//define function for ipw estimator for OE, IDE, IIE, and CDE
capture program drop ipwmed
program define ipwmed, rclass 
	syntax , d_model(string) l_model(string) m_model(string)

	gen D_orig=att22
	gen L_orig=ever_unemp_age3539
	gen M_orig=log_faminc_adj_age3539
	
	//predict D probabilities
	est restore `d_model'
	
	predict phatD1_C, pr
	gen phatD0_C=1-phatD1_C
	
	//predict L probabilities
	est restore `l_model'
	
	replace att22=1
	
	predict phatL1_D1C, pr
	gen phatL0_D1C=1-phatL1_D1C
	
	replace att22=0

	predict phatL1_D0C, pr
	gen phatL0_D0C=1-phatL1_D0C

	//predict M probabilities
	est restore `m_model'
	
	replace att22=1
	
	predict mhat_D1LC, xb
	gen phatM_D1LC=normalden(log_faminc_adj_age3539, mhat_D1LC, e(rmse))
	
	replace att22=0
	
	predict mhat_D0LC, xb
	gen phatM_D0LC=normalden(log_faminc_adj_age3539, mhat_D0LC, e(rmse))
	
	replace ever_unemp_age3539=0
	
	predict mhat_D0L0C, xb
	gen phatM_D0L0C=normalden(log_faminc_adj_age3539, mhat_D0L0C, e(rmse))

	replace ever_unemp_age3539=1
	
	predict mhat_D0L1C, xb
	gen phatM_D0L1C=normalden(log_faminc_adj_age3539, mhat_D0L1C, e(rmse))
	
	replace att22=1
	
	predict mhat_D1L1C, xb
	gen phatM_D1L1C=normalden(log_faminc_adj_age3539, mhat_D1L1C, e(rmse))
	
	replace ever_unemp_age3539=0
	
	predict mhat_D1L0C, xb
	gen phatM_D1L0C=normalden(log_faminc_adj_age3539, mhat_D1L0C, e(rmse))

	replace att22=D_orig
	replace ever_unemp_age3539=L_orig
	replace log_faminc_adj_age3539=M_orig
	
	//compute stabilized weights
	logit att22
	predict phatD1, pr
	gen phatD0=1-phatD1
	
	reg log_faminc_adj_age3539 att22 
	
	replace att22=0
	
	predict mhat_D0, xb
	gen phatM_D0=normalden(log_faminc_adj_age3539, mhat_D0, e(rmse))

	replace att22=1
	
	predict mhat_D1, xb
	gen phatM_D1=normalden(log_faminc_adj_age3539, mhat_D1, e(rmse))
	
	replace att22=D_orig
	replace ever_unemp_age3539=L_orig
	replace log_faminc_adj_age3539=M_orig

	gen _sw1 = (phatD0/phatD0_C) * (1/phatM_D0LC) * ((phatM_D0L0C * phatL0_D0C) + (phatM_D0L1C * phatL1_D0C)) if att22==0 
	gen _sw2 = (phatD1/phatD1_C) * (1/phatM_D1LC) * ((phatM_D1L0C * phatL0_D1C) + (phatM_D1L1C * phatL1_D1C)) if att22==1
	gen _sw3 = (phatD1/phatD1_C) * (1/phatM_D1LC) * ((phatM_D0L0C * phatL0_D0C) + (phatM_D0L1C * phatL1_D0C)) if att22==1
    
	gen _sw4 = (phatD1/phatD1_C) * (phatM_D1/phatM_D1LC) if att22==1
	replace _sw4 = (phatD0/phatD0_C) * (phatM_D0/phatM_D0LC) if att22==0
	
	//censor stabilized weights at 1st and 99th percentiles
	foreach i of var _sw* {
		centile `i', c(1 99) 
		replace `i'=r(c_1) if `i'<r(c_1) & `i'!=.
		replace `i'=r(c_2) if `i'>r(c_2) & `i'!=.
		}

	//compute weighted means of Y
	reg std_cesd_age40 [pw=_sw1] if att22==0
	local Ehat_Y0M0=_b[_cons]

	reg std_cesd_age40 [pw=_sw2] if att22==1
	local Ehat_Y1M1=_b[_cons]

	reg std_cesd_age40 [pw=_sw3] if att22==1
	local Ehat_Y1M0=_b[_cons]

	//compute interventional effect estimates
	return scalar OE = `Ehat_Y1M1' - `Ehat_Y0M0'
	return scalar IDE = `Ehat_Y1M0' - `Ehat_Y0M0'
	return scalar IIE = `Ehat_Y1M1' - `Ehat_Y1M0'
	
	//compute CDE estimate
	reg std_cesd_age40 i.att22##c.log_faminc_adj_age3539 [pw=_sw4]
	
	return scalar CDE = _b[1.att22] + ln(50000)*_b[1.att22#c.log_faminc_adj_age3539]
	
	drop phat* mhat* _sw* _est* D_orig L_orig M_orig
	
end

//compute ipw estimates w/ additive models
quietly logit att22 female black hispan paredu parprof parinc_prank famsize afqt3
est store DmodAdd
	
quietly logit ever_unemp_age3539 att22 female black hispan paredu parprof parinc_prank famsize afqt3
est store LmodAdd

quietly reg log_faminc_adj_age3539 att22 ever_unemp_age3539 female black hispan paredu parprof parinc_prank famsize afqt3
est store MmodAdd
	
quietly ipwmed, d_model(DmodAdd) l_model(LmodAdd) m_model(MmodAdd)

di r(OE)
di r(IDE)
di r(IIE)
di r(CDE)

//compute ipw estimates w/ interactive models
quietly logit att22 female black hispan paredu parprof parinc_prank famsize afqt3
est store DmodX

quietly logit ever_unemp_age3539 ///
	i.att22##i.female ///
	i.att22##i.black ///
	i.att22##i.hispan ///
	i.att22##c.paredu ///
	i.att22##i.parprof ///
	i.att22##c.parinc_prank ///
	i.att22##c.famsize ///
	i.att22##c.afqt3
est store LmodX

quietly reg log_faminc_adj_age3539 ///
	i.att22##i.female ///
	i.att22##i.black ///
	i.att22##i.hispan ///
	i.att22##c.paredu ///
	i.att22##i.parprof ///
	i.att22##c.parinc_prank ///
	i.att22##c.famsize ///
	i.att22##c.afqt3 ///
	i.att22##i.ever_unemp_age3539
est store MmodX

quietly ipwmed, d_model(DmodX) l_model(LmodX) m_model(MmodX)

di r(OE)
di r(IDE)
di r(IIE)
di r(CDE)

log close
