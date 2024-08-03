/*Table 5.6*/
capture clear all
capture log close
set more off

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

log using "${logdir}table_5-6.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==.  | log_faminc_adj_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

//pure regression imputation estimates
capture program drop pathimp
program define pathimp, rclass

	gen D_orig=att22

	//compute regression imputations
	reg std_cesd_age40 i.att22 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3

	replace att22 = 0

	predict	Y0hat, xb

	replace att22 = 1

	predict Y1hat, xb

	replace att22 = D_orig

	reg std_cesd_age40 i.att22 i.ever_unemp_age3539 c.log_faminc_adj_age3539 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3

	replace att22 = 1

	predict Y1M1DM2Dhat, xb

	replace att22 = D_orig

	reg Y1M1DM2Dhat ///
		i.att22##i.female ///
		i.att22##i.black ///
		i.att22##i.hispan ///
		i.att22##c.paredu ///
		i.att22##i.parprof ///
		i.att22##c.parinc_prank ///
		i.att22##c.famsize c.afqt3

	replace att22 = 0

	predict Y1M10M200hat, xb

	replace att22 = D_orig

	reg std_cesd_age40 i.att22 i.ever_unemp_age3539 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3

	replace att22 = 1

	predict Y1M1Dhat, xb

	replace att22 = D_orig

	reg Y1M1Dhat ///
		i.att22##i.female ///
		i.att22##i.black ///
		i.att22##i.hispan ///
		i.att22##c.paredu ///
		i.att22##i.parprof ///
		i.att22##c.parinc_prank ///
		i.att22##c.famsize c.afqt3

	replace att22 = 0

	predict Y1M10M210hat, xb

	replace att22 = D_orig

	//average imputations over sample members
	reg Y0hat
	local Ehat_Y0=_b[_cons]
	drop Y0hat

	reg Y1hat
	local Ehat_Y1=_b[_cons]
	drop Y1hat

	reg Y1M10M200hat
	local Ehat_Y1M10M200=_b[_cons]
	drop Y1M10M200hat

	reg Y1M10M210hat
	local Ehat_Y1M10M210=_b[_cons]
	drop Y1M10M210hat

	//compute effect estimates
	return scalar ATE=`Ehat_Y1'-`Ehat_Y0'
	return scalar PSE_DY=`Ehat_Y1M10M200'-`Ehat_Y0'
	return scalar PSE_DM2Y=`Ehat_Y1M10M210'-`Ehat_Y1M10M200'
	return scalar PSE_DM1Y=`Ehat_Y1'-`Ehat_Y1M10M210'
		
	drop D_orig Y1M1DM2Dhat Y1M1Dhat
	
end

quietly pathimp

di r(ATE)
di r(PSE_DY)
di r(PSE_DM2Y)
di r(PSE_DM1Y)

quietly bootstrap ///
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_DM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): pathimp
		
mat list e(ci_percentile)

//regression imputation + weighting estimates
capture program drop pathimpwt
program define pathimpwt, rclass

	gen D_orig=att22
	
	//compute weights
	logit att22 i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3
	
	predict phatD_C
	
	reg att22
	
	predict phatD
	
	gen wt = (1-phatD)/(1-phatD_C) if att22==0

	//compute regression imputations
	reg std_cesd_age40 i.att22 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3

	replace att22 = 0

	predict	Y0hat, xb

	replace att22 = 1

	predict Y1hat, xb

	replace att22 = D_orig
	
	reg std_cesd_age40 i.att22 i.ever_unemp_age3539 c.log_faminc_adj_age3539 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3

	replace att22 = 1

	predict Y1M1DM2Dhat, xb

	replace att22 = D_orig

	reg std_cesd_age40 i.att22 i.ever_unemp_age3539 ///
		i.female i.black i.hispan c.paredu i.parprof c.parinc_prank c.famsize c.afqt3

	replace att22 = 1

	predict Y1M1Dhat, xb

	replace att22 = D_orig

	//average imputations over sample members
	reg Y1M1DM2Dhat [pw=wt] if att22==0
	local Ehat_Y1M10M200=_b[_cons]
	drop Y1M1DM2Dhat
	
	reg Y1M1Dhat [pw=wt] if att22==0
	local Ehat_Y1M10M210=_b[_cons]
	drop Y1M1Dhat

	reg Y0hat
	local Ehat_Y0=_b[_cons]
	drop Y0hat
	
	reg Y1hat
	local Ehat_Y1=_b[_cons]
	drop Y1hat
	
	//compute effect estimates
	return scalar ATE=`Ehat_Y1'-`Ehat_Y0'
	return scalar PSE_DY=`Ehat_Y1M10M200'-`Ehat_Y0'
	return scalar PSE_DM2Y=`Ehat_Y1M10M210'-`Ehat_Y1M10M200'
	return scalar PSE_DM1Y=`Ehat_Y1'-`Ehat_Y1M10M210'
		
	drop D_orig wt phatD_C phatD
	
end

quietly pathimpwt

di r(ATE)
di r(PSE_DY)
di r(PSE_DM2Y)
di r(PSE_DM1Y)

quietly bootstrap ///
	ATE=r(ATE) PSE_DY=r(PSE_DY) PSE_DM2Y=r(PSE_DM2Y) PSE_DM1Y=r(PSE_DM1Y), ///
	reps(2000) seed(60637): pathimpwt
		
mat list e(ci_percentile)
