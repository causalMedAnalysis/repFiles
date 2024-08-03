/*Table 5.1*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

log using "${logdir}table_5-1.log", replace 

*ssc install rwrmed 

///linear model estimators for one-mediator-at-a-time approach

//unemployment as mediator

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

//additive specification
quietly rwrmed std_cesd_age40, ///
	avar(att22) mvar(ever_unemp_age3539) ///
	cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(0) mreg(regress) nointer 
		
di _b[ATE]
di _b[NDE]
di _b[NIE]

quietly bootstrap ///
	ATE=_b[ATE] NDE=_b[NDE] NIE=_b[NIE], ///
	reps(2000) seed(60637): ///
		rwrmed std_cesd_age40, ///
		avar(att22) mvar(ever_unemp_age3539) ///
		cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
		a(0) astar(1) m(0) mreg(regress) nointer 
		
mat list e(ci_percentile)


//specification with an exposure-mediator interaction
quietly rwrmed std_cesd_age40, ///
	avar(att22) mvar(ever_unemp_age3539) ///
	cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(0) mreg(regress) 
		
di _b[ATE]
di _b[NDE]
di _b[NIE]

quietly bootstrap ///
	ATE=_b[ATE] NDE=_b[NDE] NIE=_b[NIE], ///
	reps(2000) seed(60637): ///
		rwrmed std_cesd_age40, ///
		avar(att22) mvar(ever_unemp_age3539) ///
		cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
		a(0) astar(1) m(0) mreg(regress) 
		
mat list e(ci_percentile)

//income as mediator

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | log_faminc_adj_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

//additive specification
quietly rwrmed std_cesd_age40, ///
	avar(att22) mvar(log_faminc_adj_age3539) ///
	cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(0) mreg(regress) nointer 
		
di _b[ATE]
di _b[NDE]
di _b[NIE]

quietly bootstrap ///
	ATE=_b[ATE] NDE=_b[NDE] NIE=_b[NIE], ///
	reps(2000) seed(60637): ///
		rwrmed std_cesd_age40, ///
		avar(att22) mvar(log_faminc_adj_age3539) ///
		cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
		a(0) astar(1) m(0) mreg(regress) nointer 
		
mat list e(ci_percentile)

//specification with an exposure-mediator interaction
quietly rwrmed std_cesd_age40, ///
	avar(att22) mvar(log_faminc_adj_age3539) ///
	cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(0) mreg(regress) 
		
di _b[ATE]
di _b[NDE]
di _b[NIE]

quietly bootstrap ///
	ATE=_b[ATE] NDE=_b[NDE] NIE=_b[NIE], ///
	reps(2000) seed(60637): ///
		rwrmed std_cesd_age40, ///
		avar(att22) mvar(log_faminc_adj_age3539) ///
		cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
		a(0) astar(1) m(0) mreg(regress) 
		
mat list e(ci_percentile)

log close
