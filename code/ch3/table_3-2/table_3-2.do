/*Table 3.2*/
capture clear all
capture log close
set more off

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

log using "${logdir}table_3-2.log", replace 

ssc install rwrmed 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==.

egen std_cesd_age40=std(cesd_age40)

///linear model estimators

//ATEhat^lm, NDEhat^lm, NIEhat^lm, CDEhat^lm
quietly rwrmed std_cesd_age40, ///
	avar(att22) mvar(ever_unemp_age3539) cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(0) mreg(regress) nointer 
		
mat list e(b)

//ATEhat^lmi, NDEhat^lmi, NIEhat^lmi, CDE0hat^lmi
quietly rwrmed std_cesd_age40, ///
	avar(att22) mvar(ever_unemp_age3539) cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(0) mreg(regress)
		
mat list e(b)

//ATEhat^lmi+, NDEhat^lmi+, NIEhat^lmi+, CDE0hat^lmi+
quietly rwrmed std_cesd_age40, ///
	avar(att22) mvar(ever_unemp_age3539) cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(0) mreg(regress) cxa cxm
		
mat list e(b)

log close
