/*Table 4.2*/
capture clear all
capture log close
set more off

//office
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

//home
*global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

log using "${logdir}table_4-2.log", replace 

ssc install rwrmed 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

drop if ///
	cesd_age40==. | att22==. | ever_unemp_age3539==. | ///
	female==. | black==. | hispan==. | paredu==. | parprof==. | ///
	parinc_prank==. | famsize==. | afqt3==. | faminc_adj_age3539==. | ///
	log_faminc_adj_age3539==.

egen std_cesd_age40=std(cesd_age40)

///RWR estimators

//Note that rNDE = IDE, rNIE = IIE, and rATE = OE

//OEhat^rwr, IDEhat^rwr, IIEhat^rwr, CDEhat^rwr
quietly rwrmed std_cesd_age40 ever_unemp_age3539, ///
	avar(att22) mvar(faminc_adj_age3539) cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(50000) mreg(regress)  
		
mat list e(b)

//OEhat^rwr, IDEhat^rwr, IIEhat^rwr, CDEhat^rwr w/ log(income)
quietly rwrmed std_cesd_age40 ever_unemp_age3539, ///
	avar(att22) mvar(log_faminc_adj_age3539) cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(10.82) mreg(regress)
		
mat list e(b)

//OEhat^rwr+, IDEhat^rwr+, IIEhat^rwr+, CDEhat^rwr+ w/ log(income)
quietly rwrmed std_cesd_age40 ever_unemp_age3539, ///
	avar(att22) mvar(log_faminc_adj_age3539) cvar(female black hispan paredu parprof parinc_prank famsize afqt3) ///
	a(0) astar(1) m(10.82) mreg(regress) cxa cxm lxm
		
mat list e(b)

log close
