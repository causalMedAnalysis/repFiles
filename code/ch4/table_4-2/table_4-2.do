/*Table 4.2*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/cmed //module to perform causal mediation analysis

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch4\_LOGS\"

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_4-2.log", replace 

//load data
use "${datadir}NLSY79\nlsy79BK_ed2.dta"

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, log_faminc_adj_age3539, ///
	female, black, hispan, paredu, parprof, parinc_prank, famsize, afqt3)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//define macros for different variables
global C female black hispan paredu parprof parinc_prank famsize afqt3 //baseline confounders
global D att22 //exposure
global L ever_unemp_age3539 //exposure-induced confounder
global M log_faminc_adj_age3539 //mediator
global Y std_cesd_age40 //outcome

//compute RWR point esimates based on income in levels
qui cmed linear $Y faminc_adj_age3539 ($L) $D = $C
mat list e(b)

qui cmed linear $Y faminc_adj_age3539 ($L) $D = $C, m(50000)
mat list e(b)

//compute RWR point estimates based on log(income)
qui cmed linear $Y $M ($L) $D = $C
mat list e(b)

qui cmed linear $Y $M ($L) $D = $C, m(10.82)
mat list e(b)

//compute RWR point estimates from models with DxM, CxD, and CxM interactions
qui cmed linear $Y $M ($L) $D = $C, cxd cxm lxm
mat list e(b)

qui cmed linear $Y $M ($L) $D = $C, m(10.82) cxd cxm lxm
mat list e(b)

log close
