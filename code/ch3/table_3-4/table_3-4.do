/*Table 3.4*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/ipwmed, replace //module to estimate natural effects
github install causalMedAnalysis/ipwcde, replace //module to estimate controlled direct effects

//specify directories 
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

//download data
copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_3-4.log", replace 

//load data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, female, black, ///
	hispan, paredu, parprof, parinc_prank, famsize, afqt3)
	
//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//define macros for different variables
global C female black hispan paredu parprof parinc_prank famsize afqt3 //baseline confounders
global D att22 //exposure
global M ever_unemp_age3539 //mediator
global Y std_cesd_age40 //outcome

//compute point estimates using IPW
qui ipwmed $Y $M, dvar($D) d(1) dstar(0) cvars($C) censor(1 99)
mat list e(b)

qui ipwcde $Y, dvar($D) mvar($M) d(1) dstar(0) m(0) cvars($C) censor(1 99)
mat list e(b)

log close

//note the ipwmed and ipwcde estimates differ slightly from those reported in
//the text, which are based on the R implementation. This is due to slight 
//differences in how the weights are censored.
