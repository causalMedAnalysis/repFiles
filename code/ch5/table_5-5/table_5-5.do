/*Table 5.5*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/linpath, replace //module to estimate PSEs
github install causalMedAnalysis/ipwpath, replace //module to estimate PSEs

//specify directories 
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//download data
copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_5-5.log", replace 

//load data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, log_faminc_adj_age3539, ///
	female, black, hispan, paredu, parprof, parinc_prank, famsize, afqt3)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//define macros for different variables
global C female black hispan paredu parprof parinc_prank famsize afqt3 //baseline confounders
global D att22 //exposure
global M1 ever_unemp_age3539 //first mediator
global M2 log_faminc_adj_age3539 //second mediator
global Y std_cesd_age40 //outcome

//compute PSE estimates from additive linear models
qui linpath $Y $M1 $M2, dvar($D) d(1) dstar(0) cvars($C) nointer reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)
	
//compute PSE estimates from linear model estimates w/ DxM interaction
qui linpath $Y $M1 $M2, dvar($D) d(1) dstar(0) cvars($C) reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//compute PSE estimates using inverse probability weighting
qui ipwpath $Y $M1 $M2, dvar($D) d(1) dstar(0) cvars($C) censor(1 99) reps(2000) seed(60637) 

mat list e(b)
mat list e(ci_percentile)

log close

//note the ipwpath estimates differ slightly from those reported in
//the text, which are based on the R implementation. This is due to minor
//differences in how the weights are censored.
