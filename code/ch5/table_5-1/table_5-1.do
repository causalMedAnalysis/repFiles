/*Table 5.1*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/linmed, replace //module to estimate natural effects

//specify directories 
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch5\_LOGS\"

//download data
copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_5-1.log", replace 

//load data
use "${datadir}NLSY79\nlsy79BK_ed2.dta", clear

//define macros for different variables
global C female black hispan paredu parprof parinc_prank famsize afqt3 //baseline confounders
global D att22 //exposure
global M1 ever_unemp_age3539 //first mediator
global M2 log_faminc_adj_age3539 //second mediator
global Y std_cesd_age40 //outcome

//estimate natural effects through M1 (unemployment)
preserve 

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, female, black, ///
	hispan, paredu, parprof, parinc_prank, famsize, afqt3)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//compute estimates from additive linear models
qui linmed $Y $M1, dvar($D) d(1) dstar(0) cvars($C) nointer reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//compute estimates from linear models with DxM1 interaction
qui linmed $Y $M1, dvar($D) d(1) dstar(0) cvars($C) reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//estimate natural effects through M2 (income)
restore 

//keep complete cases
drop if missing(cesd_age40, att22, log_faminc_adj_age3539, female, black, ///
	hispan, paredu, parprof, parinc_prank, famsize, afqt3)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//estimate natural effects through M2 using additive linear models
qui linmed $Y $M2, dvar($D) d(1) dstar(0) cvars($C) nointer reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//estimate natural effects through M2 using linear models with DxM2 interaction
qui linmed $Y $M2, dvar($D) d(1) dstar(0) cvars($C) reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

log close
