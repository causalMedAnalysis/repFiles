/*Table 4.5*/
capture clear all
capture log close
set more off
set maxvar 50000

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/rwrlite, replace 
github install causalMedAnalysis/ventsim, replace 
github install causalMedAnalysis/ipwvent, replace

//specify directories 
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

//download data
copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}table_4-5.log", replace 

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

//compute RWR interval estimates
qui rwrlite $Y $L, dvar($D) mvar($M) cvars($C) d(1) dstar(0) m(10.82) cxd cxm lxm ///
	reps(2000) seed(8675309) saving("${datadir}\bootrwr.dta", replace)

mat list e(ci_percentile)

//compute simulation interval estimates
set seed 3308004

qui ventsim $Y, dvar($D) mvar($M) lvars($L) cvars($C) d(1) dstar(0) m(10.82) ///
	mreg(regress) yreg(regress) lregs(logit) nsim(2000) cxd cxm lxm ///
	reps(2000) seed(8675309) saving("${datadir}\bootsim.dta", replace)

mat list e(ci_percentile)

//compute IPW interval estimates
qui ipwvent $Y, dvar($D) mvar($M) lvar($L) cvars($C) d(1) dstar(0) m(10.82) ///
	mreg(regress) lreg(logit) censor(1 99) cxd lxd ///
	reps(2000) seed(8675309) saving("${datadir}\bootipw.dta", replace)

mat list e(ci_percentile)

//compute p-values for null of no effect
foreach method in rwr sim ipw {
	di "P-values for `method' estimates:"
	qui use "${datadir}\boot`method'.dta", clear
	
	foreach x in OE IDE IIE CDE {
		quietly gen `x'_lt=0 
		quietly replace `x'_lt=1 if `x'<0
		quietly sum `x'_lt
		scalar p_lt = r(mean)
	
		quietly gen `x'_rt=0 
		quietly replace `x'_rt=1 if `x'>0
		quietly sum `x'_rt
		scalar p_rt = r(mean)
	
		di "`x' p-value"
		di min(p_rt, p_lt)*2
		di "-----------"
	
		drop `x'_*
	}
	
	capture clear
	erase "${datadir}\boot`method'.dta"
}

log close

//some of the estimates differ slightly from those reported in the text, 
//which are based on the R implementation. This is variously due to slight 
//differences in how the weights are censored, to Monte Carlo error, or to 
//differences in random number seeding that influence the bootstrap samples

//note: ventsim with large nsim() and reps() can be time consuming to run in 
//Stata because it cannot parallelize the bootstrap replications. Consider 
//switching to R if computation time is a concern.
