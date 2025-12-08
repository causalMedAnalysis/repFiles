/*Table 4.5*/
capture clear all
capture log close
set more off
set maxvar 50000

//install required modules
net install cmed, from("https://raw.github.com/causalMedAnalysis/cmed/master/") replace //module for causal mediation analysis
net install parallel, from("https://raw.github.com/gvegayon/parallel/stable/") replace //module for parallelization
mata mata mlib index

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/NLSY79/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch4\_LOGS\"

//open log
log using "${logdir}table_4-5.log", replace 

//load data
use "${datadir}nlsy79BK_ed2.dta"

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

//set seed 
//note: opt parallel requires a separate seed for each child process
local myseed 3308004
qui parallel numprocessors
local ncores = max(floor(r(numprocessors) * 0.75), 1)
local parseeds 
forval i = 1/`ncores' {
    local newseed = `myseed' + `i'
    local parseeds `parseeds' " `newseed'"
}

//compute RWR interval estimates
qui cmed linear $Y $M ($L) $D = $C, cxd cxm lxm ///
	reps(2000) parallel seed(`parseeds') saving("bootrwr.dta", replace)

mat list e(ci_percentile)

qui cmed linear $Y $M ($L) $D = $C, m(10.82) cxd cxm lxm ///
	reps(2000) parallel seed(`parseeds') saving("bootrwr_cde.dta", replace)

mat list e(ci_percentile)

//compute simulation interval estimates
qui cmed sim ((regress) $Y) ((regress) $M) ((logit) $L) $D = $C, cxd cxm lxm nsim(2000) ///
	reps(2000) parallel seed(`parseeds') saving("bootsim.dta", replace)

mat list e(ci_percentile)

qui cmed sim ((regress) $Y) $M ((logit) $L) $D = $C, mvalue(10.82) cxd cxm lxm nsim(2000) ///
	reps(2000) parallel seed(`parseeds') saving("bootsim_cde.dta", replace)

mat list e(ci_percentile)

//compute IPW interval estimates
qui cmed ipw $Y ((regress) $M) ((logit) $L) $D = $C, censor(1 99) cxd lxd ///
	reps(2000) parallel seed(`parseeds') saving("bootipw.dta", replace)

mat list e(ci_percentile)

qui cmed ipw $Y ((regress) $M) ($L) $D = $C, m(10.82) censor(1 99) cxd lxd ///
	reps(2000) parallel seed(`parseeds') saving("bootipw_cde.dta", replace)

mat list e(ci_percentile)

//compute p-values for null of no effect
qui use "bootrwr.dta", clear
qui merge 1:1 _n using "bootrwr_cde.dta"
save "bootrwr.dta", replace
erase "bootrwr_cde.dta"

qui use "bootsim.dta", clear
qui merge 1:1 _n using "bootsim_cde.dta"
save "bootsim.dta", replace
erase "bootsim_cde.dta"

qui use "bootipw.dta", clear
qui merge 1:1 _n using "bootipw_cde.dta"
save "bootipw.dta", replace
erase "bootipw_cde.dta"

foreach method in rwr sim ipw {
	di "P-values for `method' estimates:"
	qui use "boot`method'.dta", clear
	
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
	erase "boot`method'.dta"
}

log close

//some of the estimates differ slightly from those reported in the text, 
//which are based on the R implementation. This is variously due to minor 
//differences in how the weights are censored, to Monte Carlo error, or to 
//differences in random number seeding that influence the bootstrap samples

//note: -cmed sim- with large nsim() and reps() can be time consuming to run;
//parallelization of the bootstrap is highly recommended in these instances
