/*Figure 3.7*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/cmed //module to perform causal mediation analysis

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch3\_LOGS\"
global figdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\figures\ch3\" 

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/NLSY79/nlsy79BK_ed2.dta" ///
	"${datadir}NLSY79\"

//open log
log using "${logdir}figure_3-7.log", replace 

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

//compute bootstrap estimates based on IPW
qui cmed ipw $Y $M $D = $C, censor(1 99) reps(2000) seed(3308004) ///
	saving("${datadir}\bootmed.dta", replace)

//plot histogram of bootstrap estimates
use "${datadir}\bootmed.dta", clear

centile NIE, c(2.5 97.5)

set scheme s2mono

hist NIE, ///
	width(0.00125) ///
	yscale(range(0 180)) ///
	xscale(range(-0.05 0.02)) ///
	ylabel(0(20)180) ///
	xlabel(-0.05(0.01)0.02) ///
	xline(-.0295328 .0033694, lcolor(black) lwidth(1pt) lpattern(dash)) ///
	ytitle("Frequency") freq /// 
	xtitle("`=ustrunescape("\u03B8\u0302")'{subscript:b}")

graph export "${figdir}\figure_3-7.pdf", replace

//erase bootstrap data	
erase "${datadir}\bootmed.dta"

log close

