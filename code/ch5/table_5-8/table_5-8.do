/*Table 5.8*/
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
copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/Brader_et_al2008/Brader_et_al2008.dta" ///
	"${datadir}Brader_et_al2008\"

//open log
log using "${logdir}table_5-8.log", replace 

//load data
use "${datadir}Brader_et_al2008\Brader_et_al2008.dta", clear

//keep complete cases
drop if missing(immigr, emo, p_harm, tone_eth, ppage, ppeducat, ppgender, ppincimp)

//standardize outcome
egen std_immigr=std(4-immigr)

//dummy code controls
tab ppeducat, gen(ppeducat_)
drop ppeducat_1
rename ppeducat_2 hs
rename ppeducat_3 sc
rename ppeducat_4 ba

tab ppgender, gen(ppgender_)
drop ppgender_1
rename ppgender_2 female

//define macros for different variables
global C ppage female hs sc ba ppincimp //baseline confounders
global D tone_eth //exposure
global M1 p_harm //first mediator
global M2 emo //second mediator
global Y std_immigr //outcome

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




