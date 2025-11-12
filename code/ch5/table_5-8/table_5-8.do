/*Table 5.8*/
capture clear all
capture log close
set more off

//install required modules
net install github, from("https://haghish.github.io/github/")
github install causalMedAnalysis/cmed //module to perform causal mediation analysis

//specify directories 
global datadir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\data\" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch5\_LOGS\"

//download data
capture copy "https://github.com/causalMedAnalysis/repFiles/raw/main/data/Brader_et_al2008/Brader_et_al2008.dta" ///
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
qui tab ppeducat, gen(ppeducat_)
drop ppeducat_1
rename ppeducat_2 hs
rename ppeducat_3 sc
rename ppeducat_4 ba

qui tab ppgender, gen(ppgender_)
drop ppgender_1
rename ppgender_2 female

//define macros for different variables
global C ppage female hs sc ba ppincimp //baseline confounders
global D tone_eth //exposure
global M1 p_harm //first mediator
global M2 emo //second mediator
global Y std_immigr //outcome

//compute PSE estimates from additive linear models
qui cmed linear $Y ($M1 $M2) $D = $C, paths nointer reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//compute PSE estimates from linear model estimates w/ Dx(M1,M2) interactions
qui cmed linear $Y ($M1 $M2) $D = $C, paths reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

//compute PSE estimates using inverse probability weighting
qui cmed ipw $Y ($M1 $M2) $D = $C, paths censor(1 99) reps(2000) seed(60637)

mat list e(b)
mat list e(ci_percentile)

log close




