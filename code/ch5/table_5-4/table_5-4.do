/*Table 5.4*/
capture clear all
capture log close
set more off

//specify directories 
global datadir "https://github.com/causalMedAnalysis/repFiles/raw/refs/heads/main/data/NLSY79/" 
global logdir "C:\Users\Geoffrey Wodtke\Dropbox\D\projects\causal_mediation_text\code\ch5\_LOGS\"

//open log
log using "${logdir}table_5-4.log", replace 

//load data
use "${datadir}nlsy79BK_ed2.dta"

//recode maternal education
recode momedu (0/12=0) (13/20=1), gen(momcol)

//recode family income
recode faminc_adj_age3539 (0/49999.99=0) (50000/9999999=1), gen(incgt50k)

//keep complete cases
drop if missing(cesd_age40, att22, ever_unemp_age3539, momcol, incgt50k)

//standardize ces-d scores
egen std_cesd_age40=std(cesd_age40)

//tabulate data
bysort momcol att22 incgt50k: tabstat std_cesd_age40, by(ever_unemp_age3539) s(mean count)

log close
