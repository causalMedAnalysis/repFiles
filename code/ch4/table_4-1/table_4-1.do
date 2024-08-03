/*Table 4.1*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch4\_LOGS\"

log using "${logdir}table_4-1.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta"

recode momedu (0/12=0) (13/20=1), gen(momcol)

recode faminc_adj_age3539 (0/49999.99=0) (50000/9999999=1), gen(incgt50k)

drop if ever_unemp_age3539==. | momcol==. | att22==. | cesd_age40==. | incgt50k==.

egen std_cesd_age40=std(cesd_age40)

//tabulate data
bysort momcol att22 incgt50k: tabstat std_cesd_age40, by(ever_unemp_age3539) s(mean count)

log close
