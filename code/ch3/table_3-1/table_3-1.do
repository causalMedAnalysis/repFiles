/*Table 3.1*/
capture clear all
capture log close
set more off

//office
*global datadir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\data\" 
*global logdir "C:\Users\Geoffrey Wodtke\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

//home
global datadir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\data\" 
global logdir "C:\Users\Geoff\Dropbox\shared\causal_mediation_text\code\ch3\_LOGS\"

log using "${logdir}table_3-1.log", replace 

//input data
use "${datadir}NLSY79\nlsy79BK_ed2.dta"

recode momedu (0/12=0) (13/20=1), gen(momcol)

drop if ever_unemp_age3539==. | momcol==. | att22==. | cesd_age40==.

egen std_cesd_age40=std(cesd_age40)

//tabulate data
bysort momcol att22: tabstat std_cesd_age40, by(ever_unemp_age3539) s(mean count)

//nonparametric estimators

//ATEhat^np
quietly reg std_cesd_age40 i.att22##i.momcol
margins, dydx(att22) 

//CDEhat^np
quietly reg std_cesd_age40 i.att22##i.ever_unemp_age3539##i.momcol
margins, dydx(att22) at(ever_unemp_age3539=(0 1)) 

//NDEhat^np and NIEhat^np
gen att22_orig=att22
gen ever_unemp_age3539_orig=ever_unemp_age3539
 
replace att22=0
replace ever_unemp_age3539=0
predict EhatY_D0M0C

replace ever_unemp_age3539=1
predict EhatY_D0M1C

replace att22=1
replace ever_unemp_age3539=0
predict EhatY_D1M0C

replace ever_unemp_age3539=1
predict EhatY_D1M1C

replace att22=att22_orig 
replace ever_unemp_age3539=ever_unemp_age3539_orig

quietly reg ever_unemp_age3539 i.att22##i.momcol

replace att22=0
predict PhatM_D0C

replace att22=1
predict PhatM_D1C

gen NDEhat_C=(EhatY_D1M0C-EhatY_D0M0C)*(1-PhatM_D0C)+(EhatY_D1M1C-EhatY_D0M1C)*PhatM_D0C
gen NIEhat_C=(EhatY_D1M0C)*[(1-PhatM_D1C)-(1-PhatM_D0C)]+(EhatY_D1M1C)*[(PhatM_D1C)-(PhatM_D0C)]

tabstat NDEhat_C NIEhat_C, s(mean) c(s)

log close
