/*=========================================================================
DO FILE NAME:			Main_analyses1_cognition.do

AUTHOR:					Rutendo Muzambi

VERSION:				v2

DATE VERSION CREATED: 	06/2020
						
DESCRIPTION OF FILE:   
1. First ever infections on cognitive impairment
*=========================================================================*/


********************************************************************************
****************************** Data analysis**********************************
********************************************************************************

use $cleandata\main_analysis_cognition, clear
stset end, failure(cognition_status=1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

list patid _t0 _t start_fup end infection_cat infection_date1 cognition_date1 cognition_status in 1/20
list patid _t0 _t start_fup end infection_cat infection_date1 cognition_date1 cognition_status if cognition_status==1 & cognition_date1>end in 1/400

sort patid _t0 _t
by patid: gen gap=1 if _t0!=_t[_n-1] & _n>1
list patid if gap==1 in 1/20
count if gap==1 & _st==1 // 
stdescribe // gaps 0.25 (3 months) 

summarize _st 
stvary
save $cleandata\cognition_analysis_cox, replace

********************************************************************************
******************************Split on infection******************************
********************************************************************************
use $cleandata\cognition_analysis_cox, clear
replace infection_date1= end+365.25 if infection_date1==. //infection date after end date
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

cap log close
log using $log/cognition_main.log, replace


strate infection_tv, per(1000)
strate infection_tv, per(1000) output ($results/strate_cognition, replace)
save $cleandata\cognition_cox, replace

capture file close textfile_cognition 
file open textfile_cognition using $results/cognitive_impairment.csv, write replace
file write textfile_cognition "sep=;" _n
file write textfile_cognition "Effect of infections on cognitive impairment" _n _n
file write textfile_cognition  ";"   "Number of cognitive impairment events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile_cognition  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

*incidence rate*
file write textfile_cognition "Type of infection" _n
stptime if infection_tv==1, title(person-years) per(1000)
return list
**infection
local persontime_inf=`r(ptime)'
local cognitionevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

stptime if infection_tv==0, title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local cognitionevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

file write textfile_cognition "No infection" ";" (`cognitionevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile_cognition %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n



********************************************************************************
********************4.Any infection Cox regression analysis ********************
********************************************************************************

*Age adjusted analysis
stcox i.infection_tv, base //
return list
ereturn list
matrix list e(b)
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*** Calendar period
stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab


**Minimally adjusted model 
stcox i.infection_tv i.sex i.imd_cat i.calperiod, base

local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	


*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile_cognition "Any Infection" ";" (`cognitionevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile_cognition %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile_cognition ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_cognition ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_cognition ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_cognition

*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001

log close


********************************************************************************
********************4.Type of infection Cox regression analysis ********************
********************************************************************************
use $cleandata\cognition_analysis_cox, replace

***set up of type of infection
tab infection_type

replace infection_date1= end+365.25 if infection_date1==.
stsplit inf_type_tv=infection_date1, at(0)
replace inf_type_tv=inf_type_tv+1
tab inf_type_tv

list patid _t0 _t start_fup end _st _d infection_date1 inf_type_tv in 1/20
replace inf_type_tv=1 if inf_type_tv==1 & infection_type ==1
replace inf_type_tv=2 if inf_type_tv==1 & infection_type ==2
replace inf_type_tv=3 if inf_type_tv==1 & infection_type ==3
replace inf_type_tv=4 if inf_type_tv==1 & infection_type ==4
replace inf_type_tv=5 if inf_type_tv==1 & infection_type ==5
replace inf_type_tv=6 if inf_type_tv==1 & infection_type ==6

label define inf_type_lab 0 "No" 1"Sepsis" 2"Pneumonia" 3"Other LRTI" 4"UTI" 5"SSTI" 6"Unknown infection"
label values inf_type_tv inf_type_lab

tab inf_type_tv //infections on same date as enddate not counted 0, infections whereby enddate occurs within 3 months not counted,
count if (end-infection_date1)<=91.3125 //
lab var inf_type_tv "Type of Infection"
list patid _t0 _t start_fup end _st infection_type inf_type_tv infection_date1 enddate cognition_date1 in 1/20

save $cleandata\cognition_type_cox, replace

****

foreach i of num 1/6 {
stptime if inf_type_tv==`i',  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local cognitionevents_`i'=`r(failures)'
local cruderate_`i'=`r(rate)'
local lci_rate_`i'=`r(lb)'
local uci_rate_`i'=`r(ub)'
}

*Age adjusted analysis
stcox i.inf_type_tv, base //
return list
ereturn list
matrix list e(b)

local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_age_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_age_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_age_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

**Minimally adjusted model 
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod, base 

local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_min_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_min_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_min_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Fully adjusted model*
stcox i.inf_type_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base //


local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_ful_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_ful_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_ful_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

capture file close textfile_cognition
*file close textfile
file open textfile_cog_type using $results/cognitive_impairment.csv, write append
file write textfile_cog_type "Sepsis" ";" (`cognitionevents_1') ";" %3.1f (`persontime_1')  ";"
file write textfile_cog_type %3.2f (`cruderate_1') "(" %3.2f (`lci_rate_1') "-" %3.2f (`uci_rate_1') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_1') " (" %3.2f (`lci_age_1')  "-" %3.2f (`uci_age_1') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_1') " (" %3.2f (`lci_min_1')  "-" %3.2f (`uci_min_1') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_1') " (" %3.2f (`lci_ful_1')  "-" %3.2f (`uci_ful_1') ")" _n

file write textfile_cog_type "Pneumonia" ";" (`cognitionevents_2') ";" %3.1f (`persontime_2')  ";"
file write textfile_cog_type %3.2f (`cruderate_2') "(" %3.2f (`lci_rate_2') "-" %3.2f (`uci_rate_2') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_2') " (" %3.2f (`lci_age_2')  "-" %3.2f (`uci_age_2') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_2') " (" %3.2f (`lci_min_2')  "-" %3.2f (`uci_min_2') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_2') " (" %3.2f (`lci_ful_2')  "-" %3.2f (`uci_ful_2') ")" _n

file write textfile_cog_type "Other LRTI" ";" (`cognitionevents_3') ";" %3.1f (`persontime_3')  ";"
file write textfile_cog_type %3.2f (`cruderate_3') "(" %3.2f (`lci_rate_3') "-" %3.2f (`uci_rate_3') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_3') " (" %3.2f (`lci_age_3')  "-" %3.2f (`uci_age_3') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_3') " (" %3.2f (`lci_min_3')  "-" %3.2f (`uci_min_3') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_3') " (" %3.2f (`lci_ful_3')  "-" %3.2f (`uci_ful_3') ")" _n

file write textfile_cog_type "UTI" ";" (`cognitionevents_4') ";" %3.1f (`persontime_4')  ";"
file write textfile_cog_type %3.2f (`cruderate_4') "(" %3.2f (`lci_rate_4') "-" %3.2f (`uci_rate_4') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_4') " (" %3.2f (`lci_age_4')  "-" %3.2f (`uci_age_4') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_4') " (" %3.2f (`lci_min_4')  "-" %3.2f (`uci_min_4') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_4') " (" %3.2f (`lci_ful_4')  "-" %3.2f (`uci_ful_4') ")" _n

file write textfile_cog_type "SSTI" ";" (`cognitionevents_5') ";" %3.1f (`persontime_5')  ";"
file write textfile_cog_type %3.2f (`cruderate_5') "(" %3.2f (`lci_rate_5') "-" %3.2f (`uci_rate_5') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_5') " (" %3.2f (`lci_age_5')  "-" %3.2f (`uci_age_5') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_5') " (" %3.2f (`lci_min_5')  "-" %3.2f (`uci_min_5') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_5') " (" %3.2f (`lci_ful_5')  "-" %3.2f (`uci_ful_5') ")" _n

file close textfile_cog_type

*Fully adjusted model*
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001


********************************************************************************
***********Cognition diagnosis**************************************************
********************************************************************************

use $cleandata\main_analysis_cognition, clear

drop if cognition_symptoms==1
tab cognition_diagnosis cognition_status
stset end, failure(cognition_status=1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

replace infection_date1= end+365.25 if infection_date1==. //infection date after end date
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

log using $log/cognition_diagnoses.log, replace

strate infection_tv, per (1000) output ($results/cognition_diagnoses, replace)

capture file close textfile_cognition 
file open textfile_cognition using $results/cognition_diagnoses.csv, write replace
file write textfile_cognition "sep=;" _n
file write textfile_cognition "Effect of infections on diagnosis of cognitive impairment" _n _n
file write textfile_cognition  ";"   "Number of cognitive impairment events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile_cognition  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

*incidence rate*
file write textfile_cognition "Type of infection" _n
stptime if infection_tv==1, title(person-years) per(1000)
return list
**infection
local persontime_inf=`r(ptime)'
local cognitionevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

stptime if infection_tv==0, title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local cognitionevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

file write textfile_cognition "No infection" ";" (`cognitionevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile_cognition %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n

*Age adjusted analysis
stcox i.infection_tv, base //
return list
ereturn list
matrix list e(b)
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*** Calendar period
stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab


**Minimally adjusted model 
stcox i.infection_tv i.sex i.imd_cat i.calperiod, base

local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	


*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile_cognition "Any Infection" ";" (`cognitionevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile_cognition %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile_cognition ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_cognition ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_cognition ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_cognition

*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001

log close

use $cleandata/cognition_diagnoses_cox, clear
tab infection_type

replace infection_date1= end+365.25 if infection_date1==.
stsplit inf_type_tv=infection_date1, at(0)
replace inf_type_tv=inf_type_tv+1
tab inf_type_tv

list patid _t0 _t start_fup end _st _d infection_date1 inf_type_tv in 1/20
replace inf_type_tv=1 if inf_type_tv==1 & infection_type ==1
replace inf_type_tv=2 if inf_type_tv==1 & infection_type ==2
replace inf_type_tv=3 if inf_type_tv==1 & infection_type ==3
replace inf_type_tv=4 if inf_type_tv==1 & infection_type ==4
replace inf_type_tv=5 if inf_type_tv==1 & infection_type ==5
replace inf_type_tv=6 if inf_type_tv==1 & infection_type ==6

label define inf_type_lab 0 "No" 1"Sepsis" 2"Pneumonia" 3"Other LRTI" 4"UTI" 5"SSTI" 6"Unknown infection"
label values inf_type_tv inf_type_lab

tab inf_type_tv //infections on same date as enddate not counted 0, infections whereby enddate occurs within 3 months not counted,
count if (end-infection_date1)<=91.3125 //
lab var inf_type_tv "Type of Infection"
list patid _t0 _t start_fup end _st infection_type inf_type_tv infection_date1 enddate cognition_date1 in 1/20

save $cleandata\cognition_type_cox, replace

****

foreach i of num 1/6 {
stptime if inf_type_tv==`i',  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local cognitionevents_`i'=`r(failures)'
local cruderate_`i'=`r(rate)'
local lci_rate_`i'=`r(lb)'
local uci_rate_`i'=`r(ub)'
}

*Age adjusted analysis
stcox i.inf_type_tv, base //
return list
ereturn list
matrix list e(b)

local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_age_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_age_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_age_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab


**Minimally adjusted model 
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod, base 

local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_min_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_min_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_min_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Fully adjusted model*
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //


local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_ful_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_ful_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_ful_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

capture file close textfile_cognition
*file close textfile
file open textfile_cog_type using $results/cognition_diagnoses.csv, write append
file write textfile_cog_type "Sepsis" ";" (`cognitionevents_1') ";" %3.1f (`persontime_1')  ";"
file write textfile_cog_type %3.2f (`cruderate_1') "(" %3.2f (`lci_rate_1') "-" %3.2f (`uci_rate_1') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_1') " (" %3.2f (`lci_age_1')  "-" %3.2f (`uci_age_1') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_1') " (" %3.2f (`lci_min_1')  "-" %3.2f (`uci_min_1') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_1') " (" %3.2f (`lci_ful_1')  "-" %3.2f (`uci_ful_1') ")" _n

file write textfile_cog_type "Pneumonia" ";" (`cognitionevents_2') ";" %3.1f (`persontime_2')  ";"
file write textfile_cog_type %3.2f (`cruderate_2') "(" %3.2f (`lci_rate_2') "-" %3.2f (`uci_rate_2') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_2') " (" %3.2f (`lci_age_2')  "-" %3.2f (`uci_age_2') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_2') " (" %3.2f (`lci_min_2')  "-" %3.2f (`uci_min_2') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_2') " (" %3.2f (`lci_ful_2')  "-" %3.2f (`uci_ful_2') ")" _n

file write textfile_cog_type "Other LRTI" ";" (`cognitionevents_3') ";" %3.1f (`persontime_3')  ";"
file write textfile_cog_type %3.2f (`cruderate_3') "(" %3.2f (`lci_rate_3') "-" %3.2f (`uci_rate_3') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_3') " (" %3.2f (`lci_age_3')  "-" %3.2f (`uci_age_3') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_3') " (" %3.2f (`lci_min_3')  "-" %3.2f (`uci_min_3') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_3') " (" %3.2f (`lci_ful_3')  "-" %3.2f (`uci_ful_3') ")" _n

file write textfile_cog_type "UTI" ";" (`cognitionevents_4') ";" %3.1f (`persontime_4')  ";"
file write textfile_cog_type %3.2f (`cruderate_4') "(" %3.2f (`lci_rate_4') "-" %3.2f (`uci_rate_4') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_4') " (" %3.2f (`lci_age_4')  "-" %3.2f (`uci_age_4') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_4') " (" %3.2f (`lci_min_4')  "-" %3.2f (`uci_min_4') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_4') " (" %3.2f (`lci_ful_4')  "-" %3.2f (`uci_ful_4') ")" _n

file write textfile_cog_type "SSTI" ";" (`cognitionevents_5') ";" %3.1f (`persontime_5')  ";"
file write textfile_cog_type %3.2f (`cruderate_5') "(" %3.2f (`lci_rate_5') "-" %3.2f (`uci_rate_5') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_age_5') " (" %3.2f (`lci_age_5')  "-" %3.2f (`uci_age_5') ")"  
file write  textfile_cog_type ";"  %3.2f (`hr_min_5') " (" %3.2f (`lci_min_5')  "-" %3.2f (`uci_min_5') ")" 
file write  textfile_cog_type ";"  %3.2f (`hr_ful_5') " (" %3.2f (`lci_ful_5')  "-" %3.2f (`uci_ful_5') ")" _n

file close textfile_cog_type

*Fully adjusted model*
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001


