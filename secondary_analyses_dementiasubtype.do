/*=========================================================================
DO FILE NAME:			Sensitivity_analysis_dementiasubtype.do

AUTHOR:					Rutendo Muzambi

VERSION:				v1

DATE VERSION CREATED: 	06/2020
						

*=========================================================================*/

**********************************************************************************************************************
********************************************1. Alzheimers disease*****************************************************
**********************************************************************************************************************

*description
use $cleandata\main_analysis1, clear
log using $log\alzheimers.log, replace

*analysis
stset end, failure(alzheimers==1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

replace infection_date1= end+365.25 if infection_date1==. //
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

list patid start_fup end infection_date1 infection_tv dementia_date1 alzheimers dementia_status if alzheimers==1 in 1/2600

*descriptive
strate infection_tv, per(1000) 
strate infection_tv, per(1000) output ($results/strate_alzheimers, replace)

capture file close textfile_alz 
file open textfile_alz using $results/dementia_subtype_cox.csv, write replace
file write textfile_alz "sep=;" _n
file write textfile_alz "Effect of common infections on dementia subtype" _n _n
file write textfile_alz  ";"   "Number of events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile_alz  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

*incidence rate*
file write textfile_alz "Type of infection" _n
file write textfile_alz "Alzheimer's Disease" _n
stptime if infection_tv==1, title(person-years) per(1000)
return list
**infection
local persontime_inf=`r(ptime)'
local dementiaevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

stptime if infection_tv==0, title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local dementiaevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

file write textfile_alz "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile_alz %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n


*Age adjusted analysis
stcox i.infection_tv, base //
return list
ereturn list
matrix list e(b)
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

**Minimally adjusted model 
stcox i.infection_tv i.sex i.imd_cat i.calperiod, base 

local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile_alz "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile_alz %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile_alz ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_alz ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_alz ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_alz
log close

use $cleandata\main_analysis1, clear

*analysis
stset end, failure(alzheimers==1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup)

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

tab inf_type_tv //
lab var inf_type_tv "Type of Infection"

foreach i of num 1/6 {
stptime if inf_type_tv==`i',  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local dementiaevents_`i'=`r(failures)'
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

capture file close textfile 
*file close textfile
file open textfile_type_alz using $results/dementia_subtype_cox.csv, write append
file write textfile_type_alz "Sepsis" ";" (`dementiaevents_1') ";" %3.1f (`persontime_1')  ";"
file write textfile_type_alz %3.2f (`cruderate_1') "(" %3.2f (`lci_rate_1') "-" %3.2f (`uci_rate_1') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_age_1') " (" %3.2f (`lci_age_1')  "-" %3.2f (`uci_age_1') ")"  
file write  textfile_type_alz ";"  %3.2f (`hr_min_1') " (" %3.2f (`lci_min_1')  "-" %3.2f (`uci_min_1') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_ful_1') " (" %3.2f (`lci_ful_1')  "-" %3.2f (`uci_ful_1') ")" _n

file write textfile_type_alz "Pneumonia" ";" (`dementiaevents_2') ";" %3.1f (`persontime_2')  ";"
file write textfile_type_alz %3.2f (`cruderate_2') "(" %3.2f (`lci_rate_2') "-" %3.2f (`uci_rate_2') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_age_2') " (" %3.2f (`lci_age_2')  "-" %3.2f (`uci_age_2') ")"  
file write  textfile_type_alz ";"  %3.2f (`hr_min_2') " (" %3.2f (`lci_min_2')  "-" %3.2f (`uci_min_2') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_ful_2') " (" %3.2f (`lci_ful_2')  "-" %3.2f (`uci_ful_2') ")" _n

file write textfile_type_alz "Other LRTI" ";" (`dementiaevents_3') ";" %3.1f (`persontime_3')  ";"
file write textfile_type_alz %3.2f (`cruderate_3') "(" %3.2f (`lci_rate_3') "-" %3.2f (`uci_rate_3') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_age_3') " (" %3.2f (`lci_age_3')  "-" %3.2f (`uci_age_3') ")"  
file write  textfile_type_alz ";"  %3.2f (`hr_min_3') " (" %3.2f (`lci_min_3')  "-" %3.2f (`uci_min_3') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_ful_3') " (" %3.2f (`lci_ful_3')  "-" %3.2f (`uci_ful_3') ")" _n

file write textfile_type_alz "UTI" ";" (`dementiaevents_4') ";" %3.1f (`persontime_4')  ";"
file write textfile_type_alz %3.2f (`cruderate_4') "(" %3.2f (`lci_rate_4') "-" %3.2f (`uci_rate_4') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_age_4') " (" %3.2f (`lci_age_4')  "-" %3.2f (`uci_age_4') ")"  
file write  textfile_type_alz ";"  %3.2f (`hr_min_4') " (" %3.2f (`lci_min_4')  "-" %3.2f (`uci_min_4') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_ful_4') " (" %3.2f (`lci_ful_4')  "-" %3.2f (`uci_ful_4') ")" _n

file write textfile_type_alz "SSTI" ";" (`dementiaevents_5') ";" %3.1f (`persontime_5')  ";"
file write textfile_type_alz %3.2f (`cruderate_5') "(" %3.2f (`lci_rate_5') "-" %3.2f (`uci_rate_5') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_age_5') " (" %3.2f (`lci_age_5')  "-" %3.2f (`uci_age_5') ")"  
file write  textfile_type_alz ";"  %3.2f (`hr_min_5') " (" %3.2f (`lci_min_5')  "-" %3.2f (`uci_min_5') ")" 
file write  textfile_type_alz ";"  %3.2f (`hr_ful_5') " (" %3.2f (`lci_ful_5')  "-" %3.2f (`uci_ful_5') ")" _n

file close textfile_type_alz

*Fully adjusted model*
stcox i.inf_type_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001


**********************************************************************************************************************
********************************************2. Vascular disease*******************************************************
**********************************************************************************************************************

use $cleandata\main_analysis1, clear
cap log close
log using $log\vascular_dem.log, replace

*analysis
stset end, failure(vascular_dem==1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

replace infection_date1= end+365.25 if infection_date1==. //infection date after end date
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

list patid _t0 _t _d end infection_date1 infection_tv dementia_date1 vascular_dem if _d==1 in 1/3000


*descriptive
strate infection_tv, per(1000) 
strate infection_tv, per(1000) output ($results/strate_vascular_dem, replace)

capture file close textfile_vasc 
file open textfile_vasc using $results/dementia_subtype_cox.csv, write append

*incidence rate*
file write textfile_vasc "Vascular Dementia" _n
stptime if infection_tv==1, title(person-years) per(1000)
return list
**infection
local persontime_inf=`r(ptime)'
local dementiaevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

stptime if infection_tv==0, title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local dementiaevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

file write textfile_vasc "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile_vasc %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n


*Age adjusted analysis
stcox i.infection_tv, base //
return list
ereturn list
matrix list e(b)
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*calender period
stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

**Minimally adjusted model 
stcox i.infection_tv i.sex i.imd_cat i.calperiod, base 

local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile_vasc "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile_vasc %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile_vasc ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_vasc ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_vasc ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_vasc
log close

*******Type - Vascular dementia

use $cleandata\main_analysis1, clear
*analysis
stset end, failure(vascular_dem==1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

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

tab inf_type_tv //
count if (end-infection_date1)<=91.3125 //
lab var inf_type_tv "Infection type"

foreach i of num 1/6 {
stptime if inf_type_tv==`i',  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local dementiaevents_`i'=`r(failures)'
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


capture file close textfile 
*file close textfile
file open textfile_type_vasc using $results/dementia_subtype_cox.csv, write append
file write textfile_type_vasc "Sepsis" ";" (`dementiaevents_1') ";" %3.1f (`persontime_1')  ";"
file write textfile_type_vasc %3.2f (`cruderate_1') "(" %3.2f (`lci_rate_1') "-" %3.2f (`uci_rate_1') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_age_1') " (" %3.2f (`lci_age_1')  "-" %3.2f (`uci_age_1') ")"  
file write  textfile_type_vasc ";"  %3.2f (`hr_min_1') " (" %3.2f (`lci_min_1')  "-" %3.2f (`uci_min_1') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_ful_1') " (" %3.2f (`lci_ful_1')  "-" %3.2f (`uci_ful_1') ")" _n

file write textfile_type_vasc "Pneumonia" ";" (`dementiaevents_2') ";" %3.1f (`persontime_2')  ";"
file write textfile_type_vasc %3.2f (`cruderate_2') "(" %3.2f (`lci_rate_2') "-" %3.2f (`uci_rate_2') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_age_2') " (" %3.2f (`lci_age_2')  "-" %3.2f (`uci_age_2') ")"  
file write  textfile_type_vasc ";"  %3.2f (`hr_min_2') " (" %3.2f (`lci_min_2')  "-" %3.2f (`uci_min_2') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_ful_2') " (" %3.2f (`lci_ful_2')  "-" %3.2f (`uci_ful_2') ")" _n

file write textfile_type_vasc "Other LRTI" ";" (`dementiaevents_3') ";" %3.1f (`persontime_3')  ";"
file write textfile_type_vasc %3.2f (`cruderate_3') "(" %3.2f (`lci_rate_3') "-" %3.2f (`uci_rate_3') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_age_3') " (" %3.2f (`lci_age_3')  "-" %3.2f (`uci_age_3') ")"  
file write  textfile_type_vasc ";"  %3.2f (`hr_min_3') " (" %3.2f (`lci_min_3')  "-" %3.2f (`uci_min_3') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_ful_3') " (" %3.2f (`lci_ful_3')  "-" %3.2f (`uci_ful_3') ")" _n

file write textfile_type_vasc "UTI" ";" (`dementiaevents_4') ";" %3.1f (`persontime_4')  ";"
file write textfile_type_vasc %3.2f (`cruderate_4') "(" %3.2f (`lci_rate_4') "-" %3.2f (`uci_rate_4') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_age_4') " (" %3.2f (`lci_age_4')  "-" %3.2f (`uci_age_4') ")"  
file write  textfile_type_vasc ";"  %3.2f (`hr_min_4') " (" %3.2f (`lci_min_4')  "-" %3.2f (`uci_min_4') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_ful_4') " (" %3.2f (`lci_ful_4')  "-" %3.2f (`uci_ful_4') ")" _n

file write textfile_type_vasc "SSTI" ";" (`dementiaevents_5') ";" %3.1f (`persontime_5')  ";"
file write textfile_type_vasc %3.2f (`cruderate_5') "(" %3.2f (`lci_rate_5') "-" %3.2f (`uci_rate_5') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_age_5') " (" %3.2f (`lci_age_5')  "-" %3.2f (`uci_age_5') ")"  
file write  textfile_type_vasc ";"  %3.2f (`hr_min_5') " (" %3.2f (`lci_min_5')  "-" %3.2f (`uci_min_5') ")" 
file write  textfile_type_vasc ";"  %3.2f (`hr_ful_5') " (" %3.2f (`lci_ful_5')  "-" %3.2f (`uci_ful_5') ")" _n

file close textfile_type_vasc

*Fully adjusted model*
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001



*******************************************************************************************************************************
*************************************************3. Unspecified Dementia*******************************************************
*******************************************************************************************************************************

*description
use $cleandata\main_analysis1, clear
cap log close
log using $log\unspecified_dem.log, replace

*analysis
stset end, failure(unspecified_dem==1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

replace infection_date1= end+365.25 if infection_date1==. //infection date after end date
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

list patid _t0 _t _d end infection_date1 infection_tv dementia_date1 unspecified_dem if _d==1 in 1/3000

strate infection_tv, per(1000) 
strate infection_tv, per(1000) output ($results/strate_unspecified_dem, replace)

capture file close textfile_unsp 
file open textfile_unsp using $results/dementia_subtype_cox.csv, write append

*incidence rate*
file write textfile_unsp "Unspecified Dementia" _n
stptime if infection_tv==1, title(person-years) per(1000)
return list
**infection
local persontime_inf=`r(ptime)'
local dementiaevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

stptime if infection_tv==0, title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local dementiaevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

file write textfile_unsp "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile_unsp %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n

*Age adjusted analysis
stcox i.infection_tv, base //
return list
ereturn list
matrix list e(b)
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab


**Minimally adjusted model 
stcox i.infection_tv i.sex i.imd_cat i.calperiod, base 

local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile_unsp "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile_unsp %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile_unsp ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_unsp ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_unsp ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_unsp

*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001

log close


use $cleandata\main_analysis1, clear

*analysis
stset end, failure(unspecified_dem==1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

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

tab inf_type_tv //
count if (end-infection_date1)<=91.3125 //
lab var inf_type_tv "Infection type"

foreach i of num 1/6 {
stptime if inf_type_tv==`i',  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local dementiaevents_`i'=`r(failures)'
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


capture file close textfile 
*file close textfile
file open textfile_type_unsp using $results/dementia_subtype_cox.csv, write append
file write textfile_type_unsp "Sepsis" ";" (`dementiaevents_1') ";" %3.1f (`persontime_1')  ";"
file write textfile_type_unsp %3.2f (`cruderate_1') "(" %3.2f (`lci_rate_1') "-" %3.2f (`uci_rate_1') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_age_1') " (" %3.2f (`lci_age_1')  "-" %3.2f (`uci_age_1') ")"  
file write  textfile_type_unsp ";"  %3.2f (`hr_min_1') " (" %3.2f (`lci_min_1')  "-" %3.2f (`uci_min_1') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_ful_1') " (" %3.2f (`lci_ful_1')  "-" %3.2f (`uci_ful_1') ")" _n

file write textfile_type_unsp "Pneumonia" ";" (`dementiaevents_2') ";" %3.1f (`persontime_2')  ";"
file write textfile_type_unsp %3.2f (`cruderate_2') "(" %3.2f (`lci_rate_2') "-" %3.2f (`uci_rate_2') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_age_2') " (" %3.2f (`lci_age_2')  "-" %3.2f (`uci_age_2') ")"  
file write  textfile_type_unsp ";"  %3.2f (`hr_min_2') " (" %3.2f (`lci_min_2')  "-" %3.2f (`uci_min_2') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_ful_2') " (" %3.2f (`lci_ful_2')  "-" %3.2f (`uci_ful_2') ")" _n

file write textfile_type_unsp "Other LRTI" ";" (`dementiaevents_3') ";" %3.1f (`persontime_3')  ";"
file write textfile_type_unsp %3.2f (`cruderate_3') "(" %3.2f (`lci_rate_3') "-" %3.2f (`uci_rate_3') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_age_3') " (" %3.2f (`lci_age_3')  "-" %3.2f (`uci_age_3') ")"  
file write  textfile_type_unsp ";"  %3.2f (`hr_min_3') " (" %3.2f (`lci_min_3')  "-" %3.2f (`uci_min_3') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_ful_3') " (" %3.2f (`lci_ful_3')  "-" %3.2f (`uci_ful_3') ")" _n

file write textfile_type_unsp "UTI" ";" (`dementiaevents_4') ";" %3.1f (`persontime_4')  ";"
file write textfile_type_unsp %3.2f (`cruderate_4') "(" %3.2f (`lci_rate_4') "-" %3.2f (`uci_rate_4') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_age_4') " (" %3.2f (`lci_age_4')  "-" %3.2f (`uci_age_4') ")"  
file write  textfile_type_unsp ";"  %3.2f (`hr_min_4') " (" %3.2f (`lci_min_4')  "-" %3.2f (`uci_min_4') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_ful_4') " (" %3.2f (`lci_ful_4')  "-" %3.2f (`uci_ful_4') ")" _n

file write textfile_type_unsp "SSTI" ";" (`dementiaevents_5') ";" %3.1f (`persontime_5')  ";"
file write textfile_type_unsp %3.2f (`cruderate_5') "(" %3.2f (`lci_rate_5') "-" %3.2f (`uci_rate_5') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_age_5') " (" %3.2f (`lci_age_5')  "-" %3.2f (`uci_age_5') ")"  
file write  textfile_type_unsp ";"  %3.2f (`hr_min_5') " (" %3.2f (`lci_min_5')  "-" %3.2f (`uci_min_5') ")" 
file write  textfile_type_unsp ";"  %3.2f (`hr_ful_5') " (" %3.2f (`lci_ful_5')  "-" %3.2f (`uci_ful_5') ")" _n


file close textfile_type_unsp

*Fully adjusted model*
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001


*******************************************************************************************************************************
*************************************************4. Exclude cause-specific dem*******************************************************
*******************************************************************************************************************************


use $cleandata\main_analysis1, clear
drop if cause_identified_dem==1

stset end, failure(dementia_status=1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 


*2.	Split on infection

replace infection_date1= end+365.25 if infection_date1==. //infection date after end date
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

*3. Analysis 

cap log close
log using $log\excl_causedem_cox.log, replace

*st
strate infection_tv, per(1000) 
strate infection_tv, per(1000) output ($results/strate_exclude_cause, replace)


label define infection_tb 0 "No infection" 1"Infection" 
label values infection_tv infection_tv

*file close textfile

*tempname myfile
capture file close textfile 
file open textfile using $results/exclude_causedem_cox.csv, write replace
file write textfile "sep=;" _n
file write textfile "Crude incidence rate and hazard ratios for the association between common infections and dementia, exclusion of cause-specific dementia" _n _n
file write textfile  ";"   "Number of events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

*incidence rate*
file write textfile "Type of infection" _n
stptime if infection_tv==1, title(person-years) per(1000)
return list
**infection
local persontime_inf=`r(ptime)'
local dementiaevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

stptime if infection_tv==0, title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local dementiaevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

file write textfile "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n

*Age adjusted analysis
stcox i.infection_tv, base //
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

**Minimally adjusted model 
stcox i.infection_tv i.sex i.imd_cat i.calperiod, base 
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	


file write textfile "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile

log close


use $cleandata\main_analysis1, clear
drop if cause_identified_dem==1

stset end, failure(dementia_status=1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

**Type***

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

tab inf_type_tv //
count if (end-infection_date1)<=91.3125 //
lab var inf_type_tv "Infection type"
list patid _t0 _t start_fup end _st infection_type inf_type_tv infection_date1 enddate dementia_date1 in 1/20


foreach i of num 1/6 {
stptime if inf_type_tv==`i',  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local dementiaevents_`i'=`r(failures)'
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

stsplit calperiod, after(mdy(1,1,2004)) at(0,4,8)
recode calperiod 0=2004 4=2008 8=2012


**Minimally adjusted model 
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod, base 

local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_min_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_min_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_min_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Fully adjusted model*
stcox i.inf_type_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //


local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_ful_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_ful_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_ful_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}


capture file close textfile 
file open textfile_type using $results/exclude_causedem_cox.csv, write append
file write textfile_type "Sepsis" ";" (`dementiaevents_1') ";" %3.1f (`persontime_1')  ";"
file write textfile_type %3.2f (`cruderate_1') "(" %3.2f (`lci_rate_1') "-" %3.2f (`uci_rate_1') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_1') " (" %3.2f (`lci_age_1')  "-" %3.2f (`uci_age_1') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_1') " (" %3.2f (`lci_min_1')  "-" %3.2f (`uci_min_1') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_1') " (" %3.2f (`lci_ful_1')  "-" %3.2f (`uci_ful_1') ")" _n


file write textfile_type "Pneumonia" ";" (`dementiaevents_2') ";" %3.1f (`persontime_2')  ";"
file write textfile_type %3.2f (`cruderate_2') "(" %3.2f (`lci_rate_2') "-" %3.2f (`uci_rate_2') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_2') " (" %3.2f (`lci_age_2')  "-" %3.2f (`uci_age_2') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_2') " (" %3.2f (`lci_min_2')  "-" %3.2f (`uci_min_2') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_2') " (" %3.2f (`lci_ful_2')  "-" %3.2f (`uci_ful_2') ")" _n  

file write textfile_type "Other LRTI" ";" (`dementiaevents_3') ";" %3.1f (`persontime_3')  ";"
file write textfile_type %3.2f (`cruderate_3') "(" %3.2f (`lci_rate_3') "-" %3.2f (`uci_rate_3') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_3') " (" %3.2f (`lci_age_3')  "-" %3.2f (`uci_age_3') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_3') " (" %3.2f (`lci_min_3')  "-" %3.2f (`uci_min_3') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_3') " (" %3.2f (`lci_ful_3')  "-" %3.2f (`uci_ful_3') ")" _n

file write textfile_type "UTI" ";" (`dementiaevents_4') ";" %3.1f (`persontime_4')  ";"
file write textfile_type %3.2f (`cruderate_4') "(" %3.2f (`lci_rate_4') "-" %3.2f (`uci_rate_4') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_4') " (" %3.2f (`lci_age_4')  "-" %3.2f (`uci_age_4') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_4') " (" %3.2f (`lci_min_4')  "-" %3.2f (`uci_min_4') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_4') " (" %3.2f (`lci_ful_4')  "-" %3.2f (`uci_ful_4') ")" _n

file write textfile_type "SSTI" ";" (`dementiaevents_5') ";" %3.1f (`persontime_5')  ";"
file write textfile_type %3.2f (`cruderate_5') "(" %3.2f (`lci_rate_5') "-" %3.2f (`uci_rate_5') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_5') " (" %3.2f (`lci_age_5')  "-" %3.2f (`uci_age_5') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_5') " (" %3.2f (`lci_min_5')  "-" %3.2f (`uci_min_5') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_5') " (" %3.2f (`lci_ful_5')  "-" %3.2f (`uci_ful_5') ")" _n

file close textfile_type


