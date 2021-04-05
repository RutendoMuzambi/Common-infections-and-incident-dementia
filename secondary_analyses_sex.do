/*=========================================================================
DO FILE NAME:			sensitivity_analyses_gender.do

AUTHOR:					Rutendo Muzambi

VERSION:				v2

DATE VERSION CREATED: 	07/2020
						
DESCRIPTION OF FILE:    Infections on dementia by sex

sex==0 //male
sex==1 //female
*=========================================================================*/


**********************************************************************************************************************
*************************************************FEMALE sex==1********************************************************
**********************************************************************************************************************

*=================================ANY INFECTION=================================

use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log\female_cox.log, replace

capture file close textfile_female
file open textfile_female using $results/sex_cox.csv, write replace
file write textfile_female "sep=;" _n
file write textfile_female "Effect of common infections on dementia, by sex" _n _n
file write textfile_female  ";"   "Number of events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile_female  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

strate sex infection_tv, per (1000) output ($results/female, replace)

file write textfile_female "Female" _n
file write textfile_female "Type of infection" _n

****diabetes
stptime if infection_tv==0 & sex==1,  title(person-years) per(1000)
return list
*no infection
local persontime_no_female=`r(ptime)'
local dementiaevents_no_female=`r(failures)'
local cruderate_no_female=`r(rate)'
local lci_rate_no_female=`r(lb)'
local uci_rate_no_female=`r(ub)'

stptime if infection_tv==1 & sex==1,  title(person-years) per(1000)
**infection
local persontime_inf_female=`r(ptime)'
local dementiaevents_inf_female=`r(failures)'
local cruderate_inf_female=`r(rate)'
local lci_rate_inf_female=`r(lb)'
local uci_rate_inf_female=`r(ub)'

file write textfile_female "No infection" ";" (`dementiaevents_no_female') ";" %3.1f (`persontime_no_female') ";"
file write textfile_female %3.2f (`cruderate_no_female') "(" %3.2f (`lci_rate_no_female') "-" %3.2f (`uci_rate_no_female') ")" _n

*age adjusted
stcox i.infection_tv  if sex==1, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.imd_cat i.calperiod  if sex==1, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==1, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_female "Any Infection" ";" (`dementiaevents_inf_female') ";" %3.1f (`persontime_inf_female')  ";"
file write textfile_female %3.2f (`cruderate_inf_female') "(" %3.2f (`lci_rate_inf_female') "-" %3.2f (`uci_rate_inf_female') ")"
file write  textfile_female ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_female ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_female ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_female

*full model
stcox i.infection_tv   i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==1, base //
estimates store a
stcox  i.ethnicity_cat i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==1, base
estimates store b
lrtest a b // P < 0.001

log close

*==============================INFECTION TYPE===================================
use $cleandata\mainanalysis_type, clear

cap log close
log using $log\female_type.log, replace
strate sex inf_type_tv, per (1000) output ($results/type_female, replace)

foreach i of num 1/6 {
stptime if inf_type_tv==`i' & sex==1,  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local dementiaevents_`i'=`r(failures)'
local cruderate_`i'=`r(rate)'
local lci_rate_`i'=`r(lb)'
local uci_rate_`i'=`r(ub)'
}

*age adjusted
stcox i.inf_type_tv  if sex==1, base 
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

*minimally adjusted
stcox i.inf_type_tv i.imd_cat i.calperiod  if sex==1, base
local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_min_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_min_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_min_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Fully adjusted 
stcox i.inf_type_tv i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==1, base 
foreach x of local infectionregression_var {
local hr_ful_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_ful_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_ful_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

capture file close textfile_type_female
*file close textfile
file open textfile_type_female using $results/sex_cox.csv, write append
file write textfile_type_female "Sepsis" ";" (`dementiaevents_1') ";" %4.1f (`persontime_1')  ";"
file write textfile_type_female %4.2f (`cruderate_1') "(" %4.2f (`lci_rate_1') "-" %4.2f (`uci_rate_1') ")" 
file write  textfile_type_female ";"  %4.2f (`hr_age_1') " (" %4.2f (`lci_age_1')  "-" %4.2f (`uci_age_1') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_min_1') " (" %4.2f (`lci_min_1')  "-" %4.2f (`uci_min_1') ")" 
file write  textfile_type_female ";"  %4.2f (`hr_ful_1') " (" %4.2f (`lci_ful_1')  "-" %4.2f (`uci_ful_1') ")" _n

file write textfile_type_female "Pneumonia" ";" (`dementiaevents_2') ";" %4.1f (`persontime_2')  ";"
file write textfile_type_female %4.2f (`cruderate_2') "(" %4.2f (`lci_rate_2') "-" %4.2f (`uci_rate_2') ")" 
file write  textfile_type_female ";"  %4.2f (`hr_age_2') " (" %4.2f (`lci_age_2')  "-" %4.2f (`uci_age_2') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_min_2') " (" %4.2f (`lci_min_2')  "-" %4.2f (`uci_min_2') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_ful_2') " (" %4.2f (`lci_ful_2')  "-" %4.2f (`uci_ful_2') ")" _n

file write textfile_type_female "Other LRTI" ";" (`dementiaevents_3') ";" %4.1f (`persontime_3')  ";"
file write textfile_type_female %4.2f (`cruderate_3') "(" %4.2f (`lci_rate_3') "-" %4.2f (`uci_rate_3') ")" 
file write  textfile_type_female ";"  %4.2f (`hr_age_3') " (" %4.2f (`lci_age_3')  "-" %4.2f (`uci_age_3') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_min_3') " (" %4.2f (`lci_min_3')  "-" %4.2f (`uci_min_3') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_ful_3') " (" %4.2f (`lci_ful_3')  "-" %4.2f (`uci_ful_3') ")" _n

file write textfile_type_female "UTI" ";" (`dementiaevents_4') ";" %4.1f (`persontime_4')  ";"
file write textfile_type_female %4.2f (`cruderate_4') "(" %4.2f (`lci_rate_4') "-" %4.2f (`uci_rate_4') ")" 
file write  textfile_type_female ";"  %4.2f (`hr_age_4') " (" %4.2f (`lci_age_4')  "-" %4.2f (`uci_age_4') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_min_4') " (" %4.2f (`lci_min_4')  "-" %4.2f (`uci_min_4') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_ful_4') " (" %4.2f (`lci_ful_4')  "-" %4.2f (`uci_ful_4') ")" _n

file write textfile_type_female "SSTI" ";" (`dementiaevents_5') ";" %4.1f (`persontime_5')  ";"
file write textfile_type_female %4.2f (`cruderate_5') "(" %4.2f (`lci_rate_5') "-" %4.2f (`uci_rate_5') ")" 
file write  textfile_type_female ";"  %4.2f (`hr_age_5') " (" %4.2f (`lci_age_5')  "-" %4.2f (`uci_age_5') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_min_5') " (" %4.2f (`lci_min_5')  "-" %4.2f (`uci_min_5') ")" 
file write  textfile_type_female ";"  %4.2f (`hr_ful_5') " (" %4.2f (`lci_ful_5')  "-" %4.2f (`uci_ful_5') ")" _n

file write textfile_type_female "Unknown" ";" (`dementiaevents_6') ";" %4.1f (`persontime_6')  ";"
file write textfile_type_female %4.2f (`cruderate_6') "(" %4.2f (`lci_rate_6') "-" %4.2f (`uci_rate_6') ")" 
file write  textfile_type_female ";"  %4.2f (`hr_age_6') " (" %4.2f (`lci_age_6')  "-" %4.2f (`uci_age_6') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_min_6') " (" %4.2f (`lci_min_6')  "-" %4.2f (`uci_min_6') ")"  
file write  textfile_type_female ";"  %4.2f (`hr_ful_6') " (" %4.2f (`lci_ful_6')  "-" %4.2f (`uci_ful_6') ")" _n

file close textfile_type_female


*Fully adjusted lrtest
stcox i.inf_type_tv  i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==1, base   
estimates store A
stcox  i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==1, base
estimates store B
lrtest a b //

log close



**********************************************************************************************************************
****************************************************MALE (sex==0)*****************************************************
**********************************************************************************************************************

*=================================ANY INFECTION=================================

use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log\male_cox.log, replace

capture file close textfile_male
file open textfile_male using $results/sex_cox.csv, write append

file write textfile_male "Male" _n

stptime if infection_tv==0 & sex==0,  title(person-years) per(1000)
return list
*no infection
local persontime_noinf_male=`r(ptime)'
local dementiaevents_noinf_male=`r(failures)'
local cruderate_noinf_male=`r(rate)'
local lci_rate_noinf_male=`r(lb)'
local uci_rate_noinf_male=`r(ub)'

stptime if infection_tv==1 & sex==0,  title(person-years) per(1000)
**infection
local persontime_inf_male=`r(ptime)'
local dementiaevents_inf_male=`r(failures)'
local cruderate_inf_male=`r(rate)'
local lci_rate_inf_male=`r(lb)'
local uci_rate_inf_male=`r(ub)'

file write textfile_male "No infection" ";" (`dementiaevents_noinf_male') ";" %3.1f (`persontime_noinf_male') ";"
file write textfile_male %3.2f (`cruderate_noinf_male') "(" %3.2f (`lci_rate_noinf_male') "-" %3.2f (`uci_rate_noinf_male') ")" _n

*age adjusted
stcox i.infection_tv  if sex==0, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.imd_cat i.calperiod  if sex==0, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==0, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_male "Any Infection" ";" (`dementiaevents_inf_male') ";" %3.1f (`persontime_inf_male')  ";"
file write textfile_male %3.2f (`cruderate_inf_male') "(" %3.2f (`lci_rate_inf_male') "-" %3.2f (`uci_rate_inf_male') ")"
file write  textfile_male ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_male ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_male ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_male

*full model
stcox i.infection_tv i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==0, base //
estimates store a
stcox i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==0, base
estimates store b
lrtest a b // P < 0.001


log close

*==============================INFECTION TYPE===================================
use $cleandata\mainanalysis_type, clear

cap log close
log using $log/typemale_cox.log, replace

foreach i of num 1/6 {
stptime if inf_type_tv==`i' & sex==0,  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local dementiaevents_`i'=`r(failures)'
local cruderate_`i'=`r(rate)'
local lci_rate_`i'=`r(lb)'
local uci_rate_`i'=`r(ub)'
}

*age adjusted
stcox i.inf_type_tv  if sex==0, base 
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

*minimally adjusted
stcox i.inf_type_tv  i.imd_cat i.calperiod  if sex==0, base
local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_min_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_min_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_min_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Fully adjusted 
stcox i.inf_type_tv  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==0, base 
foreach x of local infectionregression_var {
local hr_ful_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_ful_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_ful_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

capture file close textfile_type_male
*file close textfile
file open textfile_type_male using $results/sex_cox.csv, write append
file write textfile_type_male "Sepsis" ";" (`dementiaevents_1') ";" %4.1f (`persontime_1')  ";"
file write textfile_type_male %4.2f (`cruderate_1') "(" %4.2f (`lci_rate_1') "-" %4.2f (`uci_rate_1') ")" 
file write  textfile_type_male ";"  %4.2f (`hr_age_1') " (" %4.2f (`lci_age_1')  "-" %4.2f (`uci_age_1') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_min_1') " (" %4.2f (`lci_min_1')  "-" %4.2f (`uci_min_1') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_ful_1') " (" %4.2f (`lci_ful_1')  "-" %4.2f (`uci_ful_1') ")" _n

file write textfile_type_male "Pneumonia" ";" (`dementiaevents_2') ";" %4.1f (`persontime_2')  ";"
file write textfile_type_male %4.2f (`cruderate_2') "(" %4.2f (`lci_rate_2') "-" %4.2f (`uci_rate_2') ")" 
file write  textfile_type_male ";"  %4.2f (`hr_age_2') " (" %4.2f (`lci_age_2')  "-" %4.2f (`uci_age_2') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_min_2') " (" %4.2f (`lci_min_2')  "-" %4.2f (`uci_min_2') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_ful_2') " (" %4.2f (`lci_ful_2')  "-" %4.2f (`uci_ful_2') ")" _n

file write textfile_type_male "Other LRTI" ";" (`dementiaevents_3') ";" %4.1f (`persontime_3')  ";"
file write textfile_type_male %4.2f (`cruderate_3') "(" %4.2f (`lci_rate_3') "-" %4.2f (`uci_rate_3') ")" 
file write  textfile_type_male ";"  %4.2f (`hr_age_3') " (" %4.2f (`lci_age_3')  "-" %4.2f (`uci_age_3') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_min_3') " (" %4.2f (`lci_min_3')  "-" %4.2f (`uci_min_3') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_ful_3') " (" %4.2f (`lci_ful_3')  "-" %4.2f (`uci_ful_3') ")" _n

file write textfile_type_male "UTI" ";" (`dementiaevents_4') ";" %4.1f (`persontime_4')  ";"
file write textfile_type_male %4.2f (`cruderate_4') "(" %4.2f (`lci_rate_4') "-" %4.2f (`uci_rate_4') ")" 
file write  textfile_type_male ";"  %4.2f (`hr_age_4') " (" %4.2f (`lci_age_4')  "-" %4.2f (`uci_age_4') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_min_4') " (" %4.2f (`lci_min_4')  "-" %4.2f (`uci_min_4') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_ful_4') " (" %4.2f (`lci_ful_4')  "-" %4.2f (`uci_ful_4') ")" _n

file write textfile_type_male "SSTI" ";" (`dementiaevents_5') ";" %4.1f (`persontime_5')  ";"
file write textfile_type_male %4.2f (`cruderate_5') "(" %4.2f (`lci_rate_5') "-" %4.2f (`uci_rate_5') ")" 
file write  textfile_type_male ";"  %4.2f (`hr_age_5') " (" %4.2f (`lci_age_5')  "-" %4.2f (`uci_age_5') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_min_5') " (" %4.2f (`lci_min_5')  "-" %4.2f (`uci_min_5') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_ful_5') " (" %4.2f (`lci_ful_5')  "-" %4.2f (`uci_ful_5') ")" _n

file write textfile_type_male "Unknown" ";" (`dementiaevents_6') ";" %4.1f (`persontime_6')  ";"
file write textfile_type_male %4.2f (`cruderate_6') "(" %4.2f (`lci_rate_6') "-" %4.2f (`uci_rate_6') ")" 
file write  textfile_type_male ";"  %4.2f (`hr_age_6') " (" %4.2f (`lci_age_6')  "-" %4.2f (`uci_age_6') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_min_6') " (" %4.2f (`lci_min_6')  "-" %4.2f (`uci_min_6') ")"  
file write  textfile_type_male ";"  %4.2f (`hr_ful_6') " (" %4.2f (`lci_ful_6')  "-" %4.2f (`uci_ful_6') ")" _n

file close textfile_type_male


*Fully adjusted lrtest
stcox i.inf_type_tv  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==0, base   
estimates store A
stcox   i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if sex==0, base
estimates store B
lrtest a b //


