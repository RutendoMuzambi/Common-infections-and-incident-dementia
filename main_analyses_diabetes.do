/*=========================================================================
DO FILE NAME:			Main_analyses_diabetes.do

AUTHOR:					Rutendo Muzambi

VERSION:				v2

DATE VERSION CREATED: 	05/2020
						
DESCRIPTION OF FILE:    Type of infections on dementia
*=========================================================================*/


**********************************************************************************************************************
************************************************Diabetes Mellitus*****************************************************
**********************************************************************************************************************

*=================================ANY INFECTION=================================

use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log\diabetes_cox.log, replace

capture file close textfile_dm
file open textfile_dm using $results/diabetes_cox.csv, write replace
file write textfile_dm "sep=;" _n
file write textfile_dm "Effect of infections on dementia in people with and without diabetes mellitus" _n _n
file write textfile_dm  ";"   "Number of dementia events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile_dm  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

strate diabetes_cat infection_tv, per (1000) output ($results/diabetes, replace)

file write textfile_dm "Individuals with diabetes mellitus" _n
file write textfile_dm "Type of infection" _n

****diabetes
stptime if infection_tv==0 & diabetes_cat==1,  title(person-years) per(1000)
return list
*no infection
local persontime_noinf_dm=`r(ptime)'
local dementiaevents_noinf_dm=`r(failures)'
local cruderate_noinf_dm=`r(rate)'
local lci_rate_noinf_dm=`r(lb)'
local uci_rate_noinf_dm=`r(ub)'

stptime if infection_tv==1 & diabetes_cat==1,  title(person-years) per(1000)
**infection
local persontime_inf_dm=`r(ptime)'
local dementiaevents_inf_dm=`r(failures)'
local cruderate_inf_dm=`r(rate)'
local lci_rate_inf_dm=`r(lb)'
local uci_rate_inf_dm=`r(ub)'

file write textfile_dm "No infection" ";" (`dementiaevents_noinf_dm') ";" %3.1f (`persontime_noinf_dm') ";"
file write textfile_dm %3.2f (`cruderate_noinf_dm') "(" %3.2f (`lci_rate_noinf_dm') "-" %3.2f (`uci_rate_noinf_dm') ")" _n

*age adjusted
stcox i.infection_tv  if diabetes_cat==1, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if diabetes_cat==1, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==1, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_dm "Any Infection" ";" (`dementiaevents_inf_dm') ";" %3.1f (`persontime_inf_dm')  ";"
file write textfile_dm %3.2f (`cruderate_inf_dm') "(" %3.2f (`lci_rate_inf_dm') "-" %3.2f (`uci_rate_inf_dm') ")"
file write  textfile_dm ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_dm ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_dm ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_dm

*full model
stcox i.infection_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==1, base //
estimates store a
stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==1, base
estimates store b
lrtest a b // P < 0.001

***Effect modification diabetes fully adjusted****
stcox i.infection_tv##i.diabetes_cat i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base 
estimates store a
stcox i.infection_tv i.diabetes_cat i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base 
estimates store b
lrtest a b  //


log close

*==============================INFECTION TYPE===================================
use $cleandata\mainanalysis_type, clear

cap log close
log using $log/typedm_cox.log, replace
strate diabetes_cat inf_type_tv, per (1000) output ($results/diabetes_type, replace)

foreach i of num 1/6 {
stptime if inf_type_tv==`i' & diabetes_cat==1,  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local dementiaevents_`i'=`r(failures)'
local cruderate_`i'=`r(rate)'
local lci_rate_`i'=`r(lb)'
local uci_rate_`i'=`r(ub)'
}

*age adjusted
stcox i.inf_type_tv  if diabetes_cat==1, base 
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
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod  if diabetes_cat==1, base
local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_min_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_min_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_min_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Fully adjusted 
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==1, base 
foreach x of local infectionregression_var {
local hr_ful_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_ful_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_ful_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

capture file close textfile_type_dm
*file close textfile
file open textfile_type_dm using $results/diabetes_cox.csv, write append
file write textfile_type_dm "Sepsis" ";" (`dementiaevents_1') ";" %4.1f (`persontime_1')  ";"
file write textfile_type_dm %4.2f (`cruderate_1') "(" %4.2f (`lci_rate_1') "-" %4.2f (`uci_rate_1') ")" 
file write  textfile_type_dm ";"  %4.2f (`hr_age_1') " (" %4.2f (`lci_age_1')  "-" %4.2f (`uci_age_1') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_min_1') " (" %4.2f (`lci_min_1')  "-" %4.2f (`uci_min_1') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_ful_1') " (" %4.2f (`lci_ful_1')  "-" %4.2f (`uci_ful_1') ")" _n

file write textfile_type_dm "Pneumonia" ";" (`dementiaevents_2') ";" %4.1f (`persontime_2')  ";"
file write textfile_type_dm %4.2f (`cruderate_2') "(" %4.2f (`lci_rate_2') "-" %4.2f (`uci_rate_2') ")" 
file write  textfile_type_dm ";"  %4.2f (`hr_age_2') " (" %4.2f (`lci_age_2')  "-" %4.2f (`uci_age_2') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_min_2') " (" %4.2f (`lci_min_2')  "-" %4.2f (`uci_min_2') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_ful_2') " (" %4.2f (`lci_ful_2')  "-" %4.2f (`uci_ful_2') ")" _n

file write textfile_type_dm "Other LRTI" ";" (`dementiaevents_3') ";" %4.1f (`persontime_3')  ";"
file write textfile_type_dm %4.2f (`cruderate_3') "(" %4.2f (`lci_rate_3') "-" %4.2f (`uci_rate_3') ")" 
file write  textfile_type_dm ";"  %4.2f (`hr_age_3') " (" %4.2f (`lci_age_3')  "-" %4.2f (`uci_age_3') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_min_3') " (" %4.2f (`lci_min_3')  "-" %4.2f (`uci_min_3') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_ful_3') " (" %4.2f (`lci_ful_3')  "-" %4.2f (`uci_ful_3') ")" _n

file write textfile_type_dm "UTI" ";" (`dementiaevents_4') ";" %4.1f (`persontime_4')  ";"
file write textfile_type_dm %4.2f (`cruderate_4') "(" %4.2f (`lci_rate_4') "-" %4.2f (`uci_rate_4') ")" 
file write  textfile_type_dm ";"  %4.2f (`hr_age_4') " (" %4.2f (`lci_age_4')  "-" %4.2f (`uci_age_4') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_min_4') " (" %4.2f (`lci_min_4')  "-" %4.2f (`uci_min_4') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_ful_4') " (" %4.2f (`lci_ful_4')  "-" %4.2f (`uci_ful_4') ")" _n

file write textfile_type_dm "SSTI" ";" (`dementiaevents_5') ";" %4.1f (`persontime_5')  ";"
file write textfile_type_dm %4.2f (`cruderate_5') "(" %4.2f (`lci_rate_5') "-" %4.2f (`uci_rate_5') ")" 
file write  textfile_type_dm ";"  %4.2f (`hr_age_5') " (" %4.2f (`lci_age_5')  "-" %4.2f (`uci_age_5') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_min_5') " (" %4.2f (`lci_min_5')  "-" %4.2f (`uci_min_5') ")"  
file write  textfile_type_dm ";"  %4.2f (`hr_ful_5') " (" %4.2f (`lci_ful_5')  "-" %4.2f (`uci_ful_5') ")" _n

file close textfile_type_dm


*Fully adjusted lrtest
stcox i.inf_type_tv i.sex  i.imd_cat i.calperiod i.smokestatus i.ethnicity_cat i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==1, base   
estimates store A
stcox i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==1, base
estimates store B
lrtest a b //<0.024

*Effect modification diabetes

stcox i.inf_type_tv##i.diabetes_cat i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base   
estimates store a
stcox i.inf_type_tv i.diabetes_cat i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base
estimates store b
lrtest a b //

log close



**********************************************************************************************************************
*********************************************NO Diabetes Mellitus*****************************************************
**********************************************************************************************************************

*=================================ANY INFECTION=================================

use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log\diabetes_cox0.log, replace

capture file close textfile
file open textfile using $results/diabetes_cox.csv, write append

strate diabetes_cat infection_tv, per (1000) output ($results/diabetes, replace)

file write textfile "Individuals without diabetes mellitus" _n

****diabetes
stptime if infection_tv==0 & diabetes_cat==0,  title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local dementiaevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

stptime if infection_tv==1 & diabetes_cat==0,  title(person-years) per(1000)
**infection
local persontime_inf=`r(ptime)'
local dementiaevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

file write textfile "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n

*age adjusted
stcox i.infection_tv  if diabetes_cat==0, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if diabetes_cat==0, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==0, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile

*full model
stcox i.infection_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==0, base //
estimates store a
stcox i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==0, base
estimates store b
lrtest a b // P < 0.001

log close

*==============================INFECTION TYPE===================================
use $cleandata\mainanalysis_type, clear

cap log close
log using $log\diabetes_cox0_type.log, replace
strate diabetes_cat inf_type_tv, per (1000) output ($results/diabetes_type, replace)

foreach i of num 1/6 {
stptime if inf_type_tv==`i' & diabetes_cat==0,  title(person-years) per(1000)
local persontime_`i'=`r(ptime)'
local dementiaevents_`i'=`r(failures)'
local cruderate_`i'=`r(rate)'
local lci_rate_`i'=`r(lb)'
local uci_rate_`i'=`r(ub)'
}

*age adjusted
stcox i.inf_type_tv  if diabetes_cat==0, base 
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
stcox i.inf_type_tv i.sex i.imd_cat i.calperiod  if diabetes_cat==0, base
local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_min_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_min_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_min_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Fully adjusted 
stcox i.inf_type_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if diabetes_cat==0, base 
foreach x of local infectionregression_var {
local hr_ful_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_ful_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_ful_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

capture file close textfile_type
*file close textfile
file open textfile_type using $results/diabetes_cox.csv, write append
file write textfile_type "Sepsis" ";" (`dementiaevents_1') ";" %4.1f (`persontime_1')  ";"
file write textfile_type %4.2f (`cruderate_1') "(" %4.2f (`lci_rate_1') "-" %4.2f (`uci_rate_1') ")" 
file write  textfile_type ";"  %4.2f (`hr_age_1') " (" %4.2f (`lci_age_1')  "-" %4.2f (`uci_age_1') ")"  
file write  textfile_type ";"  %4.2f (`hr_min_1') " (" %4.2f (`lci_min_1')  "-" %4.2f (`uci_min_1') ")" 
file write  textfile_type ";"  %4.2f (`hr_ful_1') " (" %4.2f (`lci_ful_1')  "-" %4.2f (`uci_ful_1') ")" _n

file write textfile_type "Pneumonia" ";" (`dementiaevents_2') ";" %4.1f (`persontime_2')  ";"
file write textfile_type %4.2f (`cruderate_2') "(" %4.2f (`lci_rate_2') "-" %4.2f (`uci_rate_2') ")" 
file write  textfile_type ";"  %4.2f (`hr_age_2') " (" %4.2f (`lci_age_2')  "-" %4.2f (`uci_age_2') ")"  
file write  textfile_type ";"  %4.2f (`hr_min_2') " (" %4.2f (`lci_min_2')  "-" %4.2f (`uci_min_2') ")"  
file write  textfile_type ";"  %4.2f (`hr_ful_2') " (" %4.2f (`lci_ful_2')  "-" %4.2f (`uci_ful_2') ")" _n

file write textfile_type "Other LRTI" ";" (`dementiaevents_3') ";" %4.1f (`persontime_3')  ";"
file write textfile_type %4.2f (`cruderate_3') "(" %4.2f (`lci_rate_3') "-" %4.2f (`uci_rate_3') ")" 
file write  textfile_type ";"  %4.2f (`hr_age_3') " (" %4.2f (`lci_age_3')  "-" %4.2f (`uci_age_3') ")"  
file write  textfile_type ";"  %4.2f (`hr_min_3') " (" %4.2f (`lci_min_3')  "-" %4.2f (`uci_min_3') ")"  
file write  textfile_type ";"  %4.2f (`hr_ful_3') " (" %4.2f (`lci_ful_3')  "-" %4.2f (`uci_ful_3') ")" _n

file write textfile_type "UTI" ";" (`dementiaevents_4') ";" %4.1f (`persontime_4')  ";"
file write textfile_type %4.2f (`cruderate_4') "(" %4.2f (`lci_rate_4') "-" %4.2f (`uci_rate_4') ")" 
file write  textfile_type ";"  %4.2f (`hr_age_4') " (" %4.2f (`lci_age_4')  "-" %4.2f (`uci_age_4') ")"  
file write  textfile_type ";"  %4.2f (`hr_min_4') " (" %4.2f (`lci_min_4')  "-" %4.2f (`uci_min_4') ")"  
file write  textfile_type ";"  %4.2f (`hr_ful_4') " (" %4.2f (`lci_ful_4')  "-" %4.2f (`uci_ful_4') ")" _n

file write textfile_type "SSTI" ";" (`dementiaevents_5') ";" %4.1f (`persontime_5')  ";"
file write textfile_type %4.2f (`cruderate_5') "(" %4.2f (`lci_rate_5') "-" %4.2f (`uci_rate_5') ")" 
file write  textfile_type ";"  %4.2f (`hr_age_5') " (" %4.2f (`lci_age_5')  "-" %4.2f (`uci_age_5') ")"  
file write  textfile_type ";"  %4.2f (`hr_min_5') " (" %4.2f (`lci_min_5')  "-" %4.2f (`uci_min_5') ")"  
file write  textfile_type ";"  %4.2f (`hr_ful_5') " (" %4.2f (`lci_ful_5')  "-" %4.2f (`uci_ful_5') ")" _n


file close textfile_type


*Fully adjusted lrtest
stcox i.inf_type_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if diabetes_cat==0, base   
estimates store A
stcox i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat if diabetes_cat==0, base
estimates store B
lrtest a b //

log close
