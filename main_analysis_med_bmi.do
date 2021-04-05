/*=========================================================================
DO FILE NAME:			Main_analyses1_mediators.do

AUTHOR:					Rutendo Muzambi


DATE VERSION CREATED: 	05/2020
						
DESCRIPTION OF FILE:    First ever infections main analyses
*=========================================================================*/
use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log\mediators_cox.log, replace
*st
strate infection_tv, per(1000) 
strate infection_tv, per(1000) output ($results/strate_med_bmi, replace)

*by infection_tv, sort:stptime, title(person-years) per(1000)

label define infection_tb 0 "No infection" 1"Infection" 
label values infection_tv infection_tv

*file close textfile

*tempname myfile
capture file close textfile 
file open textfile using $results/med_bmi.csv, write replace
file write textfile "sep=;" _n
file write textfile "Effect of infections on dementia" _n _n
file write textfile  ";"   "Number of events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" ";" "Additionally adjusted for mediators HR*** (95% CI)" ";" "Additionally adjusted for BMI HR*** (95% CI)" _n 

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

*Additionally adjusted for potential mediators*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat af_md mi_md stroke_md, base //

local hr_med=exp(_b[1.infection_tv])
local lci_med=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_med=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*Additionally adjusted for BMI*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat bmi_cat, base //

local hr_bmi=exp(_b[1.infection_tv])
local lci_bmi=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_bmi=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	


file write textfile "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" 
file write  textfile ";"  %3.2f (`hr_med') " (" %3.2f (`lci_med')  "-" %3.2f (`uci_med') ")"  
file write  textfile ";"  %3.2f (`hr_bmi') " (" %3.2f (`lci_bmi')  "-" %3.2f (`uci_bmi') ")" _n

file close textfile


*Fully adjusted model*
stcox i.infection_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001


log close

use $cleandata\mainanalysis_type, clear
log using $log/type_med_bmi_cox.log, replace

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
stcox i.inf_type_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod, base 

local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_min_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_min_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_min_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Fully adjusted model*
stcox i.inf_type_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //


local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_ful_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_ful_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_ful_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*Additionally adjusted for potential mediators*
stcox i.inf_type_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat af_md mi_md stroke_md, base //

local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_med_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_med_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_med_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}	

*Additionally adjusted for BMI*
stcox i.inf_type_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat bmi_cat, base //

local infectionregression_var 1 2 3 4 5 6 
foreach x of local infectionregression_var {
local hr_bmi_`x'=exp(_b[`x'.inf_type_tv]) 
local lci_bmi_`x'=exp(_b[`x'.inf_type_tv]-1.96*_se[`x'.inf_type_tv]) 
local uci_bmi_`x'=exp(_b[`x'.inf_type_tv]+1.96*_se[`x'.inf_type_tv]) 
}

*file open textfile using $results/firstever_cox.csv, write append
*file write textfile "sepsis" ";" (`dementiaevents_sepsis') ";" %4.1f (`persontime_sepsis') ";"

capture file close textfile 
*file close textfile
file open textfile_type using $results/med_bmi.csv, write append
file write textfile_type "Sepsis" ";" (`dementiaevents_1') ";" %3.1f (`persontime_1')  ";"
file write textfile_type %3.2f (`cruderate_1') "(" %3.2f (`lci_rate_1') "-" %3.2f (`uci_rate_1') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_1') " (" %3.2f (`lci_age_1')  "-" %3.2f (`uci_age_1') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_1') " (" %3.2f (`lci_min_1')  "-" %3.2f (`uci_min_1') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_1') " (" %3.2f (`lci_ful_1')  "-" %3.2f (`uci_ful_1') ")" 
file write  textfile_type ";"  %3.2f (`hr_med_1') " (" %3.2f (`lci_med_1')  "-" %3.2f (`uci_med_1') ")"  
file write  textfile_type ";"  %3.2f (`hr_bmi_1') " (" %3.2f (`lci_bmi_1')  "-" %3.2f (`uci_bmi_1') ")" _n


file write textfile_type "Pneumonia" ";" (`dementiaevents_2') ";" %3.1f (`persontime_2')  ";"
file write textfile_type %3.2f (`cruderate_2') "(" %3.2f (`lci_rate_2') "-" %3.2f (`uci_rate_2') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_2') " (" %3.2f (`lci_age_2')  "-" %3.2f (`uci_age_2') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_2') " (" %3.2f (`lci_min_2')  "-" %3.2f (`uci_min_2') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_2') " (" %3.2f (`lci_ful_2')  "-" %3.2f (`uci_ful_2') ")" 
file write  textfile_type ";"  %3.2f (`hr_med_2') " (" %3.2f (`lci_med_2')  "-" %3.2f (`uci_med_2') ")"  
file write  textfile_type ";"  %3.2f (`hr_bmi_2') " (" %3.2f (`lci_bmi_2')  "-" %3.2f (`uci_bmi_2') ")" _n

file write textfile_type "Other LRTI" ";" (`dementiaevents_3') ";" %3.1f (`persontime_3')  ";"
file write textfile_type %3.2f (`cruderate_3') "(" %3.2f (`lci_rate_3') "-" %3.2f (`uci_rate_3') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_3') " (" %3.2f (`lci_age_3')  "-" %3.2f (`uci_age_3') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_3') " (" %3.2f (`lci_min_3')  "-" %3.2f (`uci_min_3') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_3') " (" %3.2f (`lci_ful_3')  "-" %3.2f (`uci_ful_3') ")" 
file write  textfile_type ";"  %3.2f (`hr_med_3') " (" %3.2f (`lci_med_3')  "-" %3.2f (`uci_med_3') ")"  
file write  textfile_type ";"  %3.2f (`hr_bmi_3') " (" %3.2f (`lci_bmi_3')  "-" %3.2f (`uci_bmi_3') ")" _n

file write textfile_type "UTI" ";" (`dementiaevents_4') ";" %3.1f (`persontime_4')  ";"
file write textfile_type %3.2f (`cruderate_4') "(" %3.2f (`lci_rate_4') "-" %3.2f (`uci_rate_4') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_4') " (" %3.2f (`lci_age_4')  "-" %3.2f (`uci_age_4') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_4') " (" %3.2f (`lci_min_4')  "-" %3.2f (`uci_min_4') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_4') " (" %3.2f (`lci_ful_4')  "-" %3.2f (`uci_ful_4') ")" 
file write  textfile_type ";"  %3.2f (`hr_med_4') " (" %3.2f (`lci_med_4')  "-" %3.2f (`uci_med_4') ")"  
file write  textfile_type ";"  %3.2f (`hr_bmi_4') " (" %3.2f (`lci_bmi_4')  "-" %3.2f (`uci_bmi_4') ")" _n

file write textfile_type "SSTI" ";" (`dementiaevents_5') ";" %3.1f (`persontime_5')  ";"
file write textfile_type %3.2f (`cruderate_5') "(" %3.2f (`lci_rate_5') "-" %3.2f (`uci_rate_5') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_5') " (" %3.2f (`lci_age_5')  "-" %3.2f (`uci_age_5') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_5') " (" %3.2f (`lci_min_5')  "-" %3.2f (`uci_min_5') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_5') " (" %3.2f (`lci_ful_5')  "-" %3.2f (`uci_ful_5') ")" 
file write  textfile_type ";"  %3.2f (`hr_med_5') " (" %3.2f (`lci_med_5')  "-" %3.2f (`uci_med_5') ")"  
file write  textfile_type ";"  %3.2f (`hr_bmi_5') " (" %3.2f (`lci_bmi_5')  "-" %3.2f (`uci_bmi_5') ")" _n

file write textfile_type "Unknown" ";" (`dementiaevents_6') ";" %3.1f (`persontime_6')  ";"
file write textfile_type %3.2f (`cruderate_6') "(" %3.2f (`lci_rate_6') "-" %3.2f (`uci_rate_6') ")" 
file write  textfile_type ";"  %3.2f (`hr_age_6') " (" %3.2f (`lci_age_6')  "-" %3.2f (`uci_age_6') ")"  
file write  textfile_type ";"  %3.2f (`hr_min_6') " (" %3.2f (`lci_min_6')  "-" %3.2f (`uci_min_6') ")" 
file write  textfile_type ";"  %3.2f (`hr_ful_6') " (" %3.2f (`lci_ful_6')  "-" %3.2f (`uci_ful_6') ")" 
file write  textfile_type ";"  %3.2f (`hr_med_6') " (" %3.2f (`lci_med_6')  "-" %3.2f (`uci_med_6') ")"  
file write  textfile_type ";"  %3.2f (`hr_bmi_6') " (" %3.2f (`lci_bmi_6')  "-" %3.2f (`uci_bmi_6') ")" _n

file close textfile_type



*Fully adjusted model*
stcox i.inf_type_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b

log close











