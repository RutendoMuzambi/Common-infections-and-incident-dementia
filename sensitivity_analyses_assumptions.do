/*=========================================================================
DO FILE NAME:			Sensitivity_analyses_assumption.do

AUTHOR:					Rutendo Muzambi


DATE VERSION CREATED: 	07/2020
						
DESCRIPTION OF FILE:    Sensitivity analyses for variables that have violated
                          the cox regression assumption
*=========================================================================*/

use $cleandata\mainanalysis_firstever, clear
cap log close
*log using $log\assumption_cox.log, replace
*st

label define infection_tb 0 "No infection" 1"Infection" 
label values infection_tv infection_tv

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

list patid _t0 _t _d start_fup end calperiod in 1/30


stcox i.infection_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat , base //
estat phtest, d //
estat phtest, plot(1.infection_tv)

stphplot, by(infection_tv) adjust(sex imd_cat calperiod ethnicity_cat smokestatus alcohol_cat anx_dep_cat smi_cat ibd_cat multiple_sclerosis_cat RA_cat psoriasis_cat asthma_cat ckd_cat cld_cat cld_cat copd_cat diabetes_cat hypertension_cat obs_sleep_cat tbi_cat benzodiazepine_cat ppi_cat steroids_cat polypharmacy_cat heart_failure_cat stroke_cat)

*/
*======================Secondary analysis YOUNGER AGE===========================

use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log\younger_cox.log, replace

capture file close textfile_younger
file open textfile_younger using $results/secondary_assump_cox.csv, write replace
file write textfile_younger "sep=;" _n
file write textfile_younger "Effect of common infections on dementia stratified according to age" _n _n
file write textfile_younger  ";"   "Number of events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile_younger  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

strate sex infection_tv, per (1000) output ($results/younger, replace)

file write textfile_younger "Less than 80 years" _n

stsplit ageband, at (80)
label define ageband_lab 0"<80 years" 80">80 years"
label values ageband ageband_lab

***split by age
stptime if infection_tv==0 & ageband==0,  title(person-years) per(1000)
return list
*no infection
local persontime_no_younger=`r(ptime)'
local dementiaevents_no_younger=`r(failures)'
local cruderate_no_younger=`r(rate)'
local lci_rate_no_younger=`r(lb)'
local uci_rate_no_younger=`r(ub)'

stptime if infection_tv==1 & ageband==0,  title(person-years) per(1000)
**infection
local persontime_inf_younger=`r(ptime)'
local dementiaevents_inf_younger=`r(failures)'
local cruderate_inf_younger=`r(rate)'
local lci_rate_inf_younger=`r(lb)'
local uci_rate_inf_younger=`r(ub)'

file write textfile_younger "No infection" ";" (`dementiaevents_no_younger') ";" %3.1f (`persontime_no_younger') ";"
file write textfile_younger %3.2f (`cruderate_no_younger') "(" %3.2f (`lci_rate_no_younger') "-" %3.2f (`uci_rate_no_younger') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==0, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==0, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.heart_failure_cat i.stroke_cat if ageband==0, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_younger "Any Infection" ";" (`dementiaevents_inf_younger') ";" %3.1f (`persontime_inf_younger')  ";"
file write textfile_younger %3.2f (`cruderate_inf_younger') "(" %3.2f (`lci_rate_inf_younger') "-" %3.2f (`uci_rate_inf_younger') ")"
file write  textfile_younger ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_younger ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_younger ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_younger

*full model
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==0, base //
estimates store a
stcox  i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==0, base
estimates store b
lrtest a b // P < 0.001

log close


*=======================Secondary analysis OLDER AGE============================

use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log\older_cox.log, replace

stsplit ageband, at (80)
label define ageband_lab 0"<80 years" 80">80 years"
label values ageband ageband_lab

capture file close textfile_older
file open textfile_older using $results/secondary_assump_cox.csv, write append

file write textfile_older "80 years and older" _n

stptime if infection_tv==0 & ageband==80,  title(person-years) per(1000)
return list
*no infection
local persontime_noinf_older=`r(ptime)'
local dementiaevents_noinf_older=`r(failures)'
local cruderate_noinf_older=`r(rate)'
local lci_rate_noinf_older=`r(lb)'
local uci_rate_noinf_older=`r(ub)'

stptime if infection_tv==1 & ageband==80,  title(person-years) per(1000)
**infection
local persontime_inf_older=`r(ptime)'
local dementiaevents_inf_older=`r(failures)'
local cruderate_inf_older=`r(rate)'
local lci_rate_inf_older=`r(lb)'
local uci_rate_inf_older=`r(ub)'

file write textfile_older "No infection" ";" (`dementiaevents_noinf_older') ";" %3.1f (`persontime_noinf_older') ";"
file write textfile_older %3.2f (`cruderate_noinf_older') "(" %3.2f (`lci_rate_noinf_older') "-" %3.2f (`uci_rate_noinf_older') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==80, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==80, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==80, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_older "Any Infection" ";" (`dementiaevents_inf_older') ";" %3.1f (`persontime_inf_older')  ";"
file write textfile_older %3.2f (`cruderate_inf_older') "(" %3.2f (`lci_rate_inf_older') "-" %3.2f (`uci_rate_inf_older') ")"
file write  textfile_older ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_older ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_older ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_older

*full model
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==80, base //
estimates store a
stcox  i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==80, base
estimates store b
lrtest a b // P < 0.001


log close

*/
********************************************************************************
******************************SENSITIVITY ANALYSES******************************
********************************************************************************

use $cleandata\mainanalysis_firstever, clear

stsplit ageband, at (80,90)
label define ageband_lab 0"<80 years" 80"80-89" 90"90+"
label values ageband ageband_lab

*sensitivity

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab


stcox i.infection_tv i.sex##ageband  i.imd_cat##ageband i.calperiod##ageband i.ethnicity_cat##ageband i.smokestatus##ageband i.alcohol_cat##ageband i.anx_dep_cat##ageband i.smi_cat##ageband i.ibd_cat##ageband i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat##ageband i.cld_cat i.copd_cat i.diabetes_cat##ageband i.hypertension_cat##ageband i.obs_sleep_cat i.tbi_cat##ageband i.benzodiazepine_cat##ageband i.ppi_cat i.steroids_cat##ageband i.polypharmacy_cat##ageband i.heart_failure_cat##ageband i.stroke_cat##ageband, base //
estat phtest, d //infection_cat P=0.08, HR 1.21 (95% CI 1.19-1.23)


strate infection_tv, per(1000) 
strate infection_tv, per(1000) output ($results/sens_assumpt, replace)

*by infection_tv, sort:stptime, title(person-years) per(1000)

label define infection_tb 0 "No infection" 1"Infection" 
label values infection_tv infection_tv

*file close textfile

*tempname myfile
capture file close textfile 
file open textfile using $results/sens_assump_cox.csv, write replace
file write textfile "sep=;" _n
file write textfile "Effect of common infections on dementia with interactions between time and covariates" _n _n
file write textfile  ";"   "Number of dementia events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
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

**Minimally adjusted model 
stcox i.infection_tv i.sex##ageband i.imd_cat##ageband i.calperiod##ageband, base 
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

*Fully adjusted model*
stcox i.infection_tv i.sex##ageband  i.imd_cat##ageband i.calperiod##ageband i.ethnicity_cat##ageband i.smokestatus##ageband i.alcohol_cat##ageband i.anx_dep_cat##ageband i.smi_cat##ageband i.ibd_cat##ageband i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat##ageband i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat##ageband i.obs_sleep_cat i.tbi_cat##ageband i.benzodiazepine_cat##ageband i.ppi_cat i.steroids_cat##ageband i.polypharmacy_cat##ageband i.heart_failure_cat##ageband i.stroke_cat##ageband, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile



*======================Secondary analysis YOUNGER AGE===========================

use $cleandata\mainanalysis_firstever, clear
cap log close

capture file close textfile_younger2
file open textfile_younger2 using $results/secondary_assump_cox.csv, write replace
file write textfile_younger2 "sep=;" _n
file write textfile_younger2 "Effect of common infections on dementia, stratified by age" _n _n
file write textfile_younger2  ";"   "Number of dementia events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile_younger2  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

file write textfile_younger2 "Less than 90 years" _n

stsplit ageband, at (90)
label define ageband_lab 0"<90 years" 90">90 years"
label values ageband ageband_lab

***split by age
stptime if infection_tv==0 & ageband==0,  title(person-years) per(1000)
return list
*no infection
local persontime_no_younger2=`r(ptime)'
local dementiaevents_no_younger2=`r(failures)'
local cruderate_no_younger2=`r(rate)'
local lci_rate_no_younger2=`r(lb)'
local uci_rate_no_younger2=`r(ub)'

stptime if infection_tv==1 & ageband==0,  title(person-years) per(1000)
**infection
local persontime_inf_younger2=`r(ptime)'
local dementiaevents_inf_younger2=`r(failures)'
local cruderate_inf_younger2=`r(rate)'
local lci_rate_inf_younger2=`r(lb)'
local uci_rate_inf_younger2=`r(ub)'

file write textfile_younger2 "No infection" ";" (`dementiaevents_no_younger2') ";" %3.1f (`persontime_no_younger2') ";"
file write textfile_younger2 %3.2f (`cruderate_no_younger2') "(" %3.2f (`lci_rate_no_younger2') "-" %3.2f (`uci_rate_no_younger2') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==0, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==0, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat if ageband==0, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_younger2 "Any Infection" ";" (`dementiaevents_inf_younger2') ";" %3.1f (`persontime_inf_younger2')  ";"
file write textfile_younger2 %3.2f (`cruderate_inf_younger2') "(" %3.2f (`lci_rate_inf_younger2') "-" %3.2f (`uci_rate_inf_younger2') ")"
file write  textfile_younger2 ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_younger2 ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_younger2 ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_younger2

/*
*full model
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==0, base //
estimates store a
stcox  i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==0, base
estimates store b
lrtest a b // P < 0.001
*/

*=======================Secondary analysis OLDER AGE============================

use $cleandata\mainanalysis_firstever, clear
cap log close

stsplit ageband, at (90)
label define ageband_lab 0"<90 years" 90">90 years"
label values ageband ageband_lab

capture file close textfile_younger2
file open textfile_older2 using $results/secondary_assump_cox.csv, write append

file write textfile_older2 "90 years and older2" _n
stptime if infection_tv==0 & ageband==90,  title(person-years) per(1000)
return list
*no infection
local persontime_noinf_older2=`r(ptime)'
local dementiaevents_noinf_older2=`r(failures)'
local cruderate_noinf_older2=`r(rate)'
local lci_rate_noinf_older2=`r(lb)'
local uci_rate_noinf_older2=`r(ub)'

stptime if infection_tv==1 & ageband==90,  title(person-years) per(1000)
**infection
local persontime_inf_older2=`r(ptime)'
local dementiaevents_inf_older2=`r(failures)'
local cruderate_inf_older2=`r(rate)'
local lci_rate_inf_older2=`r(lb)'
local uci_rate_inf_older2=`r(ub)'

file write textfile_older2 "No infection" ";" (`dementiaevents_noinf_older2') ";" %3.1f (`persontime_noinf_older2') ";"
file write textfile_older2 %3.2f (`cruderate_noinf_older2') "(" %3.2f (`lci_rate_noinf_older2') "-" %3.2f (`uci_rate_noinf_older2') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==90, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==90, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==90, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_older2 "Any Infection" ";" (`dementiaevents_inf_older2') ";" %3.1f (`persontime_inf_older2')  ";"
file write textfile_older2 %3.2f (`cruderate_inf_older2') "(" %3.2f (`lci_rate_inf_older2') "-" %3.2f (`uci_rate_inf_older2') ")"
file write  textfile_older2 ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_older2 ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_older2 ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_older2

/*full model
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==90, base //
estimates store a
stcox  i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==90, base
estimates store b
lrtest a b // P < 0.001
*/

********************************************************************************
****************80 cut off
***********************************************************************************


*======================Secondary analysis YOUNGER AGE===========================

use $cleandata\mainanalysis_firstever, clear
cap log close

capture file close textfile_older
file open textfile_younger2 using $results/secondary_assump_cox.csv, write append

file write textfile_younger2 "Less than 80 years" _n

stsplit ageband, at (80)
label define ageband_lab 0"<80 years" 80">80 years"
label values ageband ageband_lab

***split by age
stptime if infection_tv==0 & ageband==0,  title(person-years) per(1000)
return list
*no infection
local persontime_no_younger2=`r(ptime)'
local dementiaevents_no_younger2=`r(failures)'
local cruderate_no_younger2=`r(rate)'
local lci_rate_no_younger2=`r(lb)'
local uci_rate_no_younger2=`r(ub)'

stptime if infection_tv==1 & ageband==0,  title(person-years) per(1000)
**infection
local persontime_inf_younger2=`r(ptime)'
local dementiaevents_inf_younger2=`r(failures)'
local cruderate_inf_younger2=`r(rate)'
local lci_rate_inf_younger2=`r(lb)'
local uci_rate_inf_younger2=`r(ub)'

file write textfile_younger2 "No infection" ";" (`dementiaevents_no_younger2') ";" %3.1f (`persontime_no_younger2') ";"
file write textfile_younger2 %3.2f (`cruderate_no_younger2') "(" %3.2f (`lci_rate_no_younger2') "-" %3.2f (`uci_rate_no_younger2') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==0, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==0, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat if ageband==0, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_younger2 "Any Infection" ";" (`dementiaevents_inf_younger2') ";" %3.1f (`persontime_inf_younger2')  ";"
file write textfile_younger2 %3.2f (`cruderate_inf_younger2') "(" %3.2f (`lci_rate_inf_younger2') "-" %3.2f (`uci_rate_inf_younger2') ")"
file write  textfile_younger2 ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_younger2 ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_younger2 ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_younger2

/*
*full model
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==0, base //
estimates store a
stcox  i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==0, base
estimates store b
lrtest a b // P < 0.001
*/


*=======================Secondary analysis OLDER AGE============================

use $cleandata\mainanalysis_firstever, clear
cap log close

stsplit ageband, at (80)
label define ageband_lab 0"<80 years" 80">80 years"
label values ageband ageband_lab

capture file close textfile_younger2
file open textfile_older2 using $results/secondary_assump_cox.csv, write append

file write textfile_older2 "80 years and older2" _n
stptime if infection_tv==0 & ageband==80,  title(person-years) per(1000)
return list
*no infection
local persontime_noinf_older2=`r(ptime)'
local dementiaevents_noinf_older2=`r(failures)'
local cruderate_noinf_older2=`r(rate)'
local lci_rate_noinf_older2=`r(lb)'
local uci_rate_noinf_older2=`r(ub)'

stptime if infection_tv==1 & ageband==80,  title(person-years) per(1000)
**infection
local persontime_inf_older2=`r(ptime)'
local dementiaevents_inf_older2=`r(failures)'
local cruderate_inf_older2=`r(rate)'
local lci_rate_inf_older2=`r(lb)'
local uci_rate_inf_older2=`r(ub)'

file write textfile_older2 "No infection" ";" (`dementiaevents_noinf_older2') ";" %3.1f (`persontime_noinf_older2') ";"
file write textfile_older2 %3.2f (`cruderate_noinf_older2') "(" %3.2f (`lci_rate_noinf_older2') "-" %3.2f (`uci_rate_noinf_older2') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==80, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==80, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==80, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_older2 "Any Infection" ";" (`dementiaevents_inf_older2') ";" %3.1f (`persontime_inf_older2')  ";"
file write textfile_older2 %3.2f (`cruderate_inf_older2') "(" %3.2f (`lci_rate_inf_older2') "-" %3.2f (`uci_rate_inf_older2') ")"
file write  textfile_older2 ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_older2 ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_older2 ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_older2


********************************************************************************
****************3 category cut off 65-79
********************************************************************************

use $cleandata\mainanalysis_firstever, clear
cap log close

capture file close textfile_65

file open textfile_65 using $results/assump2_cox.csv, write replace
file write textfile_65 "sep=;" _n
file write textfile_65 "Effect of common infections on dementia stratified according to age categories" _n _n
file write textfile_65  ";"   "Number of dementia events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile_65  ";"   "Age-adjusted HR* (90% CI)" ";" "Minimally-adjusted HR** (90% CI)" ";" "Fully adjusted HR*** (90% CI)" _n 

stsplit ageband, at (80,90)
label define ageband_lab 0"65-79" 80"80-89" 90"90+"
label values ageband ageband_lab



file write textfile_65 "65-79" _n
stptime if infection_tv==0 & ageband==0,  title(person-years) per(1000)
return list
*no infection
local persontime_noinf_older2=`r(ptime)'
local dementiaevents_noinf_older2=`r(failures)'
local cruderate_noinf_older2=`r(rate)'
local lci_rate_noinf_older2=`r(lb)'
local uci_rate_noinf_older2=`r(ub)'

stptime if infection_tv==1 & ageband==0,  title(person-years) per(1000)
**infection
local persontime_inf_older2=`r(ptime)'
local dementiaevents_inf_older2=`r(failures)'
local cruderate_inf_older2=`r(rate)'
local lci_rate_inf_older2=`r(lb)'
local uci_rate_inf_older2=`r(ub)'

file write textfile_65 "No infection" ";" (`dementiaevents_noinf_older2') ";" %3.1f (`persontime_noinf_older2') ";"
file write textfile_65 %3.2f (`cruderate_noinf_older2') "(" %3.2f (`lci_rate_noinf_older2') "-" %3.2f (`uci_rate_noinf_older2') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==0, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==0, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==0, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_65 "Any Infection" ";" (`dementiaevents_inf_older2') ";" %3.1f (`persontime_inf_older2')  ";"
file write textfile_65 %3.2f (`cruderate_inf_older2') "(" %3.2f (`lci_rate_inf_older2') "-" %3.2f (`uci_rate_inf_older2') ")"
file write  textfile_65 ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_65 ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_65 ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_65

********************************************************************************
****************3 category cut off 80-94
********************************************************************************

use $cleandata\mainanalysis_firstever, clear
cap log close

stsplit ageband, at (80,90)
label define ageband_lab 0"65-79" 80"80-94" 90"90+"
label values ageband ageband_lab

capture file close textfile_younger2
file open textfile_80 using $results/assump2_cox.csv, write append

file write textfile_80 "80-94" _n
stptime if infection_tv==0 & ageband==80,  title(person-years) per(1000)
return list
*no infection
local persontime_noinf_older2=`r(ptime)'
local dementiaevents_noinf_older2=`r(failures)'
local cruderate_noinf_older2=`r(rate)'
local lci_rate_noinf_older2=`r(lb)'
local uci_rate_noinf_older2=`r(ub)'

stptime if infection_tv==1 & ageband==80,  title(person-years) per(1000)
**infection
local persontime_inf_older2=`r(ptime)'
local dementiaevents_inf_older2=`r(failures)'
local cruderate_inf_older2=`r(rate)'
local lci_rate_inf_older2=`r(lb)'
local uci_rate_inf_older2=`r(ub)'

file write textfile_80 "No infection" ";" (`dementiaevents_noinf_older2') ";" %3.1f (`persontime_noinf_older2') ";"
file write textfile_80 %3.2f (`cruderate_noinf_older2') "(" %3.2f (`lci_rate_noinf_older2') "-" %3.2f (`uci_rate_noinf_older2') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==80, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==80, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==80, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_80 "Any Infection" ";" (`dementiaevents_inf_older2') ";" %3.1f (`persontime_inf_older2')  ";"
file write textfile_80 %3.2f (`cruderate_inf_older2') "(" %3.2f (`lci_rate_inf_older2') "-" %3.2f (`uci_rate_inf_older2') ")"
file write  textfile_80 ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_80 ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_80 ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_80



********************************************************************************
****************3 category cut off 90+
********************************************************************************

use $cleandata\mainanalysis_firstever, clear
cap log close

stsplit ageband, at (80,90)
label define ageband_lab 0"65-79" 80"80-94" 90"90+"
label values ageband ageband_lab

capture file close textfile_younger2
file open textfile_90 using $results/assump2_cox.csv, write append

file write textfile_90 "90 years and older" _n
stptime if infection_tv==0 & ageband==90,  title(person-years) per(1000)
return list
*no infection
local persontime_noinf_older3=`r(ptime)'
local dementiaevents_noinf_older3=`r(failures)'
local cruderate_noinf_older3=`r(rate)'
local lci_rate_noinf_older3=`r(lb)'
local uci_rate_noinf_older3=`r(ub)'

stptime if infection_tv==1 & ageband==90,  title(person-years) per(1000)
**infection
local persontime_inf_older3=`r(ptime)'
local dementiaevents_inf_older3=`r(failures)'
local cruderate_inf_older3=`r(rate)'
local lci_rate_inf_older3=`r(lb)'
local uci_rate_inf_older3=`r(ub)'

file write textfile_90 "No infection" ";" (`dementiaevents_noinf_older3') ";" %3.1f (`persontime_noinf_older3') ";"
file write textfile_90 %3.2f (`cruderate_noinf_older3') "(" %3.2f (`lci_rate_noinf_older3') "-" %3.2f (`uci_rate_noinf_older3') ")" _n

*age adjusted
stcox i.infection_tv  if ageband==90, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*calender period

stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if ageband==90, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==90, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile_90 "Any Infection" ";" (`dementiaevents_inf_older3') ";" %3.1f (`persontime_inf_older3')  ";"
file write textfile_90 %3.2f (`cruderate_inf_older3') "(" %3.2f (`lci_rate_inf_older3') "-" %3.2f (`uci_rate_inf_older3') ")"
file write  textfile_90 ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_90 ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_90 ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_90


stcox i.infection_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==0, base //
estat phtest, d //

stcox i.infection_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==80, base //
estat phtest, d //

stcox i.infection_tv i.sex i.ethnicity_cat i.imd_cat i.calperiod i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if ageband==90, base //
estat phtest, d //






