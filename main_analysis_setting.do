*****************************
***Setting
*****************************

*===========================1. GP recorded infections=====================

use $cleandata\mainanalysis_gp_cox, replace
replace infection_date1= end+365.25 if infection_date1==. //infection date after end date
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

strate infection_tv, per(1000) output ($results/strate_gp, replace)
capture file close textfile 
file open textfile using $results/table3_cox.csv, write replace
file write textfile "sep=;" _n
file write textfile "Crude incidence rate and hazard ratios for the association between common infections and dementia, by GP vs HES recorded infections, diabetes status and frequency of infections " _n _n
file write textfile  ";"   "Number of events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

*incidence rate*
file write textfile "GP recorded infections" _n
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
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile

*===========================2. Hospital recorded infections=====================

use $cleandata\mainanalysis_hospital_cox, replace

replace infection_date1= end+365.25 if infection_date1==. //infection date after end date
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

strate infection_tv, per(1000) output ($results/strate_hosp, replace)

capture file close textfile_hosp 
file open textfile_hosp using $results/table3_cox.csv, write append

*incidence rate*
file write textfile_hosp "Hospital recorded infections" _n
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

file write textfile_hosp "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile_hosp %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n

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
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile_hosp "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile_hosp %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile_hosp ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_hosp ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_hosp ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_hosp

*===========================3. Primary Hospital recorded infections=====================

use $cleandata\mainanalysis_hospital_cox_primary, clear

replace infection_date1= end+365.25 if infection_date1==. //infection date after end date
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

strate infection_tv, per(1000) output ($results/strate_hosp_prim, replace)

capture file close textfile_hosp_prim 
file open textfile_hosp_prim using $results/table3_cox.csv, write append

*incidence rate*
file write textfile_hosp_prim "Primary hospital diagnosis" _n
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

file write textfile_hosp_prim "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile_hosp_prim %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n

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
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //

local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

file write textfile_hosp_prim "Any Infection" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile_hosp_prim %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile_hosp_prim ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile_hosp_prim ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile_hosp_prim ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n

file close textfile_hosp_prim

*Fully adjusted model*
stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base //
estimates store a
quietly stcox i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001





