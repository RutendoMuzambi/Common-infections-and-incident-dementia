*****************************
***Frequency
**modified 03.2021 RM
*****************************

use $cleandata\main_analysis_recurrent_test, clear
cap log close
codebook patid //
stset end, failure(dementia_status=1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 

replace infection_date1= end+365.25 if infection_date1==.
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

stptime if infection_tv==0,  title(person-years) per(1000)

list patid _t0 _t start_fup end infection_number infection_tv infection_date1 dementia_date1 in 1/20


foreach i of num 1/77 {
replace infection_date`i'=d(01jan2030) if infection_date`i'==.
}

stsplit infection_num_tv=infection_date1, at(0) 
replace infection_num_tv=infection_num_tv + 1

foreach i of num 2/77 {
stsplit infection_num_tv`i'=infection_date`i', at(0)
replace infection_num_tv`i'=infection_num_tv`i' + 1
replace infection_num_tv = infection_num_tv + infection_num_tv`i'
}

lab var infection_num_tv "Number of infections"

gen countover1=infection_num_tv-1
replace countover1=0 if countover1==-1

gen first_infection=infection_tv==1 & countover1==0

list patid _t0 _t start_fup end infection_num_tv infection_date1 first_infection in 1/20

log using $log/frequency_cox.log, replace

strate countover1, per(1000) output ($results/frequency, replace)


stptime if infection_tv==0, title(person-years) per(1000)

stptime if first_infection==1, title(person-years) per(1000)

stptime if countover1>=1, title(person-years) per(1000)

*Age adjusted analysis
stcox countover1 i.infection_tv, base //


**Minimally adjusted model 
stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab
 
stcox countover1 i.infection_tv i.sex i.imd_cat i.calperiod, base 

***Fully adjusted model
stcox countover1 i.infection_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base //


log close

stcox i.infection_tv countover1, base //
stcox i.infection_tv countover1 i.sex i.imd_cat i.calperiod, base 
stcox i.infection_tv countover1 i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base //


 
/*
stcox countover1 i.infection_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base //
estimates store a
stcox i.infection_tv i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001

stcox countover1 i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base //
estimates store a
stcox i.sex  i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat  i.stroke_cat, base
estimates store b
lrtest a b // P < 0.0001



