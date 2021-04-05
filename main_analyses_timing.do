/*=========================================================================
DO FILE NAME:			Main_analyses_Timing.do

AUTHOR:					Rutendo Muzambi


DATE VERSION CREATED: 	05/2020
						
DESCRIPTION OF FILE:    Timing of infections on dementia
*=========================================================================*/

*Time following infection diagnoses

*===============Cumulative time======================================================
use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log/timing_cox.log, replace

replace dementia_status=0 if dementia_status==.

list patid _t0 _t start_fup end _st _d infection_cat infection_date1 in 1/20
stsplit inf_time=infection_date1, at(0.25,1,2,3,4,5,6,7,8,9)
list patid _t0 _t start_fup end _st _d infection_date1 inf_time infection_cat in 1/20

replace inf_time=inf_time*100
tab inf_time

label define inf_time_lab 0 "0" 25"3months-1year" 100"1-2 years" 200"2-3 years" 300"3-4 years" 400"4-5 years" 500"5-6 years" 600"6-7 years" 700"7-8 years" 800"8-9 years" 900"9+ years"
label values inf_time inf_time_lab

*Overlapping time periods
gen inf_time0=1 if inf_time==0
gen inf_time1=1 if inf_time==0 | inf_time==25
gen inf_time2=1 if inf_time==0 | inf_time==25 | inf_time==100 
gen inf_time3=1 if inf_time==0 | inf_time==25 | inf_time==100 | inf_time==200 
gen inf_time4=1 if inf_time==0 | inf_time==25 | inf_time==100 | inf_time==200 | inf_time==300 
gen inf_time5=1 if inf_time==0 | inf_time==25 | inf_time==100 | inf_time==200 | inf_time==300 | inf_time==400
gen inf_time6=1 if inf_time==0 | inf_time==25 | inf_time==100 | inf_time==200 | inf_time==300 | inf_time==400 | inf_time==500
gen inf_time7=1 if inf_time==0 | inf_time==25 | inf_time==100 | inf_time==200 | inf_time==300 | inf_time==400 | inf_time==500 | inf_time==600
gen inf_time8=1 if inf_time==0 | inf_time==25 | inf_time==100 | inf_time==200 | inf_time==300 | inf_time==400 | inf_time==500 | inf_time==600 | inf_time==700
gen inf_time9=1 if inf_time==0 | inf_time==25 | inf_time==100 | inf_time==200 | inf_time==300 | inf_time==400 | inf_time==500 | inf_time==600 | inf_time==700 | inf_time==800
gen inf_time10=1 if inf_time==0 | inf_time==25 | inf_time==100 | inf_time==200 | inf_time==300 | inf_time==400 | inf_time==500 | inf_time==600 | inf_time==700 | inf_time==800 | inf_time==900

list patid _t0 _t start_fup end infection_date1 inf_time inf_time3 inf_time4 in 1/20

*** Calendar period
stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

capture file close textfile 
file open textfile using $results/timing_cox.csv, write replace
file write textfile "sep=;" _n
file write textfile "Effect of common infections on dementia, by time since infection" _n _n
file write textfile  ";"   "Number of dementia events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

*incidence rate*
file write textfile "Time since infection" _n

foreach var in inf_time1 inf_time2 inf_time3 inf_time4 inf_time5 inf_time6 inf_time7 inf_time8 inf_time9 inf_time10 {
*strate if `var'==1, per (1000) output ($results/timing, replace) 
stptime if infection_tv==1 & `var'==1,  title(person-years) per(1000)
**infection
local persontime_inf=`r(ptime)'
local dementiaevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

stptime if infection_tv==0 & `var'==1,  title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local dementiaevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

*age adjusted analysis
stcox i.infection_tv  if `var'==1, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if `var'==1, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if `var'==1, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n
file write textfile "`var'" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n
}

file close textfile
log close



clear all
import excel $results/timing_data, sheet("cumulative") cellrange(A1:D11) firstrow
label var Timesinceinfection "Time since infection"

gen logHR = log(HR)
gen loglci = log(lci)
gen loguci = log(uci)


metan logHR loglci loguci, eform effect(Fully adjusted HR) lcols(Timesinceinfection) astext(50) xsize(26) ysize(15) nowt nobox nooverall xlabel(1,2.0) graphregion(color(white))


*===============Discrete time======================================================
use $cleandata\mainanalysis_firstever, clear
cap log close
log using $log/timing_dis.log, replace

replace dementia_status=0 if dementia_status==.

list patid _t0 _t start_fup end _st _d infection_cat infection_date1 in 1/20
stsplit inf_time=infection_date1, at(0.25,1,2,3,4,5,6,7,8,9)
list patid _t0 _t start_fup end _st _d infection_date1 inf_time infection_cat in 1/20

replace inf_time=inf_time*100
tab inf_time

label define inf_time_lab 0 "0" 25"3months-1year" 100"1-2 years" 200"2-3 years" 300"3-4 years" 400"4-5 years" 500"5-6 years" 600"6-7 years" 700"7-8 years" 800"8-9 years" 900"9+ years"
label values inf_time inf_time_lab

*Discrete time periods
gen inf_time1=1 if inf_time==0 | inf_time==25 
gen inf_time2=1 if inf_time==0 | inf_time==100 
gen inf_time3=1 if inf_time==0 | inf_time==200 
gen inf_time4=1 if inf_time==0 | inf_time==300 
gen inf_time5=1 if inf_time==0 | inf_time==400
gen inf_time6=1 if inf_time==0 | inf_time==500
gen inf_time7=1 if inf_time==0 | inf_time==600
gen inf_time8=1 if inf_time==0 | inf_time==700
gen inf_time9=1 if inf_time==0 | inf_time==800
gen inf_time10=1 if inf_time==0 | inf_time==900

list patid _t0 _t start_fup end infection_date1 infection_tv inf_time inf_time1 inf_time2 in 1/20

*** Calendar period
stsplit calperiod, after(mdy(1,1,2004)) at(0,5,10)
label define calperiod_lab 0"2004-2008" 5"2009-2013" 10"2014-2018"
label values calperiod calperiod_lab

capture file close textfile 
file open textfile using $results/timing_dis.csv, write replace
file write textfile "sep=;" _n
file write textfile "Effect of common infections on dementia, by time since infection" _n _n
file write textfile  ";"   "Number of dementia events" ";" "Person years at risk" ";" "Crude incidence rate (per 1000 person-years)"    
file write textfile  ";"   "Age-adjusted HR* (95% CI)" ";" "Minimally-adjusted HR** (95% CI)" ";" "Fully adjusted HR*** (95% CI)" _n 

*incidence rate*
file write textfile "Time since infection" _n

foreach var in inf_time1 inf_time2 inf_time3 inf_time4 inf_time5 inf_time6 inf_time7 inf_time8 inf_time9 inf_time10 {
*strate if `var'==1, per (1000) output ($results/timing, replace) 
stptime if infection_tv==1 & `var'==1,  title(person-years) per(1000)
**infection
local persontime_inf=`r(ptime)'
local dementiaevents_inf=`r(failures)'
local cruderate_inf=`r(rate)'
local lci_rate_inf=`r(lb)'
local uci_rate_inf=`r(ub)'

stptime if infection_tv==0 & `var'==1,  title(person-years) per(1000)
return list
*no infection
local persontime_no=`r(ptime)'
local dementiaevents_no=`r(failures)'
local cruderate_no=`r(rate)'
local lci_rate_no=`r(lb)'
local uci_rate_no=`r(ub)'

*age adjusted analysis
stcox i.infection_tv  if `var'==1, base 
local hr_age=exp(_b[1.infection_tv])
local lci_age=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_age=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

*minimally adjusted 
stcox i.infection_tv i.sex i.imd_cat i.calperiod  if `var'==1, base
local hr_min=exp(_b[1.infection_tv])
local lci_min=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_min=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])	

stcox i.infection_tv i.sex i.imd_cat i.calperiod i.ethnicity_cat i.smokestatus i.alcohol_cat i.anx_dep_cat i.smi_cat i.ibd_cat i.multiple_sclerosis_cat i.RA_cat i.psoriasis_cat i.asthma_cat  i.ckd_cat i.cld_cat i.copd_cat i.diabetes_cat i.hypertension_cat i.obs_sleep_cat i.tbi_cat i.benzodiazepine_cat i.ppi_cat i.steroids_cat i.polypharmacy_cat i.heart_failure_cat i.stroke_cat if `var'==1, base //
local hr_ful=exp(_b[1.infection_tv])
local lci_ful=exp(_b[1.infection_tv]-1.96*_se[1.infection_tv]) 
local uci_ful=exp(_b[1.infection_tv]+1.96*_se[1.infection_tv])

file write textfile "No infection" ";" (`dementiaevents_no') ";" %3.1f (`persontime_no') ";"
file write textfile %3.2f (`cruderate_no') "(" %3.2f (`lci_rate_no') "-" %3.2f (`uci_rate_no') ")" _n
file write textfile "`var'" ";" (`dementiaevents_inf') ";" %3.1f (`persontime_inf')  ";"
file write textfile %3.2f (`cruderate_inf') "(" %3.2f (`lci_rate_inf') "-" %3.2f (`uci_rate_inf') ")"
file write  textfile ";"  %3.2f (`hr_age') " (" %3.2f (`lci_age')  "-" %3.2f (`uci_age') ")"  
file write  textfile ";"  %3.2f (`hr_min') " (" %3.2f (`lci_min')  "-" %3.2f (`uci_min') ")"  
file write  textfile ";"  %3.2f (`hr_ful') " (" %3.2f (`lci_ful')  "-" %3.2f (`uci_ful') ")" _n
}

file close textfile
log close


clear all
import excel $results/timing_data, sheet("Discrete") cellrange(A1:D11) firstrow
label var Timesinceinfection "Time since infection"


gen logHR = log(HR)
gen loglci = log(lci)
gen loguci = log(uci)


metan logHR loglci loguci, eform effect(Fully adjusted HR) lcols(Timesinceinfection) astext(50) xsize(26) ysize(15) nowt nobox nooverall xlabel(1,2.0) graphregion(color(white))
