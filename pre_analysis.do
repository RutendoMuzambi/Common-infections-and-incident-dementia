/*=========================================================================
DO FILE NAME:			pre_analysis.do

AUTHOR:					Rutendo Muzambi

VERSION:				v2

DATE VERSION CREATED: 	03/2020
						
*=========================================================================*/

*****************************************************************************

clear all
capture log close

use $cleandata\maincohort_final, clear
count if dementia_status==1 //
count if infection_cat==1 & (dementia_date1-infection_date1)<=91.3125 //
**

use $cleandata\maincohort_final, clear
count if dementia_status==1 //
codebook patid //
replace end=infection_date1 if infection_cat==1 //
save $cleandata\end_at_infection, replace


use $cleandata\maincohort_final, clear
gen postinfection=infection_date1+91.3125
format postinfection %td
replace start_fup=postinfection if infection_cat==1 
save $cleandata\3mon_postinf, replace

clear all
use $cleandata\3mon_postinf, clear
append using $cleandata\end_at_infection
drop if infection_date1==start_fup //
drop if start_fup>end //
bysort patid infection_cat: drop if _n==1 & infection_cat==0 //

foreach var of varlist dementia_status alzheimers vascular_dem unspecified_dem cause_identified_dem {
replace `var'=. if end<dementia_date1 & !mi(dementia_date1)
replace `var'=. if infection_date1==dementia_date1 & !mi(dementia_date1) & !mi(infection_date1)
replace `var'=0 if `var'==.
}
codebook patid //
save $cleandata\main_analysis1, replace


*****************************************************************************

use $cleandata\main_analysis1, clear
codebook patid //
stset end, failure(dementia_status=1) enter(start_fup) id(patid) origin(dob) scale(365.25) time0(start_fup) 


sort patid start_fup
list patid _t0 _t start_fup end _st _d infection_cat infection_date1 in 1/20
sort patid _t0 _t
by patid: gen gap=1 if _t0!=_t[_n-1] & _n>1
count if gap==1 & _st==1 // 
stdescribe // 

summarize _st 
stvary
count if (dementia_date1-infection_date1)<=91.3125 & dementia_date1>end & _d==1
count if infection_date1==dementia_date1 & !mi(dementia_date1) & !mi(infection_date1) & _d==1
save $cleandata\mainanalysis_cox, replace

*****************************************************************************

replace infection_date1= end+365.25 if infection_date1==. //
stsplit infection_tv=infection_date1, at(0)
replace infection_tv=infection_tv+1

save $cleandata\mainanalysis_firstever, replace

use $cleandata\mainanalysis_firstever, replace

list patid _t0 _t start_fup end _st infection_date1 infection_tv in 1/20
by infection_tv, sort:stptime, title(person-years) per(1000)

sts graph, by(infection_tv) 

stsplit ageband, at (65(5)90)

asdoc by infection_tv, sort : stptime, by(ageband) per(1000) replace dec(2)
strate infection_tv ageband, per(1000) output ($results/ageband, replace)
