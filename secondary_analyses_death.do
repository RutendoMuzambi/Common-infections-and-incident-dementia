
/*=========================================================================
DO FILE NAME:			Main_analyses_Timing.do

AUTHOR:					Rutendo Muzambi

DATE VERSION CREATED: 	07/2020
						
DESCRIPTION OF FILE:    Time period before death. Proximity of dementia diagnosis to death
*=========================================================================*/


use $cleandata\mainanalysis_firstever, clear

gen ending = min(deathdate,tod,lcd,d(31dec2018))
lab var ending "Follow-up End Date"
format ending %td

gen death_status=1 if !mi(deathdate)
keep if dementia_status==1

stset ending, failure(death_status==1) enter(dementia_date1) id(patid) origin(dementia_date1) scale(365.25) time0(start_fup) 

rename infection_tv Infection

	sts graph, by (Infection) ///
	title("Kaplan - Meier plot") ///	
	xtitle ("Time to death after dementia diagnosis (years)") ///
	risktable(, size(small)) 
	*scale (0.5) 

sts test Infection, logrank	// p < 0.0001

graph export $graphfiles/survival_terminaldecline.pdf, replace


*===============================12 months=======================================

use $cleandata\mainanalysis_firstever, clear

gen ending = min(deathdate,tod,lcd,d(31dec2018))
lab var ending "Follow-up End Date"
format ending %td

keep if dementia_status==1

gen death_status=1 if !mi(deathdate)

stset ending, failure(death_status==1) enter(dementia_date1) id(patid) origin(dementia_date1) scale(30.4) time0(start_fup)  

*keep if _t<=12.0

rename infection_tv Infection

	sts graph, by (Infection) tmax(12) ///
	title("Kaplan - Meier plot") ///	
	xtitle ("Time to death after dementia diagnosis (months)") ///
	risktable(, size(small)) 

graph export $graphfiles/survival_terminaldecline12m.pdf, replace


********************************************************************************
**************12 months closer up
********************************************************************************


use $cleandata\mainanalysis_firstever, clear

gen ending = min(deathdate,tod,lcd,d(31dec2018))
lab var ending "Follow-up End Date"
format ending %td

keep if dementia_status==1

gen death_status=1 if !mi(deathdate)

stset ending, failure(death_status==1) enter(dementia_date1) id(patid) origin(dementia_date1) scale(30.4) time0(start_fup)  

keep if _t<=12.0

rename infection_tv Infection

	sts graph, by (Infection) tmax(12) ///
	title("Kaplan - Meier plot") ///	
	xtitle ("Time to death after dementia diagnosis (months)") ///
	risktable(, size(small)) 

graph export $graphfiles/survival_terminaldecline12m_2.pdf, replace


