 
                     AERMOD Test Cases - v12345                  12/10/2012

This file describes test cases for the AERMOD dispersion model, dated 12345.
Due to the range of changes incorporated in v12345 of the AERMOD dispersion 
model (see AERMOD MCB#8) in addition to the changes incorporated in v12345 
of the AERMET meteorological processor (see AERMET MCB#3), the AERMOD test
cases include several scenarios in order to document test results for AERMOD 
changes only vs. AERMET changes only vs. AERMOD and AERMET changes combined.
In particular, the set of test scenarios provided serves to highlight the 
potential impacts of two key changes that may have widespread effects on model 
results, i.e., the convective mixing height bug fix in AERMET, and the vector 
wind speed adjustment formulation bug fix in AERMOD.

The following list summarizes the six separate scenarios included in the set
of test cases for this update to AERMOD:

                                                             VectorWS in 
Test No.   AERMOD/AERMET Runs    Description                AERMOD_v12345
-------------------------------------------------------------------------
   1       12060_Def_11059_Def   12060 AERMOD/11059 AERMET        NO               
   2       12060_Def_12345_Def   12060 AERMOD/12345 AERMET        NO
   3       12345_Def_11059_Def   12345 AERMOD/11059 AERMET        NO
   4       12345_Def_12345_Def   12345 AERMOD/12345 AERMET        NO
   5       12345_Vect_11059_Def  12345 AERMOD/11059 AERMET       YES
   6       12345_Vect_12345_Def  12345 AERMOD/12345 AERMET       YES
   
The '_Def' qualifier indicates that Default options are used, and '_Vect'
indicates that the new AERMOD VectorWS option was used.  Input and output files
associated with each of the six test cases are available for downloading
separately. 

Based on these six AERMOD/AERMET scenarios, the following comparisons have
been summarized in the 'aermod_comparisons_stats.xlsx' file included in the
aermod_testcase.zip file:

   AERMOD/AERMET            AERMOD/AERMET                  Results
-------------------       -------------------   --------------------------------
12345_Def_12345_Def   vs. 12060_Def_11059_Def   Net AERMOD/AERMET Default changes
12345_Vect_11059_Def  vs. 12060_Def_11059_Def   New AERMOD w/VectorWS vs. Old AERMOD
12345_Def_12345_Def   vs. 12060_Def_12345_Def   VectorWS Bug Fix changes with NewMet
12345_Def_11059_Def   vs. 12060_Def_11059_Def   VectorWS Bug Fix changes with OldMet
12345_Vect_12345_Def  vs. 12060_Def_11059_Def   Convective Mix Height Bug Fix changes
  
In addition to the summary of model comparison stats from these tests provided in
the 'aermod_comparisons_stats.xlsx' file, the complete test results for all of these 
scenarios have been incorporated in the 'AERMET_AERMOD_Test_Comparisons.xlsx' file
that can be used to generate detailed comparisons by averaging period for each of
the sources based on user-selected scenarios, using the yellow-highlighted drop down
list on the 'AERMOD Comparison Results' tab. These detailed comparisons are comparable
to the detailed summaries of test results provided with previous versions of AERMOD.
A summary of the scenarios included in these tests is provided below:

Test Name             Description
------------------    ---------------------------------------------------------
aertest.inp           Simple test case referenced in AERMOD User's Guide
allsrcs.inp           Test case including all standard source types
capped.inp            Capped and horizontal stack tests with BETA option
capped_nostd.inp      Capped and horizontal stack tests using MC w/NOSTD option
testgas.inp           Gas deposition test case with limited screening met data
testgas2.inp          Gas deposition test case with 1 year of Houston met data
testpart.inp          Particle deposition test with limited screening met data
testprt2.inp          Particle deposition test with 1 year of Houston met data
openpits.inp          OPENPIT source test case
lovett.inp            Lovett model evaluation data
flatelev.inp          FLAT & ELEV test case based on Lovett data
mcr.inp               Martin's Creek model evaluation data
multurb.inp           Multiple urban area test case
hrdow.inp             Hour-of-day by day-of-week EMISFACT test case
surfcoal.inp          Surface coal mine evaluation data
olm.inp               OLM test case
olmgrp.inp            OLM test case with OLMGROUP option
pvmrm.inp             PVMRM test case
psdcred.inp           PVMRM test case with PSDCREDIT BETA option
scimtest.inp          SCIM (Sampled Chronological Input Model) test case
testpm25.inp          PM-2.5 test case with 5-years of met data
testpm10_1986.inp     PM-10 test case with MULTYEAR option - Year 1
testpm10_1987.inp     PM-10 test case with MULTYEAR option - Year 2
testpm10_1988.inp     PM-10 test case with MULTYEAR option - Year 3
testpm10_1989.inp     PM-10 test case with MULTYEAR option - Year 4
testpm10_1990.inp     PM-10 test case with MULTYEAR option - Year 5
testpm10.inp          PM-10 test case with single 5-year data file (1986-1990)
