# Prenatal PFAS and Childhood Asthma – Analysis Code

This repository contains the analysis code for the manuscript:

> **Prenatal exposure to perfluoroalkyl substance (PFAS) and incidence of asthma and wheeze in childhood: A cohort study in Ronneby, Sweden**  
> Annelise J. Blomberg, Christel Nielsen, Beata Borgström Bolmsjö, Marie-Abèle Bind, Linda Hartman, Anna Saxne Jöud

## 1. Overview
The underlying individual-level data used in these analyses are not shared in this repository because they contain personally identifiable health information and are subject to ethical and legal restrictions. Intermediate analysis datasets are also not included in the repository for the same reasons. Importantly, all intermediate datasets are generated programmatically from the input data using the scripts included in this repository. No manual data cleaning or post-hoc modification were performed outside the shared scripts.

This repository therefore includes: 

- All analysis scripts used to clean data, construct variables, and fit models 
- Scripts that generate intermediate datasets
- Scripts used to generate analysis results and manuscript figures and tables 
- Exported figures and tables, which match the final manuscript
- Log of the complete final analysis run. 

The repository does *not* include: 

- De-identified individual-level data
- Intermediate analysis datasets 

Because the underlying data are not shared, users cannot directly reproduce the numerical results and manuscript figures and tables. Instead, this repository is intended to: 

- Make the logical structure of the analyses transparent 
- Document all data-processing steps and modeling decisions 
- Allow others to review or adapt the code to similar datasets.  

## 2. Analysis data (RDS files)
The data objects used and created in this analysis cannot be shared and are not included in the public repository. Below is a short description of each object. 

Study analyses are based on the following input data files: 
1. Baseline cohort file that includes exposure estimates, baseline covariates, and time-constant information:  `child_cohort_exp_time_const_covs_and_baseline.rds`  
2. Outcome data from the NPR and PDR registries:
	- `UT_R_LMED_59731_2023.sas7bdat`
	- `UT_R_PAR_OV_59731_2023.sas7bdat`
	- `UT_R_PAR_SV_59731_2023.sas7bdat`
3. Primary health records from Region Blekinge: `blekinge_care_contacts_clean.rds` 
4. Residential history information over time: `residence_flags_annual.rds` 

From these input data files, we construct:
1. **Outcome lists** (incidence of asthma/wheeze events)  
2. **Main analysis cohort** used in all primary models and KM plots  
3. **Validation cohort** for comparing register-based outcomes to primary care data  
4. **Imputation objects** for MICE-based sensitivity analyses
5. **Matched cohorts** for randomization-based inference  

## 3. Script structure
Scripts are numbered in approximate execution order. Each script contains a header describing its purpose, as well as code inputs and outputs. 

`01-flag-outcomes.R`  
  - Uses raw outcome registers (PAR, LMED) and `child_cohort_exp_time_const_covs_and_baseline.rds` to identify age at diagnosis/dispensations.  
  - Creates outcome-specific event flags, date of the event, and the subject's age for all subjects identified as having each event. 
  - Outputs:
    - `outcomes_comb_list.rds`
    - `outcomes_flag_asth.rds`

`02-create-cohort.R`  
  - Reads baseline cohort and outcome lists.  
  - Creates cohort by defining and applying inclusion criteria; identifying censoring events and calculating follow-up time; and adding covariates and quantile variables.  
  - Outputs:
    - `cohort-outcome-data.rds`

`03-cohort-description.R`  
  - Uses `cohort-outcome-data.rds` and `data_exp_validation cohort.rds` 
  - Creates descriptive statistics, cohort tables and summary figures.
  - Outputs: 
	  - `ts2-complete-vs-missing.csv`
	  - `t1-baseline-covariates.csv`
	  - `t2-exposure-verification.csv`
	  - `fs1-exp-verification.tiff`
	  - `fs2-outcome-combos.tiff`

`04-outcome-validation-against-primary.R`  
  - Builds a separate validation cohort (subset with primary care data).  
  - Creates comparison of register-based algorithms to primary care outcomes.  
  - Outputs:
    - `ts3-outcome-validation.csv`

`05-plot-KM.R`  
  - Uses `cohort-outcome-data.rds`.  
  - Creates KM plots by exposure categories.  
  - Outputs:
    - `f1-km-plots-primary.(tiff/eps)`

`06-primary-and-sex-stratified-results.R`  
  - Uses `cohort-outcome-data.rds`.  
  - Fits primary Cox models and sex-stratified models.  
  - Outputs:
    - `t3-primary-model-results.csv`
    - `results-primary-models.rds`
    - `fs4-primary-model-covariate-estimates.tiff`
    - `ts6-hazard-ratio-by-sex.csv`
    - `fs6-sex-stratified-estimates.tiff`

`07-run-MICE.R`  
  - Uses baseline cohort and outcome lists.  
  - Creates multiple imputation data for sensitivity analyses.  
  - Outputs:
    - `mice-imp-20-list.rds`

`08-MICE-models.R`
  - Uses `mice-imp-20-list.rds` 
  - Fits Cox models on the imputed datasets and pools results and creates plots comparing MICE results to primary model results 
  - Outputs:
    - `ts5-primary-vs-mice-model-results.csv`
    - `fs5-compare-sensitivity-models.tiff`

`09-RCM-matching.R`  
  - Uses `cohort-outcome-data.rds` 
  - Creates a matched cohort (Very High vs Low exposure).  
  - Outputs:
    - `dvh_matched.rds`
    - `fs7-love-plot.tiff`
    - `dvh_matched_out_list.rds`

`10-RCM-analysis.R`  
  - Uses the matched cohort (`dvh_matched_out_list.rds`)   
  - Creates matched KM curves and randomization-based p-values.  
  - Outputs:
    - `fs8-fisher-test-statistics-distributions.(tiff/png)`
    - `f2-km-plots-rcm.(tiff/eps)`
    - `fs9-adjusted-cumulative-incidence-matched-sensitivity.tiff`

## 4. Software and environment
R version and package information for each specific script are noted in the included log. 
- R version: 4.5.1 (2025-06-13)
- Platform: x86_64-w64-mingw32/x64
- Key packages: Tidyverse (2.0.0), survival (3.8-3), here (1.0.2), MatchIt (4.7.2), and mice (3.18). 

