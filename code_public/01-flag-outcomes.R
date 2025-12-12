#' ---
#' title: "01 - Define Asthma and Wheeze Outcomes"
#' author: "Annelise Blomberg"
#' date: "`r format(Sys.Date())`"
#' output:
#'      html_document:
#'           self-contained: true
#'           theme: flatly 
#'           number_sections: true
#'           toc: true
#'           toc_depth: 2
#'           toc_float: 
#'                collapsed: true
#' ---
#' 
# Introduction ------------------------------------------------------------
#' # Introduction 
#' 
#' This script constructs register-based outcome algorithms for asthma and wheeze using 
#' inpatient/outpatient visits (PAR) and dispensed medications (LMED). 
#' The outcome definitions follow published algorithms as described in the manuscript. 
#' The script identifies subjects with each outcome, the date of the first event (incidence) and age at incidence. 
#' It also creates subject-level asthma flags used to define parental asthma status. 
#' 
#' ## Inputs 
#' 
#' - Raw healthcare registers: 
#'      - `UT_R_LMED_59731_2023.sas7bdat` (dispensed medications)
#'      - `UT_R_PAR_OV_59731_2023.sas7bdat` (outpatient specialist care)
#'      - `UT_R_PAR_SV_59731_2023.sas7bdat` (inpatient specialist care)
#' - Registry cohort with baseline exposure and covariate data: 
#'      - `child_cohort_exp_time_const_covs_and_baseline.rds`
#'      
#' ## Outputs 
#' 
#' - `outcomes_comb_list.rds`: tibble with one row per outcome algorithm, containing a list-column of subject-level outcome date (ID, event, date, age).  
#' - `outcomes_flag_asth.rds`: subject-level asthma flags used to identify parental asthma status. 
#' 
#' 
# Setup -------------------------------------------------------------------
#' # Setup 

library(here)
library(tidyverse)
library(haven)

output_loc <- here("outputs_public")

sessionInfo()

# Load Data ---------------------------------------------------------------
#' # Load Data 
#' 
#' Start by loading the lmed and PAR files. 
#' 
#' * UT_R_LMED_59731_2023.sas7bdat
#' * UT_R_PAR_OV_59731_2023.sas7bdat
#' * UT_R_PAR_SV_59731_2023.sas7bdat
#' 
#' SV: inpatient register, complete national coverage starting in 1987. 
#' ov: outpatient register, began in 2001. 

dlmed  <- read_sas(here("data", "raw_data", "ut_r_lmed_59731_2023.sas7bdat")) 
dparov <- read_sas(here("data", "raw_data", "ut_r_par_ov_59731_2023.sas7bdat"))
dparsv <- read_sas(here("data", "raw_data", "ut_r_par_sv_59731_2023.sas7bdat"))

dlmed <- janitor::clean_names(dlmed)
dparov <- janitor::clean_names(dparov) %>% mutate(flag_ov = TRUE)
dparsv <- janitor::clean_names(dparsv) %>% mutate(flag_ov = FALSE)

#' Bind rows and drop unused columns 
dpar <- bind_rows(dparov, dparsv) %>% 
  rename(lopnr = lop_nr, date_event = indatum) %>% 
  select(-"sjukhus", -"mvo", -"pvard", -"drg", -"mdc", -"ktyp", -"dia_ant", -"verks_akut", -"slutrapporterad")

dlmed <- dlmed %>% 
  rename(lopnr = lop_nr, date_event = edatum) 

#' Analysis starts in 2006-01-01: first year with PAR and LMED data

start_date <- as.Date("2006-01-01")
dpar  <- dpar %>%  filter(date_event >= start_date)
dlmed <- dlmed %>% filter(date_event >= start_date)


## Create diagnosis file in long format ------------------------------------
#' ## Create diagnosis file in long format
#' 
#' We use all diagnosis codes for all analyses.

#+ pivot-data-long 
dpar_long <- dpar %>% 
  select(c("lopnr", "alder", "date_event", "hdia", "dia1", "dia2", "dia3", "dia4", 
           "dia5", "dia6", "dia7", "dia8", "dia9", "dia10", "dia11", "dia12", 
           "dia13", "dia14", "dia15", "dia16", "dia17", "dia18", "dia19", 
           "dia20", "dia21", "dia22", "dia23", "dia24", "dia25", "dia26", 
           "dia27", "dia28", "dia29", "dia30")) %>% 
  pivot_longer(cols = c("hdia", "dia1", "dia2", "dia3", "dia4", 
                        "dia5", "dia6", "dia7", "dia8", "dia9", "dia10", "dia11", "dia12", 
                        "dia13", "dia14", "dia15", "dia16", "dia17", "dia18", "dia19", 
                        "dia20", "dia21", "dia22", "dia23", "dia24", "dia25", "dia26", 
                        "dia27", "dia28", "dia29", "dia30"), 
               names_to = "diag_num", values_to = "diag") %>% 
  filter(diag != "")



## Load subject-specific data for age calculations ----------------------------------------------
#' ## Load subject-specific data for age calculations
#' 
#' Load subject birthdates so we can calculate ages for healthcare and prescriptions. 

#+ load-birthdates
dpop <- read_rds(here("data", "raw_data", "child_cohort_exp_time_const_covs_and_baseline.rds")) %>% 
  select(lopnr, fodarman)

dpar_long <- dpar_long %>% 
  left_join(dpop, by = join_by(lopnr)) %>% 
  mutate(alder_dpop = time_length(interval(fodarman, date_event), "years")) %>% 
  select(-fodarman)
  
dlmed <- dlmed %>% 
  left_join(dpop, by = join_by(lopnr)) %>% 
  mutate(alder_dpop = time_length(interval(fodarman, date_event), "years")) %>% 
  select(-fodarman)

# Create ICD10 outcome lists ----------------------------------------------
#' # Create ICD10 outcome lists 
#' 
#' 
#' Create table with outcome-specific codes. These ICD codes are specific to our received PAR datasets (all codes that meet the outcome inclusion criteria). 
icd_codes <- tribble(~ algorithm, ~ icd_codes, 
                       "asth_ort", c("J45", "J45-P", "J450", "J4501A", 
                                         "J450A", "J450B", "J450CV", "J450W", "J450Z", "J451", "J451A", 
                                         "J451B", "J451W", "J458", "J4581B", "J4582A", "J4582B", "J4582C", 
                                         "J4583C", "J458BV", "J459", "J459I", "J459O", "J45P"),
                       "wheeze", c("J45", "J45-P", "J450", "J4501A", 
                                       "J450A", "J450B", "J450CV", "J450W", "J450Z", "J451", "J451A", 
                                       "J451B", "J451W", "J458", "J4581B", "J4582A", "J4582B", "J4582C", 
                                       "J4583C", "J458BV", "J459", "J459I", "J459O", "J45P",
                                       "J20-J22"))

#' Table with codes and start dates: 
dlist <- expand_grid(algorithm = icd_codes$algorithm) %>% 
  left_join(icd_codes, by = join_by(algorithm)) 

#' Helper function: keep first ICD-10 event per subject for specified diagnosis_codes. 
CreateICD10 <- function(data, diagnosis_codes, start_date = as.Date("2006-01-01")) { 
  d2 <- data %>% 
    filter(date_event >= start_date & diag %in% diagnosis_codes) %>% 
    arrange(lopnr, date_event) %>% 
    group_by(lopnr) %>% 
    slice(1) %>% 
    ungroup()
  
  return(d2)
}

#' Add data to table. 
#' 
#' The "subj_list" has the first date of ICD10 diagnosis after 2006-01-01 for all 
#' subjects with at least one of the diagnosis codes in each algorithm.

dlist <- dlist %>% 
  mutate(data = "icd", 
         subj_list = map(icd_codes, ~ CreateICD10(data = dpar_long, diagnosis_codes = .x)))



#' Filter "wheeze" outcomes to children less than three years old (e.g., they have 3 lived years; also 36 months).
dlist[dlist$algorithm == "wheeze", "subj_list"][[1]][[1]] <- 
  filter(dlist[dlist$algorithm == "wheeze", "subj_list"][[1]][[1]], alder_dpop <= 3)


#' Finally, manually add a row for asthma 3+ where the ICD diagnosis is limited to subjects older than 36 months. 
asth_ort3_icd <- dpar_long %>% 
  filter(date_event >= "2006-01-01" &
           diag %in% c("J45", "J45-P", "J450", "J4501A", 
                       "J450A", "J450B", "J450CV", "J450W", "J450Z", "J451", "J451A", 
                       "J451B", "J451W", "J458", "J4581B", "J4582A", "J4582B", "J4582C", 
                       "J4583C", "J458BV", "J459", "J459I", "J459O", "J45P") & 
           alder_dpop > 3) %>% # this works because at age 3.1, you are greater than 36 months  
  arrange(lopnr, date_event) %>% 
  group_by(lopnr) %>% 
  slice(1) %>% 
  ungroup()

dlist_ort3 <- tibble(algorithm = "asth_ort3",
                     data = "icd", 
                     subj_list = c(list(asth_ort3_icd)))

dlist_icd <- bind_rows(dlist, dlist_ort3) %>% 
  select(-icd_codes)


# Add LMED data -------------------------------------------------------------------------
#' # Add LMED data 
#' 
#' Next, we add LMED data for each appropriate outcome algorithm. 
#' Each algorithm has different inclusion / exclusion criteria so we do them separately. 
#' 
#' For these algorithms, we want a start date in 2006 and later. 
#' We start with a tibble framework that parallels what we used in the PAR data. 

dlist_lmed <- tibble(algorithm = dlist_icd$algorithm, 
                     data = "lmed", 
                     subj_list = vector(mode = "list", length = 3))

## Asthma --------------------------------------------------------
#' ## Asthma
#' 
#' We identify LMED asthma outcomes for all subjects (asth_ort) and all subjects over age 3 (asth_ort3). 
#' 
# Function `Apply_Ortqvist_Asthma`: Apply the two algorithm criteria and return the first qualifying prescription date per subject. 
Apply_Ortqvist_Asthma <- function(data_in) {
  codes_asthma_ortqvist_c1 <- c("R03BA", "R03DC03", "R03AK", "R03AK06", "R03AK07")
  codes_asthma_ortqvist_c2 <- c("R03BA", "R03DC03", "R03AK", "R03AK06", "R03AK07", 
                                "R03AC", "R03AC02", "R03AC03", "R03AC12", "R03AC13")
  
  # Criteria 1
  d1 <- data_in %>% 
    filter(atc %in% codes_asthma_ortqvist_c1) %>% 
    group_by(lopnr) %>% 
    arrange(lopnr, date_event) %>% 
    mutate(time_since_last_prescription = as.numeric(difftime(date_event, lag(date_event), units = "days")),
           flag_4yr_less_than_2wks = alder_dpop <= 4.5 & time_since_last_prescription <= 14)
  
  #' Filter by age: 
  d2 <- d1 %>% 
    filter(flag_4yr_less_than_2wks != TRUE) %>% 
    group_by(lopnr) %>% 
    count() %>% 
    filter(n > 1) %>% 
    left_join(d1, by = join_by(lopnr))
  
  d3 <- d2 %>% 
    arrange(lopnr, date_event) %>% 
    group_by(lopnr) %>% 
    slice(1)
  
  dc1 <- d3 %>% 
    select(-n, -time_since_last_prescription, -flag_4yr_less_than_2wks)
  
  # Criteria 2
  d1 <- data_in %>% 
    filter(atc %in% codes_asthma_ortqvist_c2) %>% 
    group_by(lopnr) %>% 
    arrange(lopnr, date_event) %>% 
    # Exclude two prescriptions on the same day from counting as two separate prescriptions: 
    mutate(time_to_next_prescription =  as.numeric(difftime(lead(date_event, n = 1), date_event, units = "days"))) %>% 
    filter(time_to_next_prescription > 0) %>% 
    # Count how long until third prescription 
    mutate(time_to_third_prescription = as.numeric(difftime(lead(date_event, n = 2), date_event, units = "days"))) %>% 
    filter(time_to_third_prescription <= 365)
  
  d2 <- d1 %>% 
    arrange(lopnr, date_event) %>% 
    group_by(lopnr) %>% 
    slice(1)
  
  dc2 <- d2 %>% 
    select(-time_to_next_prescription, -time_to_third_prescription)
  
  dfinal <- bind_rows(dc1, dc2) %>% 
    arrange(lopnr, date_event) %>% 
    group_by(lopnr) %>% 
    slice(1) %>% 
    ungroup()
  
  return(dfinal)
}

asth_ort_lmed <- Apply_Ortqvist_Asthma(data_in = dlmed)
asth_ort3_lmed <- Apply_Ortqvist_Asthma(data_in = filter(dlmed, alder_dpop > 3))


#' Assign back to primary table: 
dlist_lmed[dlist_lmed$algorithm == "asth_ort",  "subj_list"][[1]] <- list(asth_ort_lmed)
dlist_lmed[dlist_lmed$algorithm == "asth_ort3", "subj_list"][[1]] <- list(asth_ort3_lmed)

## Wheeze ------------------------------------------------------------------
#' ## Wheeze 
#' 

subj_lmed_wheeze <- dlmed %>% 
  filter(alder_dpop <= 3 & 
           atc %in% c("R03BA", "R03DC03", "R03AK", "R03AK06", "R03AK07", 
                      "R03AC", "R03AC02", "R03AC03", "R03AC12", "R03AC13")) %>% 
  arrange(lopnr, date_event) %>% 
  group_by(lopnr) %>% 
  slice(1) %>% 
  ungroup()

#' Assign back to primary table: 
dlist_lmed[dlist_lmed$algorithm == "wheeze",  "subj_list"][[1]] <- list(subj_lmed_wheeze)


# Combine ICD and LMED ----------------------------------------------------
#' # Combine ICD and LMED
#' 
#' For each algorithm, combine the ICD and LMED subject lists and keep the earliest date per subject. 

CombineSubj <- function(data_icd, data_lmed) {
  # Returns one row per subject
  dcomb <- bind_rows(data_icd, data_lmed) %>% 
    arrange(lopnr, date_event) %>% 
    group_by(lopnr) %>% 
    slice(1) %>% 
    ungroup() 

  return(dcomb)
}


dlist_comb <- bind_rows(dlist_icd, dlist_lmed) %>% 
  pivot_wider(names_from = "data", values_from = "subj_list") %>% 
  mutate(data = "icd_lmed",
         subj_list = map2(icd, lmed, ~ CombineSubj(.x, .y))) %>% 
  select(-icd, -lmed)


## Finalize the Ortqvist age 3 category ------------------------------------
#' ## Finalize the Ortqvist age 3 category
#' 
#' The Ortqvist3 list includes all subjects who meet either a LMED or ICD10 criteria at age 36 months or later.  
#' 
#' But incidence is the first date of onset, which may have occured before age 36 months. 
#' So we take the incidence dates from the original Ortqvist list - but limit it to those in Ort3. 

ort3_subj <- dlist_comb %>% 
  filter(algorithm == "asth_ort3") %>% 
  .$subj_list %>% 
  .[[1]] %>% 
  select(lopnr)

ort3_dat <- dlist_comb %>% 
  filter(algorithm == "asth_ort") %>% 
  .$subj_list %>% 
  .[[1]] %>% 
  semi_join(ort3_subj, by = join_by(lopnr))

#' Merge version with correct events back into full list
dlist_comb[dlist_comb$algorithm == "asth_ort3", ]$subj_list[[1]] <- ort3_dat


# Combine all outcomes and save -------------------------------------------
#' # Combine all outcomes and save
dlist_final <- bind_rows(dlist_icd, dlist_lmed, dlist_comb) %>% 
  mutate(outcome_algo = paste0(algorithm, "_", data)) %>% 
  rename(data_source = data, 
         data = subj_list) %>% 
  filter(data_source == "icd_lmed") %>% 
  select(-data_source)

write_rds(dlist_final, file.path(output_loc, "outcomes_comb_list.rds"))


# Create Parental Asthma Data ---------------------------------------------
#' # Create Parental Asthma Data
#' 
#' Keep asthma outcomes and collapse to one row per subject with a logical flag. 

outcomes_flag_asth <- dlist_final %>% 
  filter(algorithm %in% c("asth_ort")) %>% 
  unnest(cols = c("data")) %>% 
  mutate(flag_event = TRUE)

outcomes_flag_asth <- outcomes_flag_asth %>% 
  pivot_wider(id_cols = lopnr, 
              names_from = outcome_algo, 
              values_from = flag_event, 
              names_prefix = "flag_") %>% 
  mutate(across(c("flag_asth_ort_icd_lmed"), ~ ifelse(is.na(.), FALSE, .)))

write_rds(outcomes_flag_asth, file.path(output_loc, "outcomes_flag_asth.rds"))
