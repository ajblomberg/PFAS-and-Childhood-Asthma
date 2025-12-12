#' ---
#' title: "02- Create underlying cohort"
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
#' This script constructs the main study cohort by: 
#' 
#' - Applying geographic and birth year inclusion criteria (Blekinge; births 2006-2013). 
#' - Linking the baseline cohort to the identified outcomes. 
#' - Defining follow-up time, censoring, and event indicators. 
#' - Creating covariates and quantile-based variables used in the models. 
#' 
#' Follow-up time is defined from birth until either an outcome incidence, end of follow-up, outcome-specific maximum age, 
#' death or emigration from Sweden. 
#'  
#' ## Inputs 
#' 
#' - `child_cohort_exp_time_const_covs_and_baseline.rds`
#' - `outcomes_comb_list.rds`
#' - `outcomes_flag_asth.rds`
#'      
#' ## Outputs 
#' 
#' - `cohort-outcome-data.rds`
#'      - Tibble with row per outcome algorithm
#'      - Includes a list-column with child-level cohort information for each algorithm: exposure, covariates, follow-up time, 
#'      and outcome variables.  
#' 
#' 
# Setup -------------------------------------------------------------------
#' # Setup 

library(here)
library(tidyverse)
library(lubridate)

output_loc <- here("outputs_public")

sessionInfo()

# Load and Prepare Cohort --------------------------------------------------
#' # Load and Prepare Cohort
#' 
#' We load the cohort, apply filters, and format necessary covariate and exposure variables. 

#+ load-overall-cohort
dpop <- read_rds(here("data", "raw_data", "child_cohort_exp_time_const_covs_and_baseline.rds")) 

#' Must be born in Blekinge
dpop <- dpop %>% filter(flag_born_8021_blekinge | flag_born_8021_ronneby)

#' Set NA values 
dpop$mat_educ_3cat[dpop$mat_educ_3cat == "Missing"] <- NA
dpop$rok1[dpop$rok1 == "Not specified"] <- NA
dpop$rok1_2cat[dpop$rok1_2cat == "Not specified"] <- NA

dpop$mat_educ_3cat <- fct_drop(dpop$mat_educ_3cat)
dpop$rok1 <- fct_drop(dpop$rok1)
dpop$rok1_2cat <- fct_drop(dpop$rok1_2cat)

#' Add four-category exposure variable: 
dpop <- dpop %>% 
  mutate(ecat_prenat_4 = factor(ecat_prenat_4cat, 
                                levels = c("background", "intermediate", "high", "very high"), 
                                labels = c("Low", "Intermediate", "High", "Very High")))

#' Add variable `utlsvbakg`
dpop <- dpop %>% 
  mutate(flag_parent_born_abroad = utlsvbakg != "native-born with native-born parents", 
         flag_both_parents_born_abroad = utlsvbakg == "native-born with two foreign-born parents")


#' Flag missing fathers
dpop <- dpop %>% mutate(flag_missing_father = is.na(lopnr_far))

# Merge parental asthma status  -----------------------------------------------
#' # Merge parental asthma status 
outcomes_flag_asth_mor <- read_rds(file.path(output_loc, "outcomes_flag_asth.rds")) %>% 
  rename(lopnr_mor = lopnr,
         flag_asth_ort_mor = flag_asth_ort_icd_lmed)

outcomes_flag_asth_far <- read_rds(file.path(output_loc, "outcomes_flag_asth.rds")) %>% 
  rename(lopnr_far = lopnr,
         flag_asth_ort_far = flag_asth_ort_icd_lmed)

dpop <- dpop %>% 
  left_join(outcomes_flag_asth_mor, by = join_by(lopnr_mor)) %>%
  left_join(outcomes_flag_asth_far, by = join_by(lopnr_far)) %>% 
  mutate(across(c("flag_asth_ort_mor", "flag_asth_ort_far"),  ~ replace_na(.x, FALSE)),
         flag_asth_ort_far = ifelse(flag_missing_father, NA, flag_asth_ort_far),
         flag_asth_ort_mf =  flag_asth_ort_mor  | flag_asth_ort_far,
         flag_asth_ort_mf =  if_else(
           is.na(flag_asth_ort_mor) |  is.na(flag_asth_ort_far), 
           NA, 
           flag_asth_ort_mf))


# Flag complete cases -----------------------------------------------------
#' # Flag complete cases 
covariate_variables <-  c("rok1_2cat", "paritet_2cat", "kon", "mat_educ_3cat", "flag_parent_born_abroad", 
                          "dispinkfam_comb", "mat_age_delivery", "utlsvbakg", "flag_asth_ort_mf")

dpop <- mutate(dpop, flag_missing_covariates = if_any(all_of(covariate_variables), ~ is.na(.)))

# Add Censoring Data ------------------------------------------------------
#' # Add Censoring Data
#' 
#' There are four reasons why a participant could be censored: 
#' 
#' * reach maximum age  
#' * die 
#' * move abroad
#' * study ended 
#' 
#' For moving abroad: "rtb_year_max" indicates they lived in Sweden on Dec. 31 on that reported year. 
#' Then we know that they did not live in Sweden on December 31 the following year.
#' Assume live in Sweden until reported otherwise on December 31. 
#' 
AddCensorDates <- function(data, max_age) {
  end_of_study <- ymd("2022-12-31")
  d1 <- data %>%
    mutate(date_agemax = fodarman %m+% years(max_age + 1),
           date_died = doddatum,
           date_no_rtb = ymd(paste0(rtb_year_max + 1, "-12-31"))) %>%
    mutate(
      date_censor = pmin(date_agemax, date_died, date_no_rtb, end_of_study, na.rm = T),
      censor_category = case_when(
        date_censor == end_of_study ~ "study_end",
        date_censor == date_agemax  ~ "age_max",
        date_censor == date_no_rtb  ~ "no_rtb",
        date_censor == date_died    ~ "died",
        TRUE ~ NA_character_))
  
  return(d1)
}


#' We include maximum age 3 to combine with wheezing.
dlist <- tibble(age_max = c(12, 2)) %>% # when turn 13 and when turn 3 - there is a "+1" in the function above.
  mutate(data = map(age_max, ~ AddCensorDates(data = dpop, max_age = .x)))


# Add outcome data --------------------------------------------------------
#' # Add outcome data 
dout_list <- read_rds(file.path(output_loc, "outcomes_comb_list.rds")) %>% 
  rename(data_out = data) 

#' Merge to existing data: 
dlist2 <- dout_list %>% 
  mutate(age_max = ifelse(algorithm == "wheeze", 2, 12)) %>% 
  left_join(dlist, by = join_by(age_max)) %>% 
  mutate(data = map2(data, data_out, ~ left_join(.x, .y, by = join_by(lopnr))))

#' Given cohort data with outcome and censor dates, define final event date, event indicator, and age at event/censoring. 
FlagEvents <- function(data) { 
  d1 <- data %>% 
    mutate(date_final = ifelse(is.na(date_event), date_censor, 
                               pmin(date_event, date_censor))) %>%
    mutate(date_final = as.Date(date_final)) %>% 
    mutate(flag_event = ifelse(is.na(date_event), FALSE, 
                               ifelse(date_event == date_final, TRUE, FALSE))) %>% 
    mutate(age_final = as.period(interval(fodarman, date_final), "years"),
           age_final = as.numeric(age_final, "years"))
  
  return(d1)
  }


dlist2 <- dlist2 %>% 
  mutate(data = map(data, FlagEvents))


# Create cohort years --------------------------------------------------------------
#' # Create cohort years
dlist3 <- dlist2 %>% 
  mutate(data = map(data, ~ filter(.x, fodelsear >= 2006 & fodelsear <= 2013))) %>% 
  select(outcome_algo, algorithm, age_max, data)

# Add Quartile Data -------------------------------------------------------
#' # Add Quartile Data 

AddQuantiles <- function(data) { 
  d1 <- data %>% 
    mutate(
      mat_age_delivery_quant = cut(
        mat_age_delivery, 
        breaks = quantile(mat_age_delivery, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
        include_lowest = TRUE),
      dispinkfam_quant = cut(
        dispinkfam_comb, 
        breaks = quantile(dispinkfam_comb, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
        include_lowest = TRUE))
  
  return(d1)
}

dlist4 <- dlist3 %>% 
  mutate(
    data = map(data, ~ AddQuantiles(.x)),
    algorithm_label = factor(
      algorithm, 
      levels = c("wheeze", "asth_ort", "asth_ort3"), 
      labels = c("Wheeze", "Asthma", "Asthma (3+)")))

# Save Data ---------------------------------------------------------------
#' # Save Data 

write_rds(dlist4, file.path(output_loc, "cohort-outcome-data.rds"))

