#' ---
#' title: "07 - MICE imputation"
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
#' This script constructs outcome-specific analysis datasets and runs 
#' multiple imputation using `MICE` to create 20 imputed datasets for asthma, asthma (3+) and wheeze. 
#' These are later used in the MICE sensitivity analyses.  
#'  
#' ## Inputs 
#' 
#' - `child_cohort_exp_time_const_covs_and_baseline.rds`
#' - `outcomes_flag_asth.rds`
#' - `outcomes_comb_list.rds` (main cohort with outcome-specific follow-up time)
#' 
#' ## Outputs
#' - `mice-imp-20-list.rds` (list of 3 `mids` objects)
#' 
#' 
#' *NOTE* that this script can take several hours to run.


# Setup -------------------------------------------------------------------
#' # Setup

#+ load-libraries, message = F
library(survival)
library(here)
library(tidyverse)
library(mice)

output_loc <- here("outputs_public")

sessionInfo()


# Process Cohort Data -----------------------------------------------------
#' # Process Cohort Data

#' Load population dataset and apply filters. 
#+ load-overall-cohort
dpop <- read_rds(here("data", "raw_data", "child_cohort_exp_time_const_covs_and_baseline.rds"))  %>% 
  filter(fodlandgr_eu27_2020 == "Sverige") %>% 
  filter(fodelsear == rtb_year_min) %>% 
  filter(flag_born_8021_blekinge | flag_born_8021_ronneby) %>% 
  filter(fodelsear >= 2006 & fodelsear <= 2013) %>% 
  mutate(flag_parent_born_abroad = utlsvbakg != "native-born with native-born parents")

#' Set levels for exposure: 
dpop <- dpop %>% 
  mutate(ecat_prenat = factor(ecat_prenat, 
                              levels = c("background", "intermediate", "high"), 
                              labels = c("Low", "Intermediate", "High")))

#' Add four-category exposure variable: 
dpop <- dpop %>% 
  mutate(ecat_prenat_4 = factor(ecat_prenat_4cat, 
                                levels = c("background", "intermediate", "high", "very high"), 
                                labels = c("Low", "Intermediate", "High", "Very High")))

#' Create some LISA variables
dpop <- dpop %>% 
  mutate(yseg  = factor(yseg),
         fodlandgr_eu27_2020 = factor(fodlandgr_eu27_2020),
         paritet_cat = factor(paritet), 
         paritet_cat = fct_collapse(paritet_cat, "5+" = c("5", "6", "7", "8", "9", "10", "12")))

#' Set NA values 
dpop$mat_educ_3cat[dpop$mat_educ_3cat == "Missing"] <- NA
dpop$sun_old_cat[dpop$sun_old_cat == "Missing"] <- NA
dpop$rok0[dpop$rok0 == "Not specified"] <- NA
dpop$rok1[dpop$rok1 == "Not specified"] <- NA
dpop$rok2[dpop$rok2 == "Not specified"] <- NA

dpop$rok1_2cat[dpop$rok1_2cat == "Not specified"] <- NA
dpop$yseg[dpop$yseg == ""] <- NA
dpop$yseg[dpop$yseg == "**"] <- NA

dpop$mat_educ_3cat <- fct_drop(dpop$mat_educ_3cat)
dpop$sun_old_cat <- fct_drop(dpop$sun_old_cat)
dpop$rok0        <- fct_drop(dpop$rok0)
dpop$rok1        <- fct_drop(dpop$rok1)
dpop$rok2        <- fct_drop(dpop$rok2)
dpop$rok1_2cat   <- fct_drop(dpop$rok1_2cat)
dpop$yseg        <- fct_drop(dpop$yseg)
dpop$flag_missing_father <- is.na(dpop$lopnr_far)
dpop$utlsvbakg <- fct_drop(dpop$utlsvbakg)

## Order factors -----------------------------------------------------------
#' ## Order factors
dpop <- dpop %>% 
  mutate(rok0 = ordered(rok0, levels = c("Non-smoker", "1-9 cig/day","10+ cig/day")),
         rok1 = ordered(rok1, levels = c("Non-smoker", "1-9 cig/day","10+ cig/day")),
         rok2 = ordered(rok2, levels = c("Non-smoker", "1-9 cig/day","10+ cig/day")),
         sun_old_cat = ordered(sun_old_cat, 
                               levels = c("Pre-secondary education shorter than 9 years", 
                                          "Pre-secondary education equivalent to 9 years", "Upper-secondary education, no more than 2 years", 
                                          "Secondary education, 3 years", "Post-secondary education shorter than 3 years", 
                                          "Post-secondary education 3 years or longer", "Postgraduate education")))

# Add parental asthma status -----------------------------------------
#' ## Add parental asthma status 
outcomes_flag_asth_mor <- read_rds(file.path(output_loc, "outcomes_flag_asth.rds")) %>% 
  rename(lopnr_mor = lopnr,
         flag_asth_ort_mor = flag_asth_ort_icd_lmed)

outcomes_flag_asth_far <- read_rds(file.path(output_loc, "outcomes_flag_asth.rds")) %>% 
  rename(lopnr_far = lopnr,
         flag_asth_ort_far = flag_asth_ort_icd_lmed)

dpop <- dpop %>% 
  left_join(outcomes_flag_asth_mor, by = join_by(lopnr_mor)) %>%
  left_join(outcomes_flag_asth_far, by = join_by(lopnr_far)) %>% 
  mutate(across(c("flag_asth_ort_mor", "flag_asth_ort_far"),  ~ ifelse(is.na(.), FALSE, .))) %>% 
  mutate(flag_asth_ort_far = ifelse(flag_missing_father == TRUE, NA, flag_asth_ort_far)) %>%
  mutate(flag_asth_ort_mf =  flag_asth_ort_mor  | flag_asth_ort_far,
         flag_asth_ort_mf =  ifelse(is.na(flag_asth_ort_mor)   |  is.na(flag_asth_ort_far), NA, flag_asth_ort_mf))

# Select variables --------------------------------------------------------
#' # Select variables 
#' 
#' Start by selecting relevant variables at the most original level. 
dpop_lim <- dpop %>% 
  select(c("lopnr", "lopnr_mor", 
           "rok0", "rok1", "rok2", "paritet", "kon", 
           "utlsvbakg", "sun_old_cat", "yseg", "dispinkfam_comb", 
           "mat_age_delivery", "ecat_prenat_4", 
           "dispink04", "dispinkke04", "flag_asth_ort_mf")) 

#' Add censor dates 
ddates <- dpop %>% 
  select(lopnr, fodarman, doddatum, rtb_year_max)

dpop_lim2 <- dpop_lim %>% 
  left_join(ddates, by = join_by(lopnr))


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

dpop_lim2_12 <- AddCensorDates(dpop_lim2, max_age = 12) 
dpop_lim2_2  <- AddCensorDates(dpop_lim2,  max_age = 2) 


# Add Outcome Data --------------------------------------------------------
#' # Add Outcome Data

dout_list <- read_rds(file.path(output_loc, "outcomes_comb_list.rds"))

dpop_out_asthma <- dout_list[dout_list$algorithm == "asth_ort", ]$data[[1]] %>%  
  right_join(dpop_lim2_12, by = join_by(lopnr))

dpop_out_asthma3 <- dout_list[dout_list$algorithm == "asth_ort3", ]$data[[1]] %>%  
  right_join(dpop_lim2_12, by = join_by(lopnr))

dpop_out_wheeze <- dout_list[dout_list$algorithm == "wheeze", ]$data[[1]] %>%  
  right_join(dpop_lim2_2, by = join_by(lopnr))

dpop_out_list <- list(asthma  = dpop_out_asthma,
                      asthma3 = dpop_out_asthma3,
                      wheeze  = dpop_out_wheeze)

#' Flag outcome events and add dates
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

dpop_out_list <- map(dpop_out_list, FlagEvents)
dpop_out_list <- map(dpop_out_list, ~ mutate(.x, nelson_aalen = mice::nelsonaalen(.x, timevar = "age_final", statusvar = "flag_event")))
dpop_out_list <- map(dpop_out_list, ~ select(.x, c(colnames(dpop_lim), c("flag_event", "nelson_aalen", "age_final"))))


# Impute Missings ---------------------------------------------------------
#' # Impute Missings
#' 
#' 
#' Use the asthma dataset to derive a common predictor matrix used across all three outcomes. 
ini <- dpop_out_list[[1]] %>% 
  mice(maxit = 0, seed = 123)

#' Double check correct imputation methods depending on variable type: no changes needed. 
ini$method

#' Change predictor matrix for variables that we are not imputing. 
pred2 <- ini$predictorMatrix
pred2[, "lopnr"] <- 0
pred2[, "lopnr_mor"] <- 0
pred2[, "age_final"] <- 0

#' Impute with continuous parity
#+ impute-data, message = FALSE, warning = FALSE
imp_list <- map(dpop_out_list, 
                ~ mice(data = .x, 
                       maxit = 20, 
                       m = 20, 
                       pred = pred2, 
                       seed = 123,
                       print = TRUE))


write_rds(imp_list, file.path(output_loc, "mice-imp-20-list.rds"))

#' Plot imputations (DQ check only)
# names(imp_list)
# plot(imp_list[[1]])
# plot(imp_list[[2]])
# plot(imp_list[[3]])
