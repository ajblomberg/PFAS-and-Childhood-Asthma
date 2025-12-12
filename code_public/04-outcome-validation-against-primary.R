#' ---
#' title: "04 - Outcome Validation against Primary Care"
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



# Introduction ------------------------------------------------------------
#' # Introduction 
#' 
#' This script validates register-based asthma and wheeze outcome algorithms 
#' against primary care data from Blekinge. It: 
#' 
#' - Constructs a validation cohort of children born in Blekinge county between 2010-2021 who are followed as long as they are resident in the county.  
#' - Harmonizes follow-up time and censoring across data sources
#' - Applies asthma and wheeze algorithms to NPR and LMED data, and also identifies outcomes in primary care
#' - Computes sensitivity, specificity, PPV and NPV for each algorithm. 
#' 
#' ## Inputs 
#' 
#' - `child_cohort_exp_time_const_covs_and_baseline.rds`
#' - `residence_flags_annual.rds`
#' - `outcome_comb_list.rds`
#' - `blekinge_care_contacts_clean.rds`
#' 
#' ## Outputs
#' - `ts3-outcome-validation.csv`: sensitivity and specificity table 
#' - Inline text for inclusion in the manuscript 


# Setup -------------------------------------------------------------------
#' # Setup 

#+ load-libraries, message = F
library(here)
library(tidyverse)
library(gt)

output_loc <- here("outputs_public")

sessionInfo()

# Load Data ---------------------------------------------------------------
#' # Load Data 

#+ placeholder

## Population Data ---------------------------------------------------------
#' ## Population Data 
dpop <- read_rds(here("data", "raw_data", "child_cohort_exp_time_const_covs_and_baseline.rds")) %>% 
  filter(flag_born_8021_blekinge | flag_born_8021_ronneby) %>% 
  filter(fodelsear >= 2010 & fodelsear <= 2021) %>% 
  select(lopnr, fodarman, fodelsear)

#' Limit to subjects while they live in Blekinge
dres <- read_rds(here("data", "raw_data", "residence_flags_annual.rds")) %>% 
  semi_join(dpop)

pop_move <- dres %>% 
  filter(flag_res_xblekinge) %>% 
  group_by(lopnr) %>% 
  arrange(lopnr, year) %>% 
  slice(1) %>% 
  select(lopnr, year_move = year)

#' Merge back to primary population: 
dpop2 <- dpop %>% 
  left_join(pop_move, by = join_by(lopnr)) %>% 
  mutate(
    year_move = as.integer(year_move),
    
    # if there is a move year, censor at Dec. 31 of previous year; 
    # otherwise at 2021-12-31
    
    censor_year = ifelse(is.na(year_move), 2021L, year_move - 1L), 
    
    # build date string and parse 
    censor_date_str = paste0(censor_year, "-12-31"),
    date_censor = ymd(censor_date_str)
  ) %>% 
  select(-censor_year, -censor_date_str)


## NPR Outcome Data ------------------------------------------------------------
#' ## NPR Outcome Data
dout <- read_rds(file.path(output_loc, "outcomes_comb_list.rds")) %>% 
  mutate(data_source = "icd_lmed")

#' Filter to subjects who had the outcome before the maximum age 
dout <- dout %>% 
  mutate(max_age = ifelse(algorithm == "wheeze", 3, 13)) %>% 
  mutate(data = map2(data, max_age, ~ filter(.x, alder_dpop <= .y)))

#' Rename wheeze to note that it is J20-J22.  
dlist_out <- dout %>% 
  mutate(algorithm = ifelse(algorithm == "wheeze", "wheeze_j20_j21_j22", algorithm)) %>% 
  select(-outcome_algo)

## Primary Care Data -------------------------------------------------------
#' ## Primary Care Data 
dprim <- read_rds(here("data", "raw_data", "blekinge_care_contacts_clean.rds")) %>%
  filter(year >= 2010) %>% 
  mutate(indatum = as.Date(indatum))

dprim_long <- dprim %>% 
  select(lopnr, indatum, utdatum, year, diagnos1:diagnos8) %>% 
  pivot_longer(
    cols = starts_with("diagnos"),
    names_to = "diag_num", 
    values_to = "diag"
  ) %>% 
  filter(!is.na(diag), diag != "")

#' Limit primary care data to population in the primary registry cohort and add ages   
dprim_long <- inner_join(dprim_long, dpop, by = join_by(lopnr)) %>% 
  mutate(
    alder_dpop = as.period(interval(fodarman, indatum), "years"),
    alder_dpop = as.numeric(alder_dpop, "years")
    )


### Format Primary care -----------------------------------------------------
#' ### Format Primary Care
dprim_long <- dprim_long %>% 
  mutate(diag_letter = str_sub(diag, 1, 1), 
         diag_3 = str_sub(diag, 1, 3))

dprim_long2 <- dprim_long %>% 
  mutate(
    diag_final = ifelse(diag %in% c("L308C", "J310"), diag, diag_3)
  ) %>% 
  filter(diag_final %in% c("J45", "J20", "J21", "J22")) %>% 
  mutate(
    outcome = ifelse(diag_final == "J45", "Asthma", "Wheeze"),
    date_event = as.Date(indatum)
  )

#' Table with ICD-10 codes for each algorithm: 
dlist_prim <- tribble(~ algorithm, ~ icd_codes, 
                     "asth_ort",           c("J45"),
                     "asth_ort3",          c("J45"),
                     "wheeze_j20_j21_j22", c("J45", "J20", "J21", "J22"))


#' Function to create data set wih first event per subject: 
CreateICD10 <- function(data, diagnosis_codes) { 
  d2 <- data %>% 
    filter(diag_final %in% diagnosis_codes) %>% 
    arrange(lopnr, date_event) %>% 
    group_by(lopnr) %>% 
    slice(1) %>% 
    ungroup()
  
  return(d2)
}

#' Add data to table. 
#' The "subj_list" has the first date of ICD10 diagnosis for ALL subjects starting with diagnoses in 2010. 
dlist_prim <- dlist_prim %>% 
  mutate(data_source = "prim") %>% 
  mutate(data = map(icd_codes, ~ CreateICD10(data = dprim_long2, 
                                             diagnosis_codes = .x)))

#' Fix data for ort3: first outcome after 3rd birthday
dlist_prim$data[dlist_prim$algorithm == "asth_ort3"][[1]] <- dprim_long2 %>%
  filter(diag_final == "J45" & alder_dpop > 3) %>%
  arrange(lopnr, date_event) %>%
  group_by(lopnr) %>%
  slice(1) %>%
  ungroup()

#' Filter outcomes to the right ages: 
dlist_prim <- dlist_prim %>% 
  mutate(max_age = ifelse(algorithm == "wheeze_j20_j21_j22", 3, 13)) %>% 
  mutate(data = map2(data, max_age, ~ filter(.x, alder_dpop <= .y))) %>% 
  select(-icd_codes)


# Combine data sets -------------------------------------------------------
#' # Combine data sets

dcomb <- bind_rows(dlist_out, dlist_prim) 

#' Merge in our population list
dcomb2 <- dcomb %>% 
  mutate(data = map(data, ~ select(.x, lopnr, date_event)),
         data = map(data, ~ mutate(.x, date_event = as.Date(as.character(date_event)))))

dcomb3 <- dcomb2 %>% 
  unnest(cols = data) %>% 
  pivot_wider(values_from = date_event, names_from = data_source) %>% 
  group_by(algorithm, max_age) %>% 
  nest() %>% 
  # merge in full population list and remove any outcomes that happen after the censor date (move or data)
  mutate(data = map(data, ~ left_join(dpop2, .x, by = join_by(lopnr))))

dfull <- dcomb3 %>% 
  unnest(cols = data) 

#' Drop events that happen after the censor date
dfull2 <- dfull %>% 
  mutate(
    across(
      c(icd_lmed, prim),
      ~if_else(
        !is.na(.x) & .x > date_censor, 
        as.Date(NA), 
        .x
      )
    )
  )

#' Convert from dates to flags
dfull3 <- dfull %>% 
  mutate(across(c(icd_lmed, prim), 
                ~ !is.na(.x)))


# Calculate Sensitivity and Specificity -----------------------------------
#' # Calculate Sensitivity and Specificity

dconf <- dfull3 %>% 
  pivot_longer(cols = c("icd_lmed"), 
               values_to = "prediction", 
               names_to = "data_source") %>% 
  group_by(algorithm, data_source) %>% 
  nest() 

#' Add confusion matrix information 
dconf <- dconf %>% 
  mutate(conf = map(
    data, 
    ~ DescTools::Conf(x = .x$prediction, ref = .x$prim, pos = "TRUE") # TRUE = event present
  ),          
  conf_table = map(
    conf, ~ .x$table
  ),
  conf_table = map(
    conf_table, as_tibble
  ),
  conf_sum   = map(
    conf, ~.x$byclass %>% as_tibble(., rownames = "metric")
  )
  )

#' Unnest outcomes
dconf2 <- dconf %>% 
  select(algorithm, data_source, conf_sum) %>% 
  unnest(cols = c(conf_sum)) %>% 
  rename(value = `TRUE`)

dconf3 <- dconf2 %>% 
  filter(data_source == "icd_lmed") %>% 
  pivot_wider(names_from = metric, 
              values_from = value)

#' Format outcomes
sens_final <- dconf3 %>% 
  ungroup() %>% 
  rename(outcome = algorithm) %>% 
  mutate(outcome = factor(outcome, 
                          levels = c("asth_ort", "asth_ort3", "wheeze_j20_j21_j22"),
                          labels = c("Asthma", "Asthma (3+)", "Wheeze"))) %>% 
  select(outcome, 
         Sensitivity = sens, 
         Specificity = spec, 
         "Positive Predictive Value" = ppv, 
         "Negative Predictive Value" = npv) %>% 
  mutate(across(where(is.double), round, 2)) 

gt(sens_final)
write_csv2(sens_final, file.path(output_loc, "ts3-outcome-validation.csv"))


## Paper Text --------------------------------------------------------------
#' ## Paper Text 

#' In our outcome validation cohort of 
#' `r length(unique(dfull2$lopnr))` subjects born in Blekinge between 
#' `r min(dfull2$fodelsear)` and 
#' `r max(dfull2$fodelsear)`
#' 
#' Our algorithms performed well with estimated sensitivities between
#' `r round(min(dconf3$sens),2)` and `r round(max(dconf3$sens),2)` 
#' and specificity between 
#' `r round(min(dconf3$spec),2)` and `r round(max(dconf3$spec),2)`



# Check for wheeze_j20_j21_j22 --------------------------------------------
#' # Check for wheeze_j20_j21_j22
#' Reviewer-requested checks not used directly in the manuscript 
# test <- dfull2 %>% filter(algorithm == "wheeze_j20_j21_j22")
# 
# sum(test$prim)
# sum(test$icd_lmed)
# 
# 
# t2 <- filter(test, prim == TRUE)
# sum(t2$icd_lmed) / nrow(t2) * 100
# 
# 
# #' What type of diagnosis codes are missing from t2 (i.e., driving low sensitivity)? 
# t3 <- filter(t2, icd_lmed == FALSE)
# t4 <- semi_join(dprim_long2, t3)
# 
# t4 %>% group_by(diag_final) %>% count() %>% arrange(-n)
# 
# 
# t4 %>% 
#   group_by(diag_final, lopnr) %>% 
#   count() %>% 
#   pivot_wider(names_from = diag_final, values_from = n) %>% 
#   mutate(Jall = sum(J20, J21, J22, na.rm = T)) %>% 
#   mutate(across(everything(), ~ !is.na(.x))) %>% 
#   ungroup() %>% 
#   summarize(J20 = sum(J20),
#             J21 = sum(J21),
#             J22 = sum(J22),
#             J45 = sum(J45),
#             Jall = sum(Jall),
#             n = n())
