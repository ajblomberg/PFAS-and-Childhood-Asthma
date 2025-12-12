#' ---
#' title: "06 - Primary and sex-stratified survival models"
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
#' This constructs the main models used in this paper, including overall (unadjusted and fully adjusted), sex-stratified and sex-interaction models.
#' It creates the summary tables and figures and also runs the proportional-hazard assumptions.  
#' 
#' ## Inputs 
#' 
#' - `cohort-outcome-data.rds` (main cohort with outcome-specific follow-up time)
#' 
#' ## Outputs
#' - `t3-primary-model-results.csv`
#' - `results-primary-models.rds`
#' - `ts6-hazard-ratio-by-sex.csv`


# Setup -------------------------------------------------------------------
#' # Setup 

#+ load-libraries, message = F
library(here)
library(tidyverse)
library(gt)
library(survival)
library(broom) 
library(car)  

output_loc <- here("outputs_public")

sessionInfo()


# Load Data ---------------------------------------------------------------
#' # Load Data 
dlist <- read_rds(file.path(output_loc, "cohort-outcome-data.rds")) %>% 
  mutate(data = map(data, ~ filter(.x, flag_missing_covariates == FALSE)))

# Create Models  ----------------------------------------------------------
#' # Create Models 
dlist_male <- dlist %>% 
  mutate(sex = "Male") %>% 
  mutate(data = map(data, ~ filter(.x, kon == "male")))

dlist_female <- dlist %>% 
  mutate(sex = "Female") %>% 
  mutate(data = map(data, ~ filter(.x, kon == "female")))

dlist_comb <- dlist %>% 
  mutate(sex = "Combined")

dlist_sex <- bind_rows(dlist_male, dlist_female, dlist_comb)

#' Note on sex-interaction model: 
#' 
#' The "sex_int" specification uses ecat_prenat_4:strata(kon) as the sex interaction. 
#'  
#' In coxph, interacting with strata(kon) is a supported way to obtain sex-specific exposure coefficients 
#' while keeping sex-specific baseline hazards, without creating aliased interaction terms. 

model_specs_sex <- tribble(~ sex, ~ model_label, ~ covariates, 
                           "Male",     "full",    "ecat_prenat_4 + rok1_2cat + mat_educ_3cat + mat_age_delivery_quant + dispinkfam_quant + strata(paritet_2cat) + flag_parent_born_abroad + flag_asth_ort_mf",
                           "Female",   "full",    "ecat_prenat_4 + rok1_2cat + mat_educ_3cat + mat_age_delivery_quant + dispinkfam_quant + strata(paritet_2cat) + flag_parent_born_abroad + flag_asth_ort_mf",
                           "Combined", "unadj",   "ecat_prenat_4 + strata(kon)",
                           "Combined", "full",    "ecat_prenat_4 + rok1_2cat + mat_educ_3cat + mat_age_delivery_quant + dispinkfam_quant + strata(paritet_2cat) + flag_parent_born_abroad + flag_asth_ort_mf + strata(kon)",
                           "Combined", "sex_int", "ecat_prenat_4 + rok1_2cat + mat_educ_3cat + mat_age_delivery_quant + dispinkfam_quant + strata(paritet_2cat) + flag_parent_born_abroad + flag_asth_ort_mf + strata(kon) + ecat_prenat_4:strata(kon)")

dlist_sex <- left_join(dlist_sex, model_specs_sex, by = join_by(sex), relationship = "many-to-many")

#' Add survival models and PYs
mlist <- dlist_sex %>% 
  mutate(
    model = map2(
      data, covariates, 
      ~ coxph(
        as.formula(paste("Surv(age_final, flag_event) ~", .y)),
        cluster = lopnr_mor,
        data = .x
      )
    ),
    model_est = map(model, ~ broom::tidy(.x, conf.int = T, exp = T)),
    pyears    = map(
      data, 
      ~ pyears(Surv(time = .x$age_final, event = .x$flag_event) ~ .x$ecat_prenat_4, data = .x, scale = 1)
    ),
    pyears = map(pyears, broom::tidy),
    pyears = map(pyears, ~ mutate(.x, exposure = c("Low", "Intermediate", "High", "Very High"))
    )
  ) 


# Primary Model Effects ---------------------------------------------------
#' # Primary Model Effects 
mlist1 <- mlist %>% 
  filter(sex == "Combined" & model_label %in% c("full", "unadj")) %>% 
  mutate(model_label = factor(model_label, levels = c("unadj", "full"), labels = c("Unadjusted", "Adjusted"))) %>% 
  arrange(algorithm_label)

mlist2 <- mlist1 %>% 
  select(algorithm_label, model_label, model_est) %>% 
  unnest(cols = c(model_est)) %>% 
  filter(term %in% c("ecat_prenat_4Intermediate", "ecat_prenat_4High", "ecat_prenat_4Very High")) %>% 
  mutate(exposure = str_extract(term, "Intermediate|High|Very High"), 
         exposure = factor(exposure, levels = c("Intermediate", "High", "Very High"))) 

mlist3 <- mlist1 %>% 
  select(algorithm_label, model_label, pyears) %>% 
  unnest(cols = c(pyears)) %>% 
  select(-model_label) %>% 
  unique()

#' Summary Table 
t1 <- mlist2 %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ round(., 2))) %>% 
  mutate(est_format = paste0(estimate, " (", conf.low, ", ", conf.high, ")"),
         est_format = ifelse(is.na(estimate), "-", est_format)) %>% 
  select(algorithm_label, model_label, exposure, est_format) %>% 
  pivot_wider(values_from = est_format, names_from = model_label) 

t2 <- mlist3 %>% 
  full_join(t1, by = join_by(algorithm_label, exposure)) %>% 
  mutate(across(c(Unadjusted, Adjusted), ~ ifelse(is.na(.x), "-", .x))) %>% 
  select(algorithm_label, exposure, n, event, pyears, Unadjusted, Adjusted) %>% 
  mutate(pyears = round(pyears, 0), 
         nevents_perc = round(event / n * 100,0),
         event = paste0(event, " (", nevents_perc, "%)")) %>%
  select(-nevents_perc) %>% 
  group_by(algorithm_label) %>% 
  rename(Exposure = exposure, 
         N = n,
         Events = event, 
         'PYs' = pyears) 

gt(t2)
write_csv2(t2,    file.path(output_loc, "t3-primary-model-results.csv"))
write_rds(mlist2, file.path(output_loc, "results-primary-models.rds"))

## Covariate Effect Estimates  ---------------------------------------------
#' ## Covariate Effect Estimates
covs <- mlist %>% 
  filter(sex == "Combined", 
         model_label == "full") %>% 
  select(algorithm_label, model_est) %>% 
  unnest(cols = model_est) %>% 
  filter(! (term %in% c("ecat_prenat_4Intermediate", "ecat_prenat_4High", "ecat_prenat_4Very High"))) %>% 
  mutate(term = factor(term, 
                       levels = c("rok1_2catSmoker", "mat_educ_3catUpper secondary", "mat_educ_3catPost secondary", 
                                  "mat_age_delivery_quant(26.3,30]", "mat_age_delivery_quant(30,33.7]", "mat_age_delivery_quant(33.7,50.7]", 
                                  "dispinkfam_quant(2.81e+03,3.72e+03]", "dispinkfam_quant(3.72e+03,4.51e+03]", "dispinkfam_quant(4.51e+03,3.21e+04]", 
                                  "paritet_2catMultiparous", "flag_parent_born_abroadTRUE", "flag_asth_ort_mfTRUE"
                       ),
                       labels = c("Maternal smoking in early pregnancy", "Maternal education: Upper Secondary", "Maternal education: Post Secondary",
                                  "Maternal age: Q2 (26.3,30]", "Maternal age: Q3 (30,33.7]", "Maternal age: Q4 (33.7,50.7]", 
                                  "Income: Q2 (2.81e+03,3.72e+03]", "Income: Q3 (3.72e+03,4.51e+03]", "Income: Q4 (4.51e+03,3.21e+04]",
                                  "Multiparous", "Parent born abroad", "Parental asthma")))

ggplot(covs, aes(x = term, y = estimate, color = algorithm_label, shape = algorithm_label)) + 
  geom_point(position = position_dodge(width = 0.6)) + 
  geom_linerange(aes(ymin = conf.low, ymax = conf.high),
                 position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_color_brewer(type = "discrete", palette = "Dark2") + 
  labs(x = "Covariate", 
       y = "Hazard ratio (95% CI)",
       color = "Outcome", 
       shape = "Outcome") + 
  scale_y_log10(breaks = c(0.75, 1.0, 1.5, 2.0, 2.5),
                labels = c(0.75, 1.0, 1.5, 2.0, 2.5)) + 
  coord_flip() + 
  theme_bw() +
  guides(color = guide_legend(reverse = T),
         shape = guide_legend(reverse = T))

#' Save: 
ggsave(filename = "fs4-primary-model-covariate-estimates.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       path = output_loc,
       dpi = 600,
       width = 1877*2, 
       height = 1877*2, 
       units = "px")


# Sex Models --------------------------------------------------------------
#' # Sex Models 
CalculateRobustWald <- function(m_full, term){
  terms_full <- names(coef(m_full))

  L_idx <- which(str_detect(terms_full, term))

  car::linearHypothesis(
    model = m_full,
    hypothesis.matrix = diag(length(coef(m_full)))[L_idx, , drop = FALSE],
    vcov. = vcov(m_full)
  )
}


list_mwald <- mlist %>%
  filter(model_label == "sex_int")

# Terms with ":strata" are the sex interaction terms
list_robust_wald <- map(list_mwald$model, CalculateRobustWald, term = ":strata")
names(list_robust_wald) <- list_mwald$outcome_algo
list_robust_wald


## Sex-Stratified Estimates  -----------------------------------------------
#' ## Sex-Stratified Estimates 
p1 <- mlist %>% 
  filter(model_label == "full") %>% 
  select(c("algorithm_label", "outcome_algo", "model_est", "sex")) %>% 
  mutate(sex = factor(sex, levels = c("Combined", "Male", "Female"),
                      labels = c("Overall", "Male", "Female"))) %>% 
  unnest(cols = c(model_est)) %>% 
  filter(term %in% c("ecat_prenat_4Intermediate", "ecat_prenat_4High", "ecat_prenat_4Very High")) %>%
  mutate(term = factor(term, levels = c("ecat_prenat_4Intermediate", "ecat_prenat_4High", "ecat_prenat_4Very High"),
                       labels = c("Intermediate", "High", "Very High"))) %>% 
  ggplot(aes(x = term, y = estimate, color = sex, shape = sex)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(x = term, ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  scale_color_brewer(type = "discrete", palette = "Dark2") + 
  facet_wrap(~ algorithm_label, ncol = 3) + 
  labs(x = "Prenatal Exposure Group",
       y = "Hazard Ratio (95% CI)", 
       color = "Stratified models",
       shape = "Stratified models") + 
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1.0, 2.0),
                labels = c(0.1, 0.2, 0.5, 1.0, 2.0)) + 
  theme_bw() + 
  theme(legend.position = "bottom")

p1

ggsave(filename = "fs6-sex-stratified-estimates.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       path = output_loc,
       dpi = 600,
       width = 1877*2, 
       height = 940*2, 
       units = "px")

# Create a table with sex-stratified results 
df_est <- mlist %>% 
  filter(model_label =="full" & sex %in% c("Male", "Female")) %>%
  arrange(algorithm_label) %>% 
  mutate(algorithm_label = paste0(algorithm_label, "-", sex)) %>% 
  select(algorithm_label, model_est) %>% 
  unnest(cols = c(model_est)) %>% 
  filter(term %in% c("ecat_prenat_4Intermediate", "ecat_prenat_4High", "ecat_prenat_4Very High")) %>% 
  mutate(exposure = str_extract(term, "Intermediate|High|Very High"), 
         exposure = factor(exposure, levels = c("Intermediate", "High", "Very High"))) %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~ round(., 2))) %>% 
  mutate(est_format = paste0(estimate, " (", conf.low, ", ", conf.high, ")"),
         est_format = ifelse(is.na(estimate), "-", est_format)) %>% 
  select(algorithm_label, exposure, est_format)


df_py <- mlist %>% 
  filter(model_label =="full" & sex %in% c("Male", "Female")) %>%
  arrange(algorithm_label) %>% 
  mutate(algorithm_label = paste0(algorithm_label, "-", sex)) %>% 
  select(algorithm_label, model_label, pyears) %>% 
  unnest(cols = c(pyears)) %>% 
  select(-model_label) %>% 
  unique()

df_comb <- df_py %>% 
  full_join(df_est, by = join_by(algorithm_label, exposure)) %>% 
  select(algorithm_label, exposure, n, event, pyears, est_format) %>% 
  mutate(pyears = round(pyears, 0), 
         nevents_perc = round(event / n * 100,0),
         event = paste0(event, " (", nevents_perc, "%)")) %>%
  select(-nevents_perc) %>% 
  group_by(algorithm_label) %>% 
  rename(Exposure = exposure, 
         N = n,
         Events = event, 
         'PYs' = pyears) 

gt(df_comb)
write_csv2(df_comb, file.path(output_loc, "ts6-hazard-ratio-by-sex.csv"))

# Proportional Hazards Assumption -----------------------------------------------------------
#' # Proportional Hazards Assumption

#' Check interaction between ecat_prenat_4 and time in combined models: 

#+ identify-models-for-ph-check
mlist_ph <- mlist %>% 
  filter(sex == "Combined" & model_label == "full")

#' Convert factor variable to binary indicators: 
#+ create-time-interaction
mlist_ph$data <- map(mlist_ph$data, 
                     ~ mutate(.x, 
                              e4_int = ecat_prenat_4 == "Intermediate",
                              e4_hi  = ecat_prenat_4 == "High",
                              e4_vh  = ecat_prenat_4 == "Very High"))

#+ fit-time-interaction-models
list_fit_tv <- map(mlist_ph$data, 
                   ~ coxph(Surv(time = age_final, event = flag_event) ~ 
                             e4_int + e4_hi + e4_vh + 
                             tt(e4_int) + tt(e4_hi) + tt(e4_vh) + 
                             rok1_2cat + mat_educ_3cat + mat_age_delivery_quant + dispinkfam_quant +  
                             flag_parent_born_abroad + flag_asth_ort_mf + strata(kon) + strata(paritet_2cat),
                           cluster = lopnr_mor,
                           data = .x,
                           tt = function(x, t, ...) x * t))
 
#' Test the interaction terms ("tt") using a Robust Wald test here as well:  
list_robust_wald <- map(list_fit_tv, CalculateRobustWald, term = "tt")
names(list_robust_wald) <- mlist_ph$outcome_algo
list_robust_wald
