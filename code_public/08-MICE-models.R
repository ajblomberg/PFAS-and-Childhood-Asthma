#' ---
#' title: "08 - MICE models"
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
#' This script fits Cox models in the MICE imputed datasets (asthma, asthma 3+, wheeze), 
#' pools the results, and compares them to the primary complete-case models. 
#'  
#' ## Inputs 
#' 
#' - `mice-imp-20-list.rds` (list of 3 `mids` objects)
#' - `results-primary-models.rds` (results from primary models)
#' 
#' ## Outputs
#'
#' - `ts5-primary-vs-mice-model-results.csv`
#' - `fs5-compare-sensitivity-models.tiff`


# Setup -------------------------------------------------------------------
#' # Setup

#+ load-libraries, message = F
library(tidyverse)
library(survival)
library(here)
library(mice)

output_loc <- here("outputs_public")

sessionInfo()


# Run Models --------------------------------------------------------------
#' # Run Models 
#' 
#' Imputations have already been completed - we load theme here. 
#+ save or load imputations 
imp_list <- read_rds(file.path(output_loc, "mice-imp-20-list.rds"))
length(imp_list)


#' Add important variables using a helper function: 
imp_long_list <- map(imp_list, ~ complete(.x, "long")) %>% 
  map(., ~ as_tibble(.x)) %>% 
  map(., ~ mutate(.x, rok1_2cat = fct_collapse(rok1, "Smoker" = c("1-9 cig/day", "10+ cig/day")), 
                  mat_educ_3cat = fct_collapse(sun_old_cat, 
                                               "Primary and lower secondary" = c("Pre-secondary education shorter than 9 years", 
                                                                                 "Pre-secondary education equivalent to 9 years",
                                                                                 "Upper-secondary education, no more than 2 years"),
                                               "Upper secondary" = c("Secondary education, 3 years"),
                                               "Post secondary"  = c("Post-secondary education shorter than 3 years",
                                                                     "Post-secondary education 3 years or longer", 
                                                                     "Postgraduate education")),
                  paritet_2cat = factor(paritet), 
                  paritet_2cat = fct_collapse(paritet_2cat, 
                                              "Multiparous" = c(2, 4, 6, 7, 3, 5, 9, 10, 8, 12)),
                  paritet_2cat = fct_recode(paritet_2cat, 
                                            Primiparous = "1"),
                  flag_parent_born_abroad = utlsvbakg != "native-born with native-born parents",
                  mat_age_delivery_quant = cut(mat_age_delivery, 
                                               include.lowest = TRUE, 
                                               breaks = quantile(mat_age_delivery, 
                                                                 probs = c(0, 0.25, 0.5, 0.75, 1))),
                  dispinkfam_quant = cut(dispinkfam_comb, 
                                         include.lowest = TRUE,
                                         breaks = quantile(dispinkfam_comb,
                                                           probs = c(0, 0.25, 0.5, 0.75, 1)))))


## Run Models ----------------------------------------------------
#' ## Run Models 

#' Helper functions: fit one Cox model per imputed dataset (.imp) 
#' and pool results with mice::pool(). 
#' Returns a pooled model object. 
RunAdjSurvival <- function(data) { 
  est <- data %>% 
    group_by(.imp) %>% 
    do(model = 
         coxph(Surv(time = age_final, event = flag_event) ~ 
                 ecat_prenat_4 + rok1_2cat + mat_educ_3cat + mat_age_delivery_quant + dispinkfam_quant + flag_parent_born_abroad + flag_asth_ort_mf + strata(kon, paritet_2cat), 
               cluster = lopnr_mor, 
               data = .)) %>% 
    as.list() %>% 
    .[[-1]] %>% 
    pool()
}

RunUnadjSurvival <- function(data) { 
  est <- data %>% 
    group_by(.imp) %>% 
    do(model = 
         coxph(Surv(time = age_final, event = flag_event) ~ 
                       ecat_prenat_4 + strata(kon), 
                     cluster = lopnr_mor, 
                     data = .)) %>% 
    as.list() %>% 
    .[[-1]] %>% 
    pool()
}

list_adj_model   <- map(imp_long_list, ~ RunAdjSurvival(.x))
list_unadj_model <- map(imp_long_list, ~ RunUnadjSurvival(.x))


list_adj_model_summary <- map(list_adj_model, 
                              ~ summary(.x, conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
                                as_tibble() %>% 
                                janitor::clean_names())

list_unadj_model_summary <- map(list_unadj_model, 
                                ~ summary(.x, conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
                                  as_tibble() %>% 
                                  janitor::clean_names())



adj_model_summary   <- bind_rows(list_adj_model_summary, .id = "algorithm")
unadj_model_summary <- bind_rows(list_unadj_model_summary, .id = "algorithm")


## Summarize MICE ----------------------------------------------------------
#' ## Summarize MICE

model_summary_mice <- bind_rows(unadj = unadj_model_summary, 
                                adj   = adj_model_summary, 
                                .id = "model_label") %>% 
  filter(term %in% c("ecat_prenat_4Intermediate", "ecat_prenat_4High", "ecat_prenat_4Very High")) %>% 
  mutate(exposure = str_extract(term, "Intermediate|High|Very High")) %>% 
  mutate(model_label_mice = "MICE") 


#' Counts of events and person-years are identical across imputations, 
#' so use the first completed dataset (.imp == 1) for pyears. 
pyear_list <- map(imp_long_list, ~ filter(.x, .imp == 1)) %>% 
  map(., ~ pyears(Surv(time = .x$age_final, event = .x$flag_event) ~ .x$ecat_prenat_4, 
                  scale = 1, 
                  data = .x)) %>% 
  map(., broom::tidy) %>% 
  map(., ~ mutate(.x, exposure = c("Background", "Intermediate", "High", "Very High")))

pyears <- bind_rows(pyear_list, .id = "outcome") %>% 
  mutate(pyears = round(pyears, 0)) %>% 
  select(algorithm = outcome, exposure, event, pyears) %>% 
  mutate(algorithm = factor(algorithm, levels = c("wheeze", "asthma", "asthma3"))) %>% 
  arrange(algorithm)


model_summary_mice2 <- model_summary_mice %>% 
  mutate(across(c(estimate, conf_low, conf_high), ~ round(., 2))) %>% 
  mutate(est_format = paste0(estimate, " (", conf_low, ", ", conf_high, ")"),
         est_format = ifelse(is.na(estimate), "-", est_format)) %>% 
  mutate(algorithm = factor(algorithm, 
                            levels = c("wheeze", "asthma", "asthma3"))) %>% 
  select(algorithm, model_label, exposure, est_format) %>% 
  pivot_wider(values_from = est_format, names_from = model_label) %>% 
  right_join(pyears, by = join_by(algorithm, exposure)) %>% 
  arrange(algorithm, exposure) %>% 
  select(algorithm, exposure, event, pyears, unadj, adj)
  

gt(model_summary_mice2)
write_csv2(model_summary_mice2, 
           file.path(output_loc, "ts5-primary-vs-mice-model-results.csv"))



# Plot Model Comparison  --------------------------------------------------
#' # Plot Model Comparison 
results_primary <- read_rds(file.path(output_loc, "results-primary-models.rds")) %>% 
  mutate(model_type = "Primary") %>% 
  select(model_adj = model_label, model_type, algorithm_label, exposure, estimate, conf_low= conf.low, conf_high = conf.high)


results_mice <- model_summary_mice %>% 
  rename(model_adj = model_label) %>% 
  mutate(model_adj = factor(model_adj, 
                            levels = c("unadj", "adj"),
                            labels = c("Unadjusted", "Adjusted"))) %>% 
  mutate(algorithm_label = factor(algorithm, 
                                  levels = c("wheeze", "asthma", "asthma3"),
                                  labels = c("Wheeze", "Asthma", "Asthma (3+)"))) %>% 
  mutate(model_type = "MICE") %>% 
  select(model_adj, model_type, algorithm_label, exposure, estimate, conf_low, conf_high)


bind_rows(results_primary, results_mice) %>% 
  filter(model_adj == "Adjusted") %>% 
  mutate(exposure = factor(exposure, 
                           levels = c("Intermediate", "High", "Very High")),
         algorithm_label = factor(algorithm_label, 
                                  levels = c("Wheeze", "Asthma", "Asthma (3+)"))) %>% 
  ggplot(aes(x = exposure, y = estimate, color = model_type, shape = model_type)) + 
  geom_point(position = position_dodge(width = 0.6)) + 
  geom_linerange(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  facet_wrap( ~ algorithm_label) + 
  scale_color_brewer(type = "discrete", palette = "Dark2") +
  scale_y_log10(breaks = c(0.5, 0.75, 1.0, 1.5, 2.0),
                labels = c(0.5, 0.75, 1.0, 1.5, 2.0)) + 
  labs(x = "Prenatal Exposure Group",
       y = "Hazard Ratio (95% CI)",
       color = "Model Type", 
       shape = "Model Type") +
  theme_bw()  + 
  theme(legend.position = "bottom")

ggsave(filename = "fs5-compare-sensitivity-models.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       path = output_loc,
       dpi = 600,
       width = 1877*2, 
       height = 940*2, 
       units = "px")
  
  
  