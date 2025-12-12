#' ---
#' title: "09 - Matching for RCM"
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
#' This script constructs a matched cohort comparing children in the Very High vs. Low (background) prenatal exposure groups, 
#' using the same complete-cases as the primary analysis. 
#' It then restricts the outcome-specific datasets to the matched subjects 
#' and produces a covariate balance plot (Love plot).   
#'  
#' ## Inputs 
#' 
#' - `cohort-outcome-data.rds`
#' 
#' ## Outputs
#' 
#' - `fs7-love-plot.tiff`
#' - `dvh_matched_out_list.rds` - list-column tibble with one matched dataset per outcome algorithm. 
#' 

# Setup -------------------------------------------------------------------
#' # Setup

#+ load-libraries, message = F
library(here)
library(tidyverse)
library(gt)
library(survival)
library(survminer)
library(MatchIt)

output_loc <- here("outputs_public")

sessionInfo()

set.seed(314)

# Create Matching Data ----------------------------------------------------
#' # Create Matching Data 

dlist <- read_rds(file.path(output_loc, "cohort-outcome-data.rds"))
dlist

#' For convenience, take the first outcome's dataset as the base cohort. 
#' All outcomes have the same base cohort information. 

dpop <- dlist$data[[1]] %>% 
  filter(flag_missing_covariates == FALSE)

#' Matching Datasets 
dvh <- dpop %>% 
  filter(ecat_prenat_4 %in% c("Low", "Very High")) %>% 
  mutate(ecat_prenat_4 = fct_drop(ecat_prenat_4))


nn_rmah_31_dvh <- matchit(ecat_prenat_4 ~ rok1_2cat + mat_educ_3cat + mat_age_delivery + dispinkfam_comb + 
                            paritet_3cat + flag_parent_born_abroad + kon + flag_asth_ort_mf,
                          data = dvh, 
                          method = "nearest",
                          distance = "robust_mahalanobis",
                          ratio = 3)

#' Tibble with matched subjects and matching group
dvh_matched <- match.data(nn_rmah_31_dvh) %>%
  select(lopnr, subclass, weights, ecat_prenat_4)

#' Limit full cohort data for each outcome to the matched subjects
dlist_vh <- dlist %>% 
  mutate(
    data = map(data, 
               ~ inner_join(.x, dvh_matched, by = join_by(lopnr, ecat_prenat_4)))
  )

write_rds(dlist_vh, file.path(output_loc, "dvh_matched_out_list.rds"))


# Love Plot ---------------------------------------------------------------
#' # Love Plot

#' Variable labels: 
var_labels <- tribble(~ var, ~ var_label,
          "rok1_2catNon-smoker",                      "Maternal smoking: Non-smoker", 
          "rok1_2catSmoker",                          "Maternal smoking: Smoker", 
          "mat_educ_3catPrimary and lower secondary", "Maternal education: Primary and lower secondary",
          "mat_educ_3catUpper secondary",             "Maternal education: Upper secondary",
          "mat_educ_3catPost secondary",              "Maternal education: Post secondary",
          "mat_age_delivery",                         "Maternal age at delivery",    
          "dispinkfam_comb",                          "Family disposable income",  
          "paritet_3cat1",                            "Parity: first-born",             
          "paritet_3cat2",                            "Parity: second-born",  
          "paritet_3cat3+",                           "Parity: third-born or later",    
          "flag_parent_born_abroadTRUE",              "At least one parent born abroad",         
          "konmale",                                  "Male",          
          "konfemale",                                "Female",       
          "flag_asth_ort_mfTRUE",                     "Parental Asthma",
          "flag_asth_ort_mf",                         "Parental Asthma",
          "flag_parent_born_abroad",                  "At least one parent born abroad", 
          "kon",                                      "Child Sex", 
          "mat_educ_3cat",                            "Maternal Education", 
          "paritet_3cat",                             "Parity", 
          "rok1_2cat",                                "Maternal Smoking in Early Pregnancy")

#' Standardized Mean Difference: 
msum_vh <- summary(nn_rmah_31_dvh)$sum.matched %>%
  as_tibble(rownames = "var") %>% 
  janitor::clean_names() %>% 
  mutate(abs_smd = abs(std_mean_diff))

#' Love Plot:
MakeLovePlot <- function(matchit_object) {
  
  msum <- summary(matchit_object)
  
  
  msum_prematch <- as_tibble(msum$sum.all, rownames = "var")
  msum_prematch <- janitor::clean_names(msum_prematch) %>% 
    mutate(dataset = "Pre-Match",
           abs_smd = abs(std_mean_diff))  %>% 
    left_join(var_labels, by = join_by(var)) 
  
  msum_match <- as_tibble(msum$sum.matched, rownames = "var")
  msum_match <- janitor::clean_names(msum_match) %>% 
    mutate(dataset = "Matched",
           abs_smd = abs(std_mean_diff))  %>% 
    left_join(var_labels, by = join_by(var))
  
  #' Order factor variable: 
  ordered_var <- msum_prematch %>% 
    filter(!is.na(var_label)) %>% 
    arrange(abs_smd) %>% 
    pull(var_label) %>% 
    unique()
  
  msum_all <- bind_rows(msum_prematch, msum_match) %>% 
    mutate(var_label = factor(var_label, levels = unique(ordered_var)))
  
  p1 <- ggplot(msum_all) + 
    geom_point(aes(x = var_label, y = abs_smd, shape = dataset, fill = dataset)) +
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = 0.1) + 
    geom_hline(yintercept = 0.05, linetype = "dashed") + 
    coord_flip() + 
    scale_shape_manual(values = c(21, 21)) + 
    scale_fill_manual(values= c("black", "white")) + 
    labs(y = "Absolute Standardized Mean Difference",
         shape = "", 
         fill = "") + 
    theme_bw() + 
    theme(axis.title.y = element_blank(),
          legend.position = "bottom")
  
  return(p1)
  
}

p1 <- MakeLovePlot(nn_rmah_31_dvh) + 
  labs(title = "Pre- and post-match covariate balance") + 
  ylim(c(0, 0.52))

p1

#' Save: 
ggsave(filename = "fs7-love-plot.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       path = output_loc,
       dpi = 600,
       width = 1877*2, 
       height = 1410*2, 
       units = "px")

# Paper Text --------------------------------------------------------------
#' # Paper Text

#' The matched dataset included 
#' `r sum(dvh_matched$ecat_prenat_4 == "Very High")` 
#' very highly exposed and 
#' `r sum(dvh_matched$ecat_prenat_4 == "Low")` 
#' background exposed.
#' 
#' 
#' The cohort showed excellent covariate balance, with a maximum standardized mean difference of 
#' `r signif(max(msum_vh$abs_smd), 2)`.

