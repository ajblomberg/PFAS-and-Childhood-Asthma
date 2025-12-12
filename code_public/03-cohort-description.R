#' ---
#' title: "03 - Cohort Description and Methods Summaries"
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
#' This script produces descriptive summaries and figures used in the manuscript's Methods and Results sections. 
#' It: 
#' 
#' - Summarizes the study cohort 
#' - Describes covariate completeness and baseline characteristics 
#' - Verifies address-based exposure categories using biomonitoring data
#' - Summarizes outcome prevalence and overlap between outcomes. 
#' - Provides manuscript-ready text via inline R and software version information. 
#' 
#' ## Inputs 
#' 
#' - `cohort-outcome-data.rds`: main analysis cohort with one row per outcome algorithm 
#' and a list-column of subject-level data. 
#' - `data_exp_validation_cohort.rds`: biomonitoring validation cohort for exposure verification. 
#' 
#' ## Outputs 
#' 
#' - `ts2-complete-vs-missing-covariates.csv`
#' - `t1-baseline-covariates.csv`
#' - `t2-exposure-verification.csv`
#' - `fs1-exp-verification.tiff`
#' - `fs2-outcome-combos.tiff`
#' - Inline text used directly in the manuscript. 


# Setup -------------------------------------------------------------------
#' # Setup 

#' Load core packages 
#+ load-libraries, message = F
library(here)
library(tidyverse)
library(gt)
library(tableone)

#' Optional visualization packages 
#+ load-opt-libraries, message = F
library(ggrepel)
library(DescTools)
library(UpSetR)

#+ load-internal-libraries, message = F
path.library <- "\\\\lces1113cs.srv.lu.se/epidemiologi/Fogrupp/rlibrary"
library(ggbeeswarm, lib.loc = path.library)

output_loc <- here("outputs_public")

sessionInfo()


# Load Main Cohort ---------------------------------------------------------------
#' # Load Main Cohort

dlist <- read_rds(file.path(output_loc, "cohort-outcome-data.rds"))
dlist

#' For convenience, take the first outcome's dataset as the base cohort. 
#' All outcomes have the same base cohort information. 

dpop <- dlist$data[[1]] %>% 
  mutate(dispinkfam_comb = dispinkfam_comb / 1000)


# Study Cohort & Deaths ---------------------------------------------------
#' # Study Cohort & Deaths
#' 
#' Count of deaths in the cohort: 
deaths_all <- map_dbl(dlist$data, ~ sum(.x$censor_category == "died"))
deaths_complete_cases <- dlist$data %>% 
  map( ~ filter(.x, flag_missing_covariates == FALSE)) %>% 
  map_dbl(., ~ sum(.x$censor_category == "died"))

names(deaths_all) <- dlist$algorithm_label
names(deaths_complete_cases) <- dlist$algorithm_label

deaths_all
deaths_complete_cases


# Missing Covariate Data  -------------------------------------------------
#' # Missing Covariate Data 
#' 
#' Summarize extent of missingness in covariates and describe characteristics of children 
#' with and without missing covariate data. 
#' 
#' Compute overall and complete-case counts: 
nfull <- nrow(dpop)
nlim  <- nrow(filter(dpop, flag_missing_covariates == FALSE))


#' Build dataset for Table S2 (complete vs. missing covariates)
dt <- dpop %>%
  mutate(across(
    c(rok1_2cat, paritet_2cat, kon, mat_educ_3cat), 
    ~ fct_na_value_to_level(., level = "Missing"))) %>% 
  select(c("lopnr", 
           "lopnr_mor",
           `Maternal smoking status` = "rok1_2cat", 
           `Parity` = "paritet_2cat", 
           `Child sex` = "kon", 
           `Maternal education` = "mat_educ_3cat", 
           `Maternal age at delivery` = "mat_age_delivery", 
           `Family disposable income` = "dispinkfam_comb",
           `Parent born abroad`    = "flag_parent_born_abroad",
           `Parental asthma` = "flag_asth_ort_mf", 
           "Prenatal exposure group" = "ecat_prenat_4",
           "Missing covariates" = flag_missing_covariates)
         )

covariates <- c(
  "Maternal smoking status", 
  "Parity", 
  "Child sex", 
  "Maternal education", 
  "Maternal age at delivery", 
  "Parent born abroad", 
  "Family disposable income", 
  "Parental asthma", 
  "Prenatal exposure group"
)

t1_missing_overall <- CreateTableOne(
  vars = covariates,
  data = dt,
  test = FALSE,
  includeNA = TRUE
) %>% 
  print(
    showAllLevels = TRUE,
    nonnormal = c("Maternal age at delivery", "Family disposable income"),
    quote = FALSE,
    noSpaces = TRUE,
    printToggle = FALSE,
    explain = FALSE,
    catDigits = 1,
    contDigits = 1
  ) %>%
  as_tibble(rownames = "Variable")


t1_missing_strat <- CreateTableOne(
  vars = covariates,
  data = dt,
  test = FALSE,
  includeNA = TRUE,
  strata = "Missing covariates"
) %>% 
  print(
    showAllLevels = TRUE,
    nonnormal = c("Maternal age at delivery", "Family disposable income"),
    quote = FALSE,
    noSpaces = TRUE,
    printToggle = FALSE,
    explain = FALSE,
    catDigits = 1, 
    contDigits = 1
  ) %>%
  as_tibble(rownames = "Variable")

t1_missing <- bind_cols(t1_missing_overall, t1_missing_strat[, 3:4])

#' Create and save Table S2
gt(t1_missing) %>% 
  tab_header(
    title = "Covariate Summary",
    subtitle = "Complete vs missing covariates"
  )

write_csv2(t1_missing, file.path(output_loc, "ts2-complete-vs-missing-covariates.csv"))

#' Summarize variable-specific missingness
dmiss <- dpop %>% 
  summarize(across(
    c(
      "rok1_2cat", "paritet_2cat", "kon", 
      "mat_educ_3cat", "dispinkfam_comb", 
      "mat_age_delivery", "flag_asth_ort_mf"
    ), 
    ~ sum(is.na(.)))) %>% 
  pivot_longer(
    cols = everything(), 
    names_to = "variable", 
    values_to = "nmiss"
  ) %>% 
  arrange(desc(nmiss)) %>% 
  mutate(perc_miss = round(nmiss / nfull * 100, 2))

gt(dmiss, caption = "Summary of missing variables")

# Baseline Characteristics - Complete Cases ------------------------------------------------
#' # Baseline Characteristics - Complete Cases

dt_lim <- dt %>% 
  filter(`Missing covariates` == FALSE) %>% 
  mutate(across(
    c("Maternal smoking status", "Parity", "Child sex", 
      "Maternal education","Prenatal exposure group"), 
    ~ fct_drop(.)
  ))

covariates <- c("Maternal smoking status", 
                "Parity", 
                "Child sex", 
                "Maternal education", 
                "Maternal age at delivery", 
                "Parent born abroad", 
                "Family disposable income", 
                "Parental asthma", 
                "Prenatal exposure group")

t1_overall <- CreateTableOne(
  vars = covariates,
  data = dt_lim,
  test = FALSE
) %>% 
  print(
    showAllLevels = FALSE,
    missing = FALSE,
    nonnormal = c("Maternal age at delivery", "Family disposable income"),
    quote = FALSE,
    noSpaces = TRUE,
    printToggle = FALSE,
    explain = FALSE,
    catDigits = 1,
    contDigits = 1
  ) %>%
  as_tibble(
    rownames = "Variable"
  )


t1_strat <- CreateTableOne(
  vars = covariates,
  data = dt_lim,
  test = FALSE, 
  strata = "Prenatal exposure group"
) %>% 
  print(
    showAllLevels = FALSE,
    nonnormal = c("Maternal age at delivery", "Family disposable income"),
    quote = FALSE,
    noSpaces = TRUE,
    printToggle = FALSE,
    explain = FALSE,
    catDigits = 1,
    contDigits = 1
  ) %>%
  as_tibble(
    rownames = "Variable"
  )

t1 <- bind_cols(t1_overall, t1_strat[, 2:5])

gt(t1)  %>% 
  tab_header(
    title = "Covariate Summary",
    subtitle = "Limited to complete cases"
  ) 

write_csv2(t1, file.path(output_loc, "t1-baseline-covariates.csv"))


# Exposure Groups: Counts and Moving --------------------------------------
#' # Exposure Groups: Counts and Moving 

n_exp <- dt_lim %>% 
  group_by(`Prenatal exposure group`) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(percent = round(n / nrow(dt_lim) * 100, 0))

gt(n_exp)



## Exposure over Childhood -------------------------------------------------
#' ## Exposure over Childhood 
#' Limited to children in the complete-case analysis who fall in the "High" and "Very High" exposure groups. 
dlim_exp <- dpop %>%
  filter(flag_missing_covariates == FALSE) %>% 
  filter(ecat_prenat_4 %in% c("High", "Very High"))


#' Convert ecat_yr from reference/low/high to true (lived high) and false (not lived high): 
dlim_exp <- dlim_exp %>% 
  select(
    lopnr, fodelsear, 
    ecat_yr0:ecat_yr12, 
    ecat_prenat_4
  ) %>% 
  mutate(across(
    c(ecat_yr0:ecat_yr12),
    ~ .x == "high"
  ))

n_by_age_exp <- dlim_exp %>% 
  pivot_longer(
    cols = c(ecat_yr0:ecat_yr12),
    names_to = "child_age", 
    values_to = "flag_high_exposure"
  ) %>% 
  mutate(
    child_age = str_remove(child_age, "ecat_yr"),
    child_age = as.numeric(child_age)
  ) %>% 
  group_by(ecat_prenat_4, child_age) %>% 
  summarize(count_high_exposure = sum(flag_high_exposure)) 

n_by_age_exp <- dlim_exp %>% 
  group_by(ecat_prenat_4) %>% 
  count() %>% 
  left_join(n_by_age_exp, by = join_by(ecat_prenat_4)) %>% 
  arrange(child_age)

n_by_age_exp <- n_by_age_exp %>%  
  mutate(
    percent = count_high_exposure / n * 100
  )

# Exposure Verification ---------------------------------------------------
#' # Exposure Verification 
#' 
#' Summarize PFAS concentrations by address-based exposure categories 
#' and perform Jonckheere-Terpstra trend tests. 

dexp <- read_rds(here("data", "data_exp_validation_cohort.rds"))

dt_exp <- dexp %>%
  select(
    ban, 
    `PFHxS` = "pfhxs", 
    `PFOS` = "pfos", 
    `PFOA` = "pfoa",
    `Exposure Categories` = exp_4cat
  )

exp_summary <- CreateTableOne(
  vars = c("PFOS", "PFHxS", "PFOA"),
  data = dt_exp,
  test = FALSE,
  strata = "Exposure Categories") %>% 
  print(
    showAllLevels = TRUE,
    nonnormal = c("PFHxS", "PFOS", "PFOA"),
    quote = FALSE,
    noSpaces = TRUE,
    printToggle = FALSE,
    explain = FALSE,
    catDigits = 1,
    contDigits = 1
  ) %>%
  as_tibble(
    rownames = "Variable"
  ) %>% 
  select(-level) %>% 
  pivot_longer(
    cols = c("Low", "Intermediate", "High", "Very High"), 
    names_to = "Exposure Group", 
    values_to = "conc"
  ) %>% 
  pivot_wider(
    names_from = Variable, 
    values_from = conc
  ) %>% 
  arrange(desc(row_number()))

gt(exp_summary) %>% 
  tab_header(title = "PFAS Concentrations by Exposure Category",
             subtitle = "Women of childbearing age in the Ronneby Biomonitoring Cohort (N = 209)")

#' Save: 
write_csv2(exp_summary, file.path(output_loc, "t2-exposure-verification.csv"))


## Jonckheere-Terpstra test per PFAS ------------------------------------------------
#' ## Jonckheere-Terpstra test per PFAS

tdat <- dexp %>% 
  select(ban, exp_4cat, pfhxs, pfoa, pfos) %>% 
  mutate(exp_4cat = factor(exp_4cat, levels = c("Low", "Intermediate", "High", "Very High"))) %>% 
  pivot_longer(
    cols = c("pfhxs", "pfoa", "pfos"), 
    values_to = "conc", 
    names_to = "pfas"
  ) %>% 
  group_by(pfas) %>% 
  nest()

jt_tests <- map(
  tdat$data, ~ JonckheereTerpstraTest(.x$conc, .x$exp_4cat, alternative = "increasing", nperm = 100000)
)

names(jt_tests) <- tdat$pfas

jt_tests

## Plot Exposure Categories ------------------------------------------------
#' ## Plot Exposure Categories 
pdat <- dexp %>% 
  rename(exp_group = exp_4cat) %>% 
  mutate(exp_group = factor(exp_group, levels = c("Very High", "High", "Intermediate", "Low"))) %>% 
  pivot_longer(
    cols = c("pfhxs", "pfoa", "pfos"), 
    names_to = "pfas",
    values_to = "conc"
  ) %>% 
  mutate(pfas = factor(
    pfas, 
    levels = c("pfhxs", "pfos", "pfoa"),
    labels = c("PFHxS", "PFOS", "PFOA")))

comb_n <- pdat %>% 
  group_by(exp_group, pfas) %>% 
  summarize(
    y = max(conc, na.rm = T), 
    n = n(),
    .groups = "drop"
  )

ggplot(pdat, aes(x = exp_group, y = conc)) +
  facet_wrap(~ pfas, scales = "free", ncol = 1) + 
  geom_boxplot(outliers = FALSE) + 
  geom_quasirandom(color = "blue", alpha = 0.3) + 
  scale_y_continuous(expand = c(0.25,0)) +
  geom_text(data = comb_n, aes(x = exp_group, y = y, label = n), vjust = -0.5) + 
  labs(
    y = "Measured PFAS Concentration (ng/mL)",
    x = "Address-based exposure category"
  ) + 
  theme_bw() + 
  theme(text = element_text(family = "sans", size = 12))

ggsave(
  filename = "fs1-exp-verification.tiff", 
  plot = last_plot(), 
  device = "tiff", 
  path = output_loc,
  dpi = 600,
  width = 1877*2, 
  height = 2600*2, 
  units = "px"
)

# Disease Prevalence ------------------------------------------------------
#' # Disease Prevalence 

prev <- dlist %>%
  select(outcome_algo, data) %>%
  unnest(cols = c(data)) %>%
  filter(flag_missing_covariates == FALSE) %>%
  group_by(outcome_algo) %>%
  summarize(
    nsubj = n(),
    ncases = sum(flag_event),
    prevalence = round(ncases / nsubj * 100,1),
    .groups = "drop"
  ) %>%
  arrange(desc(ncases))

gt(prev)


## Outcome Combinations ----------------------------------------------------
#' ## Outcome Combinations 

dout_wide <- dlist %>%
  mutate(
    data = map(data, ~ filter(.x, flag_missing_covariates == FALSE & flag_event == TRUE)),
    lopnr = map(data, ~ .x %>%
                  select(lopnr) %>% 
                  mutate(flag = TRUE))
  ) %>% 
  select(algorithm, lopnr) %>% 
  unnest(cols = lopnr) %>% 
  pivot_wider(names_from = algorithm, values_from = flag) %>% 
  mutate(across(c("asth_ort", "wheeze", "asth_ort3"), ~ ifelse(is.na(.x), 0, 1))) %>% 
  rename(
    "Asthma" = "asth_ort", 
    "Asthma (3+)" = asth_ort3,
    "Wheeze" = "wheeze"
  )

# Create plot and save
tiff(file.path(output_loc, "fs2-outcome-combos.tiff"), 
     width = 2250, 
     height = 1800, 
     units = "px", 
     res = 300, 
     compression = "lzw")

UpSetR::upset(as.data.frame(dout_wide),
              sets = c("Asthma", "Asthma (3+)", "Wheeze"), 
              order.by = "freq",
              text.scale = 1.6)

dev.off()




# Paper Text --------------------------------------------------------------
#' # Paper Text 
#' 
#' This section contains inline R code to generate descriptive text for the manuscript. 
#' 


#' ## Methods 
#' We limited our validation cohort to women of childbearing age 
#' (`r min(dexp$age)` - `r max(dexp$age)` years). 



## Complete Cases ----------------------------------------------------------
#' ## Complete Cases 
#' 
#' Of the `r nfull` children born in Blekinge between 2006-2013, 
#' `r nlim` (`r round(nlim / nfull * 100 ,1)`%) 
#' had complete covariate information and were included in our final study population.


## Baseline Characteristics ------------------------------------------------
#' ## Baseline Characteristics
#' 
#' Most children in the study had older siblings 
#' (`r round(sum(dt_lim$Parity == "Multiparous") / nrow(dt_lim) * 100, 0)`%), 
#' 
#' were born to non-smoking mothers 
#' (`r round(sum(dt_lim$"Maternal smoking status" == "Non-smoker")/ nrow(dt_lim) * 100, 0)`%)
#' 
#' and had parents who were born in Sweden  
#' (`r round(sum(dt_lim$"Parent born abroad" == FALSE)/ nrow(dt_lim) * 100, 0)`%). 
#' 
#' The median maternal age at delivery was 
#' `r median(dt_lim$"Maternal age at delivery")` years (
#' IQR: 
#' `r round(quantile(dt_lim$"Maternal age at delivery", c(0.25, 0.75)), 0)`). 
#' 
#' Overall, 
#' `r round(sum(dt_lim$"Parental asthma" == TRUE)/ nrow(dt_lim) * 100, 0)`% of children in the study 
#' had at least one parent with asthma and nearly one-quarter 
#' (`r round(nrow(dt_lim %>% group_by(lopnr_mor) %>% count() %>% filter(n > 1)) / nrow(dt_lim) * 100, 0)`%) 
#' had at least one maternal sibling also included in the study population. 
#' 
#' The study included 
#' `r n_exp[n_exp$"Prenatal exposure group" == "Very High", "n"][[1]]` children 
#' (`r n_exp[n_exp$"Prenatal exposure group" == "Very High", "percent"][[1]]`%) 
#' who were classified in the very high prenatal exposure group,
#' 
#' `r n_exp[n_exp$"Prenatal exposure group" == "High", "n"][[1]]` children 
#' (`r n_exp[n_exp$"Prenatal exposure group" == "High", "percent"][[1]]`%) 
#' who were classified in the high prenatal exposure group,
#' 
#' `r n_exp[n_exp$"Prenatal exposure group" == "Intermediate", "n"][[1]]` children 
#' (`r n_exp[n_exp$"Prenatal exposure group" == "Intermediate", "percent"][[1]]`%) 
#' in the intermediate prenatal exposure group, and
#' 
#' `r n_exp[n_exp$"Prenatal exposure group" == "Low", "n"][[1]]` children 
#' (`r n_exp[n_exp$"Prenatal exposure group" == "Low", "percent"][[1]]`%) in the low exposure group.
#' 

## Outcomes ----------------------------------------------------------------
#' ## Outcomes 

dout_wide_all <- filter(dout_wide, Asthma == 1 & Wheeze == 1 & `Asthma (3+)` == 1)
dout_wide_wheeze <- filter(dout_wide, Asthma == 0 & Wheeze == 1 & `Asthma (3+)` == 0)


#' Of the three outcomes evaluated in this study, wheeze had the highest
#' overall prevalence in the study population (`r prev$prevalence[[1]]`%).
#' 
#' The prevalence of asthma was `r prev$prevalence[[2]]`% 
#' and the prevalence of asthma (3+) was `r prev$prevalence[[3]]`%.
#' 
#' 
#' 
#' A total of 
#' `r nrow(dout_wide)` study subjects 
#' (`r round(nrow(dout_wide) / prev$nsubj[[1]] * 100, 1)`%)
#' had at least one outcome.  
#' 
#' 
#' For individuals with at least one outcome, it was most common to have all all three outcomes  
#' `r nrow(dout_wide_all)`; 
#' (`r round(nrow(dout_wide_all) / nrow(dout_wide) * 100, 1)`) 
#' or to just have wheeze: 
#' `r nrow(dout_wide_wheeze)`
#' (`r round(nrow(dout_wide_wheeze) / nrow(dout_wide) * 100, 1)`)


## Exposure over Childhood -------------------------------------------------
#' ## Exposure over Childhood

#' At age 3, most children in the very high prenatal exposure group 
#' (`r round(n_by_age_exp[n_by_age_exp$ecat_prenat_4 == "Very High" & n_by_age_exp$child_age == 3,]$percent, 1)`%) 
#' still lived at a high-exposed address. 





