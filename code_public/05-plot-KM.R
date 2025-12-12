#' ---
#' title: "05 - Plot primary KM models"
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
#' This script generates Kaplan-Meier plots (cumulative incidence curves) 
#' for the three main outcomes by prenatal exposure group, 
#' and saves the combined figure used as Figure 1 in the manuscript. 
#' 
#' ## Inputs 
#' 
#' - `cohort-outcome-data.rds` (main cohort with outcome-specific follow-up time)
#' 
#' ## Outputs
#' - `f1-km-plots-primary.tiff`
#' - `f1-km-plots-primary.eps`


# Setup -------------------------------------------------------------------
#' # Setup 

#+ load-libraries, message = F
library(here)
library(tidyverse)
library(survival)
library(survminer)
library(cowplot)

output_loc <- here("outputs_public")

sessionInfo()

# Load Data ---------------------------------------------------------------
#' # Load Data 
dlist <- read_rds(file.path(output_loc, "cohort-outcome-data.rds")) %>% 
  mutate(data = map(data, ~ filter(.x, flag_missing_covariates == FALSE)))

# KM Plots ----------------------------------------------------------------
#' # KM Plots
palette <- RColorBrewer::brewer.pal(name = "Dark2", n = 4)
rev_palette <- rev(palette)

theme_font <- theme(legend.text = element_text(family = "sans", size = 8),
                    legend.title = element_text(family = "sans", size = 8),
                    axis.title.x = element_text(family = "sans", size = 10),
                    axis.title.y = element_text(family = "sans", size = 10),
                    axis.text.x = element_text(family = "sans", size = 8),
                    axis.text.y = element_text(family = "sans", size = 8),
                    plot.title = element_text(family = "sans", size = 12),
                    plot.subtitle = element_text(family = "sans", size = 8),
                    axis.line.x.bottom = element_line(color = "black", linewidth = 0.3),
                    axis.line.y.left = element_line(color = "black", linewidth = 0.3),
                    panel.grid.minor = element_blank(), 
                    panel.grid.major.y = element_line(color = "grey", linewidth = 0.3),
                    axis.ticks = element_line(color = "black", linewidth = 0.3))

# Plot cumulative incidence with exposure ordered Background -> Very high 
PlotCI_survplot <- function(km_dat, outcome_label) {
  
  km_fit <- survfit(Surv(time = age_final, event = flag_event) ~ ecat_prenat_4, data = km_dat)
  # ggsurvplot(km_fit, data = km_dat)
  
  if(outcome_label != "Wheeze") {
    p1 <- ggsurvplot(fit = km_fit, 
                     data = km_dat, 
                     fun = function(s) 1-s,
                     palette = palette,
                     break.x.by =  4,
                     break.y.by = 0.2,
                     xlab = "Time to diagnosis (years)",
                     ylab = "Cumulative Incidence",
                     title = outcome_label,
                     legend = c(0.1, 0.85),
                     censor = FALSE,
                     size = 0.3,
                     legend.title = "Prenatal exposure group", 
                     legend.labs =c("Background", "Intermediate", "High", "Very High"),
                     risk.table = TRUE,
                     pval = FALSE,
                     pval.size = 5, 
                     pval.method = FALSE,
                     pval.method.size = 5,
                     fontsize = 2.8,
                     tables.theme = theme_survminer(font.tickslab = 8, font.main = 10))  
  } else {
    p1 <- ggsurvplot(fit = km_fit, 
                     data = km_dat, 
                     fun = function(s) 1-s,
                     palette = palette,
                     break.x.by =  1,
                     break.y.by = 0.2,
                     xlim = c(0,3),
                     xlab = "Time to diagnosis (years)",
                     ylab = "Cumulative Incidence",
                     title = outcome_label,
                     legend = c(0.1, 0.85),
                     censor = FALSE,
                     size = 0.3,
                     legend.title = "Prenatal exposure group", 
                     legend.labs =c("Background", "Intermediate", "High", "Very High"),
                     risk.table = TRUE,
                     pval = FALSE,
                     pval.size = 5, 
                     pval.method = FALSE,
                     pval.method.size = 5,
                     fontsize = 2.8,
                     tables.theme = theme_survminer(font.tickslab = 8, font.main = 10))
  }
  
  
  return(p1)
}

# Plot cumulative incidence with exposure reversed (Very high -> Background)
# used for number-at-risk tables so Very High is the top row
PlotCI_survplot_rev <- function(km_dat, outcome_label) {
  
  km_dat <- mutate(km_dat, ecat_prenat_4 = fct_rev(ecat_prenat_4))
  
  km_fit <- survfit(Surv(time = age_final, event = flag_event) ~ ecat_prenat_4, data = km_dat)
  # ggsurvplot(km_fit, data = km_dat)
  
  if(outcome_label != "Wheeze") {
    p1 <- ggsurvplot(fit = km_fit, 
                     data = km_dat, 
                     fun = function(s) 1-s,
                     palette = rev_palette,
                     break.x.by =  4,
                     break.y.by = 0.2,
                     xlab = "Time to diagnosis (years)",
                     ylab = "Cumulative Incidence",
                     title = outcome_label,
                     legend = c(0.1, 0.85),
                     censor = FALSE,
                     size = 0.3,
                     legend.title = "Prenatal exposure group", 
                     legend.labs = c("Very High", "High", "Intermediate", "Background"),
                     risk.table = TRUE,
                     fontsize = 2.8,
                     tables.theme = theme_survminer(font.tickslab = 8, font.main = 10))  
  } else {
    p1 <- ggsurvplot(fit = km_fit, 
                     data = km_dat, 
                     fun = function(s) 1-s,
                     palette = rev_palette,
                     break.x.by =  1,
                     break.y.by = 0.2,
                     xlim = c(0,3),
                     xlab = "Time to diagnosis (years)",
                     ylab = "Cumulative Incidence",
                     title = outcome_label,
                     legend = c(0.1, 0.85),
                     censor = FALSE,
                     size = 0.3,
                     legend.title = "Prenatal exposure group", 
                     legend.labs = c("Very High", "High", "Intermediate", "Background"),
                     risk.table = TRUE,
                     fontsize = 2.8,
                     tables.theme = theme_survminer(font.tickslab = 8, font.main = 10))
  }
  
  
  return(p1)
}


# Basic plots: 
ci_list_survplot <-     map2(dlist$data, dlist$algorithm_label, PlotCI_survplot)
ci_list_survplot_rev <- map2(dlist$data, dlist$algorithm_label, PlotCI_survplot_rev)

names(ci_list_survplot)     <- dlist$algorithm_label
names(ci_list_survplot_rev) <- dlist$algorithm_label

# Edit plots: 
plot_asthma <- ci_list_survplot$Asthma$plot +
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_x_continuous(limits = c(0, 13), 
                     breaks = c(0, 4, 8, 12), 
                     expand = expansion(mult = c(0.1, 0.1))) + 
  theme(panel.grid.major.y = element_line(color = "grey")) + 
  theme_font

plot_asthma3 <- ci_list_survplot$`Asthma (3+)`$plot +
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_x_continuous(limits = c(0, 13), 
                     breaks = c(0, 4, 8, 12), 
                     expand = expansion(mult = c(0.1, 0.1))) + 
  theme(panel.grid.major.y = element_line(color = "grey")) + 
  theme_font

plot_wheeze <- ci_list_survplot$Wheeze$plot +
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_x_continuous(limits = c(0, 4), 
                     breaks = c(0, 1, 2, 3),
                     expand = expansion(mult = c(0.1, 0.1))) + 
  theme(panel.grid.major.y = element_line(color = "grey")) + 
  theme_font

# Edit tables: 
table_asthma <- ci_list_survplot_rev$Asthma$table + 
  scale_x_continuous(limits = c(0, 13), 
                     breaks = c(0, 4, 8, 12), 
                     expand = expansion(mult = c(0.1, 0.1))) + 
  labs(title = "Number at risk") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "grey"),
        axis.line.x.top = element_line(color = "black"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.3),
        axis.line.y.left = element_line(color = "black", linewidth = 0.3),
        axis.line.y.right = element_blank())

table_asthma3 <- ci_list_survplot_rev$`Asthma (3+)`$table + 
  scale_x_continuous(limits = c(0, 13), 
                     breaks = c(0, 4, 8, 12), 
                     expand = expansion(mult = c(0.1, 0.1))) + 
  labs(title = "Number at risk") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "grey"),
        axis.line.x.top = element_line(color = "black"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.3),
        axis.line.y.left = element_line(color = "black", linewidth = 0.3),
        axis.line.y.right = element_blank()) 

table_wheeze <- ci_list_survplot_rev$Wheeze$table + 
  scale_x_continuous(limits = c(0, 4), 
                     breaks = c(0, 1, 2, 3),
                     expand = expansion(mult = c(0.1, 0.1))) + 
  labs(title = "Number at risk") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),        
        axis.text.x = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "grey"),
        axis.line.x.top = element_line(color = "black"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.3),
        axis.line.y.left = element_line(color = "black", linewidth = 0.3),
        axis.line.y.right = element_blank()) 

cp_wheeze <- cowplot::plot_grid(
  plot_wheeze + theme(legend.position = "none"), 
  table_wheeze,
  nrow = 2, align = 'vh', axis = 'l',
  rel_heights = c(1, 0.5))

cp_asthma <- cowplot::plot_grid(
  plot_asthma + 
    theme(legend.position = "none", axis.title.y = element_blank()), 
  table_asthma,
  nrow = 2, align = 'vh', axis = 'l',
  rel_heights = c(1, 0.5))

cp_asthma3 <- cowplot::plot_grid(
  plot_asthma3 + 
    theme(legend.position = "none", axis.title.y = element_blank()), 
  table_asthma3,
  nrow = 2, align = 'vh', axis = 'l',
  rel_heights = c(1, 0.5))

cp_legend <- cowplot::get_plot_component(plot_asthma +
                                           theme(legend.position = "bottom",
                                                 legend.text = element_text(size = 12),
                                                 legend.title = element_text(size = 12),
                                                 legend.box.margin = margin(0,0,0,0.1)),
                                         "guide-box-bottom")

# Figure: Three columns ------------------------------------------------
#' # Figure: Three columns 
#' 
#' Maximum width: 2250 pixesl at 300 dpi
cp_comb <- cowplot::plot_grid(cp_wheeze, cp_asthma, cp_asthma3, ncol = 3, rel_widths = c(1.1, 1, 1), align = 'vh')

# Plot 
cowplot::plot_grid(cp_comb, cp_legend, ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "f1-km-plots-primary.tiff", 
       plot = last_plot(), 
       path = output_loc,
       dpi = 300,
       device = "tiff", 
       width = 7.5, 
       height = 5,
       units = "in")

ggsave(filename = "f1-km-plots-primary.eps", 
       plot = last_plot(), 
       path = output_loc,
       dpi = 300,
       device = cairo_ps, 
       width = 7.5, 
       height = 5,
       family = "Arial", 
       units = "in")




