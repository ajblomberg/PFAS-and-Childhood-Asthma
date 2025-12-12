#' ---
#' title: "10 - Rubin Causal Model"
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
#' This script implements the Rubin Causal Model analysis for the matched 
#' comparison between the "Very High" and "Background" prenatal exposure groups. 
#' 
#' It performs: 
#' 
#' - Kaplan-Meier estimation of cumulative incidence in the matched sample 
#' - Calculation of the observed test statistic (difference in cumulative incidence)
#' - Fisher randomization test using 100,000 reassignments of treatment within matched sets 
#' 
#' ## Inputs 
#' 
#' - `dvh_matched_out_list.rds`: Outcome-specific datasets for the matched cohort 
#' 
#' ## Outputs
#' 
#' - `f2-km-plots-rcm.tiff`
#' - `fs8-fisher-test-statistics-distributions.tiff`
#' 

# Setup -------------------------------------------------------------------
#' # Setup
#' 
#+ load-libraries, message = F
library(here)
library(tidyverse)
library(gt)
library(survival)
library(survminer)
library(cowplot)

output_loc <- here("outputs_public")

set.seed(314)

dlist <- read_rds(file.path(output_loc, "dvh_matched_out_list.rds"))


# KM Plots ----------------------------------------------------------------
#' # KM Plots
palette <- RColorBrewer::brewer.pal(name = "Dark2", n = 4)[c(1,4)]
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
                     xlab = "Time to diagnosis (age)",
                     ylab = "Cumulative Incidence",
                     title = outcome_label,
                     legend = c(0.1, 0.85),
                     censor = FALSE,
                     size = 0.3,
                     legend.title = "Prenatal exposure group", 
                     legend.labs =c("Background", "Very High"),
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
                     xlab = "Time to diagnosis (age)",
                     ylab = "Cumulative Incidence",
                     title = outcome_label,
                     legend = c(0.1, 0.85),
                     censor = FALSE,
                     size = 0.3,
                     legend.title = "Prenatal exposure group", 
                     legend.labs =c("Background", "Very High"),
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
                     xlab = "Time to diagnosis (age)",
                     ylab = "Cumulative Incidence",
                     title = outcome_label,
                     legend = c(0.1, 0.85),
                     censor = FALSE,
                     size = 0.3,
                     legend.title = "Prenatal exposure group", 
                     legend.labs = c("Very High", "Background"),
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
                     xlab = "Time to diagnosis (age)",
                     ylab = "Cumulative Incidence",
                     title = outcome_label,
                     legend = c(0.1, 0.85),
                     censor = FALSE,
                     size = 0.3,
                     legend.title = "Prenatal exposure group", 
                     legend.labs = c("Very High", "Background"),
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
  rel_heights = c(1, 0.35))

cp_asthma <- cowplot::plot_grid(
  plot_asthma + theme(legend.position = "none",
                      axis.title.y = element_blank()),
  table_asthma,
  nrow = 2, align = 'vh', axis = 'l',
  rel_heights = c(1, 0.35))

cp_asthma3 <- cowplot::plot_grid(
  plot_asthma3 + theme(legend.position = "none",
                       axis.title.y = element_blank()),
  table_asthma3,
  nrow = 2, align = 'vh', axis = 'l',
  rel_heights = c(1, 0.35))

cp_legend <- cowplot::get_plot_component(plot_asthma +
                                           theme(legend.position = "bottom",
                                                 legend.text = element_text(size = 12),
                                                 legend.title = element_text(size = 12),
                                                 legend.box.margin = margin(0,0,0,0.1)),
                                         "guide-box-bottom")

cp_comb <- cowplot::plot_grid(cp_wheeze, cp_asthma, cp_asthma3, ncol = 3, rel_widths = c(1.1, 1, 1), align = 'vh')

# Plot 
cowplot::plot_grid(cp_comb, cp_legend, ncol = 1, rel_heights = c(1, 0.05))

# Maximum width: 2250 pixels at 300 dpi
ggsave(filename = "f2-km-plots-rcm.tiff",
       plot = last_plot(),
       device = "tiff",
       path = output_loc,
       dpi = 600,
       width = 2250*2,
       height = 1360*2,
       units = "px")


# Cumulative Incidence  -------------------------------------------------------
#' # Cumulative Incidence 

list_km_est <- dlist %>% 
  mutate(km_fit = map(data, ~ survfit(Surv(time = age_final, event = flag_event) ~ ecat_prenat_4 + cluster(subclass), data = .x)),
         km_sum = map2(km_fit, age_max + 1, ~ summary(.x, times = .y))) %>% 
  mutate(km_table_cum = map(km_sum, ~ tibble(exp = .x$strata, 
                                             probs = 1 - .x$surv,
                                             ci_lo = 1 - .x$lower,
                                             ci_hi = 1- .x$upper))) 

t1 <- list_km_est %>% 
  select(algorithm, km_table_cum) %>%
  unnest(cols = c(km_table_cum)) %>% 
  select(-ci_lo, -ci_hi) %>% 
  mutate(exp = str_sub(exp, 15, -1)) %>% 
  pivot_wider(names_from = exp, values_from = probs) %>% 
  mutate(ci_diff = `Very High` - Low)

gt(t1)


# Fisher Exact p-values ---------------------------------------------------
#' # Fisher Exact p-values 

list_data <- map(dlist$data, ~ select(.x, age_final, flag_event, ecat_prenat_4, subclass) %>% 
                   arrange(subclass))

names(list_data) <- dlist$algorithm_label

list_age_final  <- map(list_data, ~ .x$age_final)
list_flag_event <- map(list_data, ~ .x$flag_event)
list_wobs       <- map(list_data, ~ (as.numeric(.x$ecat_prenat_4 == "Very High")))


list_event_times <- map2(list_age_final, list_flag_event, ~ sort(unique(.x[.y == 1])))

list_K <- map(list_event_times, ~ length(.x))
list_n <- map(list_wobs, ~ length(.x))

list_Mrisk  <- map2(list_K, list_n, ~ matrix(0L, nrow = .x, ncol = .y))
list_Mevent <- map2(list_K, list_n, ~ matrix(0L, nrow = .x, ncol = .y))


#' Create matrices for all three datasets:
for(alg in names(list_n)) {
  n_i <- list_n[[alg]]
  K_i <- list_K[[alg]]
  
  for (i in seq_len(n_i)) { 
    for (k in seq_len(K_i)) {
      
      list_Mrisk[[alg]][k, i] <-  
        list_age_final[[alg]][i] >= list_event_times[[alg]][k]
      
      list_Mevent[[alg]][k,i] <- 
        list_flag_event[[alg]][i] == 1 && abs(list_age_final[[alg]][i] - list_event_times[[alg]][k]) <1e-15
    }
  }
}



# for (i in seq_len(list_n[[1]])) { 
#   for (k in seq_len(list_K[[1]])) { 
#     list_Mrisk[[1]][k, i] <- list_age_final[[1]][i] >= list_event_times[[1]][k]
#     list_Mevent[[1]][k, i] <- list_flag_event[[1]][i] == 1 && abs(list_age_final[[1]][i] - list_event_times[[1]][k]) <1e-15
#   }
# }
# 
# for (i in seq_len(list_n[[2]])) { 
#   for (k in seq_len(list_K[[2]])) { 
#     list_Mrisk[[2]][k, i] <- list_age_final[[2]][i] >= list_event_times[[2]][k]
#     list_Mevent[[2]][k, i] <- list_flag_event[[2]][i] == 1 && abs(list_age_final[[2]][i] - list_event_times[[2]][k]) <1e-15
#   }
# }
# 
# for (i in seq_len(list_n[[3]])) { 
#   for (k in seq_len(list_K[[3]])) { 
#     list_Mrisk[[3]][k, i] <- list_age_final[[3]][i] >= list_event_times[[3]][k]
#     list_Mevent[[3]][k, i] <- list_flag_event[[3]][i] == 1 && abs(list_age_final[[3]][i] - list_event_times[[3]][k]) <1e-15
#   }
# }

#' Function for cumulative incidence 
km_diff <- function(W, M_risk, M_event, K) { 
  # For each distinct time k up to 13, get # at risk and # events
  
  # group = 1
  Y1_k <- M_risk  %*% W
  d1_k <- M_event %*% W 
  
  # group = 0
  W0 <- 1-W
  Y0_k <- M_risk %*% W0
  d0_k <- M_event %*% W0
  
  # Convert to vectors
  Y1_k <- as.numeric(Y1_k)
  d1_k <- as.numeric(d1_k)
  Y0_k <- as.numeric(Y0_k)
  d0_k <- as.numeric(d0_k)
  
  # KM for group 1 at time <= 13
  S1 <- 1
  S0 <- 1
  
  for(k in seq_len(K)) {
    if(Y1_k[k] > 0) S1 <- S1 * (1 - d1_k[k]/Y1_k[k])
    if(Y0_k[k] > 0) S0 <- S0 * (1 - d0_k[k]/Y0_k[k])
  }
  return(S0 - S1) 
}



#' Manual Tobs matches Tobs from the KM model.  
Tobs <- vector("list", length(list_n))
for (i in seq_along(list_n)) {
  Tobs[i] <- km_diff(list_wobs[[i]], list_Mrisk[[i]], list_Mevent[[i]], list_K[[i]])
}

Tobs


## Sampling Structure  -----------------------------------------------------
#' ## Sampling Structure 
# Dimensions
nsample <- 100000
nrow <- nrow(dlist$data[[1]]) # same for all three data tables
nsub <- nrow / 4

# Create Matrix
all_samples <- matrix(rep(c(1,0,0,0), times = nsub), nrow = 4, ncol = nsub * nsample)
Wmatrix <- apply(all_samples, 2, sample)
Wmatrix <- matrix(as.vector(Wmatrix), 
                  nrow = nrow, 
                  ncol = nsample)


## Cumulative Incidence from Sampled Values -----------------------------
#' ## Cumulative Incidence from Sampled Values 

Trep <- vector("list", length(list_n))

for(i in seq_along(list_n)) {

  Trep[[i]] <- matrix(ncol = 1, nrow = nsample)

  for(j in 1:nsample) {
    Trep[[i]][j] <- km_diff(Wmatrix[,j], list_Mrisk[[i]], list_Mevent[[i]], list_K[[i]])
  }
}

names(Trep) <- dlist$algorithm_label


#' Create basic plots 
plots_trep <- vector("list", length(list_n))
names(plots_trep) <- dlist$algorithm_label

for (i in seq_along(list_n)) {
  plots_trep[[i]] <- ggplot(tibble(Trep = Trep[[i]]), aes(x = Trep)) + 
    geom_histogram() +
    geom_vline(xintercept = Tobs[[i]], color = "red") + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) + 
    labs(title = dlist$algorithm_label[[i]], y = "Frequency", x = "Test Statistic") + 
    theme_bw() + 
    theme_font
  
}


#' Calculate p-value: 
list_pval <- vector("list", length(list_n))
names(list_pval) <- dlist$algorithm_label

for(i in seq_along(list_n)){ 
  list_pval[[i]] <- sum(Trep[[i]] >= Tobs[[i]]) / nsample
}


#' Turn list pval into numeric vector: 
pvals <- unlist(list_pval)

FormatPvals <- function(p) {
  if (p < 0.001) "Fisherian p-value < 0.001" 
  else paste0("Fisherian p-value =", signif(p,1))
}

subtitles <- tibble(algorithm_label = names(pvals),
                    subtitle = vapply(pvals,FormatPvals, character(1)))

GetSubtitle <- function(name) {
  subtitles$subtitle[subtitles$algorithm_label == name]
}

#' Plot with subtitles: note that the p-values are manual here
p_wheeze <- plots_trep[["Wheeze"]] + 
  theme_font + 
  scale_x_continuous(breaks = c(-0.12, -0.06, 0, 0.06, 0.12)) + 
  labs(subtitle = GetSubtitle("Wheeze"))

p_asthma <- plots_trep[["Asthma"]] + 
  theme_font + 
  theme(axis.title.y = element_blank()) + 
  scale_x_continuous(breaks = c(-0.12, -0.06, 0, 0.06, 0.12)) + 
  labs(subtitle = GetSubtitle("Asthma"))

p_asthma3 <- plots_trep[["Asthma (3+)"]] + 
  theme_font + 
  theme(axis.title.y = element_blank()) +
  labs(subtitle = GetSubtitle("Asthma (3+)"))

cowplot::plot_grid(p_wheeze, p_asthma, p_asthma3, 
                   ncol = 3,
                   rel_widths = c(1, 0.96, 0.96))

ggsave(filename = "fs8-fisher-test-statistics-distributions.tiff",
       plot = last_plot(),
       device = "tiff",
       path = output_loc,
       dpi = 600,
       width = 1877*2,
       height = 900*2,
       units = "px")

ggsave(filename = "fs8-fisher-test-statistics-distributions.png", 
       plot = last_plot(), 
       device = "png", 
       path = output_loc,
       dpi = 600,
       width = 1877*2,
       height = 900*2,
       units = "px")


# Add Fisherian p-values to KM plots --------------------------------------
#' # Add Fisherian p-values to KM plots 

cp_wheeze <- cowplot::plot_grid(
  plot_wheeze + 
    theme(legend.position = "none") + 
    labs(subtitle = GetSubtitle("Wheeze")),
  table_wheeze,
  nrow = 2, align = 'vh', axis = 'l',
  rel_heights = c(1, 0.5))

cp_asthma <- cowplot::plot_grid(
  plot_asthma + 
    theme(legend.position = "none", axis.title.y = element_blank()) + 
    labs(subtitle = GetSubtitle("Asthma")),
  table_asthma,
  nrow = 2, align = 'vh', axis = 'l',
  rel_heights = c(1, 0.5))

cp_asthma3 <- cowplot::plot_grid(
  plot_asthma3 + 
    theme(legend.position = "none", axis.title.y = element_blank()) + 
    labs(subtitle = GetSubtitle("Asthma (3+)")),
  table_asthma3,
  nrow = 2, align = 'vh', axis = 'l',
  rel_heights = c(1, 0.5))

cp_legend <- cowplot::get_plot_component(plot_asthma +
                                           theme(legend.position = "bottom",
                                                 legend.text = element_text(size = 12),
                                                 legend.title = element_text(size = 12),
                                                 legend.box.margin = margin(0,0,0,0.1)),
                                         "guide-box-bottom")

#' Maximum width: 2250 pixesl at 300 dpi
cp_comb <- cowplot::plot_grid(cp_wheeze, cp_asthma, cp_asthma3, ncol = 3, rel_widths = c(1.1, 1, 1), align = 'vh')
cowplot::plot_grid(cp_comb, cp_legend, ncol = 1, rel_heights = c(1, 0.05))

#' Save 
ggsave(filename = "f2-km-plots-rcm.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       path = output_loc,
       dpi = 600,
       width = 7.5, 
       height = 5,
       units = "in")

ggsave(filename = "f2-km-plots-rcm.eps", 
       plot = last_plot(), 
       path = output_loc,
       dpi = 300,
       device = cairo_ps, 
       width = 7.5, 
       height = 5,
       family = "Arial", 
       units = "in")
