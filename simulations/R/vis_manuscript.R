library(here)
library(biotmle)
library(data.table)
library(tidyverse)
library(latex2exp)
library(patchwork)
library(ggsci)
library(wesanderson)
library(conflicted)
conflict_prefer("filter", "dplyr")
pd_error <- position_dodge(0.5)
pd_power <- position_dodge(0.9)
pal <- wes_palette("Zissou1")

##############################################################################
# Utility functions and useful constants
##############################################################################

# nominal control rates
fwer_cntrl <- c(0.01, 0.05, 0.1)
fdr_cntrl <- c(0.01, 0.05, 0.1)
fdr_cutoff <- 0.05

# simulation parameters across various settings
nonnull_genes <- 150 * c(0.1, 0.3, 0.5)
results_label <- c("min", "mod", "max")
true_psi <- 0 # A has no effect on Y => ATE = 0

# function to create useful summary measures for computing FDP, NPV, etc.
make_stat_table <- function(sim_data, fdr_cntrl_level) {
  # useful quantities for summary measures
  stat_table <- sim_data %>%
    group_by(sim_id, n_samp, method) %>%
    summarize(
      num_neg = sum(label == "null", na.rm = TRUE),
      num_pos = sum(label != "null", na.rm = TRUE),
      num_neg_fdr = sum(p_adj_fdr > fdr_cntrl_level, na.rm = TRUE),
      num_pos_fdr = sum(p_adj_fdr < fdr_cntrl_level, na.rm = TRUE),
      num_neg_fwer = sum(p_adj_fwer > fdr_cntrl_level, na.rm = TRUE),
      num_pos_fwer = sum(p_adj_fwer < fdr_cntrl_level, na.rm = TRUE)
    )
  # spit out useful table
  return(stat_table)
}

# make FDR table by properly computing rate across simulations
make_summary_table <- function(sim_data, fdr_cntrl) {
  # loop over FDR cutoffs
  table_fdr <- lapply(fdr_cntrl, function(fdr_cntrl_level) {
    # compute useful quantities for FDP, NPV, etc.
    metric_summary <- sim_data %>%
      filter(!is.na(psi)) %>%
      make_stat_table(fdr_cntrl_level = fdr_cntrl_level)

    # combine with metric summary with simulation results
    sim_data_metrics <- merge(
      sim_data, metric_summary,
      by = c("sim_id", "n_samp", "method")
    )

    # make FDR table
    table_fdr <- sim_data_metrics %>%
      filter(!is.na(psi)) %>%
      group_by(sim_id, n_samp, method) %>%
      mutate(
        fdp = if_else(num_pos_fdr != 0,
            sum(p_adj_fdr < fdr_cntrl_level & label == "null", na.rm = TRUE) /
              num_pos_fdr, 0
          ),
        npv = if_else(num_neg_fdr != 0,
            sum(p_adj_fdr > fdr_cntrl_level & label == "null", na.rm = TRUE) /
              num_neg_fdr, 1
          ),
        tpr = sum(p_adj_fdr < fdr_cntrl_level & label != "null",
                  na.rm = TRUE) / num_pos,
        tnr = sum(p_adj_fdr > fdr_cntrl_level & label == "null",
                  na.rm = TRUE) / num_neg
      )

    # return FDR table
    return(table_fdr)
  })

  # create summary table
  summary_fdr <- lapply(table_fdr, function(this_fdr_table) {
    # summarize by marginalizing over simulation iterations
    fdr_summary <- this_fdr_table %>%
      ungroup(sim_id) %>%
      #group_by(n_samp, method) %>%
      summarize(
        fdr = mean(fdp),
        fdp_worst = max(fdp),
        tnr_worst = min(tnr),
        npv_worst = min(npv),
        tpr_worst = min(tpr)
      )
    return(fdr_summary)
  })

  # NOTE: older implementation, manually checked against
  # make FDR table
  #table_fdr <- sim_data %>%
    #filter(!is.na(psi)) %>%
    #mutate(
      #reject_low = (p_adj_fdr < fdr_cntrl[1]),
      #reject_med = (p_adj_fdr < fdr_cntrl[2]),
      #reject_high = (p_adj_fdr < fdr_cntrl[3])
    #) %>%
    #mutate(
      #fdp_low = ifelse(sum(reject_low) > 0,
                       #sum(label == "null" & reject_low) /
                           #max(1, sum(reject_low)), 0),
      #fdp_med = ifelse(sum(reject_med) > 0,
                       #sum(label == "null" & reject_med) /
                           #max(1, sum(reject_med)), 0),
      #fdp_high = ifelse(sum(reject_high) > 0,
                        #sum(label == "null" & reject_high) /
                            #max(1, sum(reject_high)), 0)
    #) %>%
    #ungroup(sim_id) %>%
    #group_by(method, n_samp) %>%
    #summarise(
      #fdr_low = mean(fdp_low),
      #fdr_med = mean(fdp_med),
      #fdr_high = mean(fdp_high)
    #)
  return(list(table = table_fdr, summary = summary_fdr))
}

# convenience function for worst-case power
make_power_table <- function(sim_data, err_type = c("fdr", "fwer"), p_cntrl) {
  table_power <- sim_data %>%
    filter(!is.na(psi), label == "hit") %>%
    group_by(n_samp, method, gene_id) %>%
    summarise(
      reject_err_low = mean(get(paste0("p_adj_", err_type)) < p_cntrl[1]),
      reject_err_med = mean(get(paste0("p_adj_", err_type)) < p_cntrl[2]),
      reject_err_high = mean(get(paste0("p_adj_", err_type)) < p_cntrl[3]),
    ) %>%
    summarize(
      power_err_low = min(reject_err_low),
      power_err_med = min(reject_err_low),
      power_err_high = min(reject_err_low)
    )
  return(table_power)
}

##############################################################################
# SIMULATION #1 - global null
##############################################################################

# results files and settings
for (sim_dir in c("g0prob_no", "g0prob_yes")) {
  if (sim_dir == "g0prob_no") {
    sim_files <- paste0("sl-global_null-pval_logistic-std_eif-",
                       c("nsim_300-cvfolds_1-2020-04-15_08:56:56.rds",
                         "nsim_300-cvfolds_2-2020-04-20_18:33:13.rds"))
  } else if (sim_dir == "g0prob_yes") {
    sim_files <- paste0("sl-global_null-pval_logistic-std_eif-",
                       c("nsim_300-cvfolds_1-2020-04-15_20:21:05.rds",
                         "nsim_300-cvfolds_2-2020-04-16_20:31:34.rds"))
  }

  # read in results
  sim_results <- readRDS(
    here("data", paste0("main_", sim_dir), sim_files[1])
  )
  sim_results_cv <- readRDS(
    here("data", paste0("main_", sim_dir), sim_files[2])
  )

  # properly compute adjusted p-values across methods and sample sizes
  sim_data <- sim_results %>%
    bind_rows(.id = "n_samp") %>%
    mutate(
      n_samp = as.numeric(str_remove(n_samp, "n_")),
    ) %>%
    group_by(method, n_samp, sim_id) %>%
    mutate(
      p_adj_fwer = p.adjust(p_val, method = "bonferroni"),
      p_adj_fdr = p.adjust(p_val, method = "BH"),
      label = rep("null", n())
    )
  no_sim <- length(unique(sim_data$sim_id))

  sim_data_cv <- sim_results_cv %>%
    bind_rows(.id = "n_samp") %>%
    mutate(
      n_samp = as.numeric(str_remove(n_samp, "n_")),
    ) %>%
    group_by(method, n_samp, sim_id) %>%
    mutate(
      p_adj_fwer = p.adjust(p_val, method = "bonferroni"),
      p_adj_fdr = p.adjust(p_val, method = "BH"),
      label = rep("null", n())
    ) %>%
    ungroup() %>%
    mutate(
      method = paste(method, "cv", sep = "_")
    ) %>%
    group_by(method, n_samp, sim_id)

  # SUMMARIZE FDR - global null
  table_summary_fdr <- make_summary_table(sim_data, fdr_cntrl)
  table_summary_fdr_cv <- make_summary_table(sim_data_cv, fdr_cntrl)

  # combine tables from non-sample-split and CV results
  table_fdr_all <- lapply(seq_along(fdr_cntrl), function(idx) {
    # combine full-sample and CV results
    table_fdr_c <- rbind(
      table_summary_fdr$table[[idx]],
      table_summary_fdr_cv$table[[idx]]
    ) %>%
    filter(method != "limma_cv") %>%
    as.data.table() %>%
    melt(id.vars = 1:3, measure.vars = c("fdp", "npv", "tpr", "tnr")) %>%
    mutate(
      n_samp = n_samp,
      method_label = case_when(
        method == "biotmle" ~ "EIF Moderated Variance",
        method == "classic_tmle" ~ "EIF Empirical Variance",
        method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
        method == "classic_tmle_cv" ~ "CV-EIF Empirical Variance",
        method == "limma" ~ "LIMMA Moderated Variance"
      ),
      metric_label = case_when(
        variable == "fdp" ~ "False Discovery Proportions",
        variable == "npv" ~ "Negative Predictive Values",
        variable == "tpr" ~ "True Positive Rates",
        variable == "tnr" ~ "True Negative Rates"
      ),
      est_type = ifelse(method == "limma", "Parametric", "Semiparametric"),
      fdr_label = paste("FDR =", fdr_cntrl[idx]),
      fdr_line = as.numeric(str_remove(fdr_label, "FDR = "))
    )
  }) %>% rbindlist()

  # combine summaries from non-sample-split and CV results
  summary_fdr_all <- lapply(seq_along(fdr_cntrl), function(idx) {
    # combine full-sample and CV results
    summary_fdr_c <- rbind(
      table_summary_fdr$summary[[idx]],
      table_summary_fdr_cv$summary[[idx]]
    ) %>%
    filter(method != "limma_cv") %>%
    as.data.table() %>%
    melt(id.vars = 1:2, measure.vars = "fdr") %>%
    mutate(
      n_samp = n_samp,
      method_label = case_when(
        method == "biotmle" ~ "EIF Moderated Variance",
        method == "classic_tmle" ~ "EIF Empirical Variance",
        method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
        method == "classic_tmle_cv" ~ "CV-EIF Empirical Variance",
        method == "limma" ~ "LIMMA Moderated Variance"
      ),
      est_type = ifelse(method == "limma", "Parametric", "Semiparametric"),
      fdr_label = paste("FDR =", fdr_cntrl[idx]),
      fdr_line = as.numeric(str_remove(fdr_label, "FDR = "))
    )
  }) %>% rbindlist()

  ## make FDR visualization
  p_fdr_main <- summary_fdr_all %>%
    filter(fdr_line == fdr_cutoff) %>%
    ggplot(aes(x = as.factor(n_samp), y = value, group = method_label,
               shape = method_label, fill = method_label)) +
    geom_line(linetype = "dotted", color = "black", size = 1,
              position = pd_error) +
    geom_point(size = 6, alpha = 0.8, position = pd_error) +
    geom_hline(aes(yintercept = fdr_line), colour = "red",
               linetype = "dashed") +
    labs(
      x = "",
      y = "False discovery rate",
      title = "FDR control of all candidate estimators",
      subtitle = paste0("(Nominal FDR = ", fdr_cutoff, ")")
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
    scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
    theme_bw() +
    theme(
      legend.background = element_rect(
        fill = "gray90", size = 0.25, linetype = "dotted"
      ),
      legend.position = "none",
      legend.title = element_blank(),
      text = element_text(size = 25),
      axis.text.x = element_text(angle = 25, hjust = 1)
    ) +
    guides(fill = guide_legend(nrow = 2))
  p_fdr_main_zoom <- summary_fdr_all %>%
    filter(fdr_line == fdr_cutoff, est_type == "Semiparametric") %>%
    ggplot(aes(x = as.factor(n_samp), y = value, group = method_label,
               shape = method_label, fill = method_label)) +
    geom_line(linetype = "dotted", color = "black", size = 1,
              position = pd_error) +
    geom_point(size = 6, alpha = 0.8, position = pd_error) +
    geom_hline(aes(yintercept = fdr_line), colour = "red",
               linetype = "dashed") +
    labs(
      x = "",
      y = "",
      title = "FDR control of efficient estimators",
      subtitle = "(zoomed in from left panel)"
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
    scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
    scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
    theme_bw() +
    theme(
      legend.background = element_rect(fill = "gray90", size = 0.25,
                                       linetype = "dotted"),
      legend.position = "none",
      legend.title = element_blank(),
      text = element_text(size = 25),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )

  # plots of empirical distributions of FDP, NPV, TNR, TPR
  p_emp_main <- table_fdr_all %>%
    filter(fdr_line == fdr_cutoff) %>%
    ggplot(aes(x = as.factor(n_samp), y = value,
               shape = method_label, fill = method_label)) +
    geom_point(size = 4, alpha = 0.01, position = pd_power) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, position = pd_power,
                 show.legend = FALSE) +
    facet_wrap(~ metric_label, nrow = 2, ncol = 2, scales = "fixed") +
    labs(
      x = "Sample size",
      y = "Empirical discovery rates",
      title = ""
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
    scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
    theme_bw() +
    theme(
      legend.background = element_rect(
        fill = "gray90", size = 0.25, linetype = "dotted"
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      text = element_text(size = 25),
      axis.text.x = element_text(angle = 25, hjust = 1)
    ) +
    guides(fill = guide_legend(
      override.aes = list(alpha = 1, size = 7))
    )

  # combine and save plots
  p_main <- (p_fdr_main | p_fdr_main_zoom) / p_emp_main +
    plot_annotation(
      title = paste("Variance moderation of efficient estimators",
                    "enhances control of FDR"),
      theme = theme(plot.title = element_text(size = 36))
    ) +
    plot_layout(heights = c(1, 2))
  ggsave(filename = here("graphs", "manuscript",
                         paste0("fdr_sl_null_logistic_", sim_dir, ".png")),
         plot = p_main, width = 22, height = 16, dpi = 400)
}


##############################################################################
# SIMULATION #2 - true findings, strong effect
##############################################################################

# results files and settings
for (sim_dir in c("g0prob_no", "g0prob_yes")) {
  if (sim_dir == "g0prob_no") {
    sim_files <- list(
      min = c(paste0("0.1-pval_logistic-std_eif-",
                     c("nsim_300-cvfolds_1-2020-04-15_09:03:36.rds",
                       "nsim_300-cvfolds_2-2020-04-16_21:12:54.rds"))),
      mod = c(paste0("0.3-pval_logistic-std_eif-",
                     c("nsim_300-cvfolds_1-2020-04-15_08:57:27.rds",
                       "nsim_300-cvfolds_2-2020-04-16_20:46:49.rds"))),
      max = c(paste0("0.5-pval_logistic-std_eif-",
                     c("nsim_300-cvfolds_1-2020-04-15_08:48:51.rds",
                       "nsim_300-cvfolds_2-2020-04-17_03:49:59.rds")))
    )
  } else if (sim_dir == "g0prob_yes") {
    sim_files <- list(
      min = c(paste0("0.1-pval_logistic-std_eif-",
                     c("nsim_300-cvfolds_1-2020-04-15_10:17:46.rds",
                       "nsim_300-cvfolds_2-2020-04-15_23:17:30.rds"))),
      mod = c(paste0("0.3-pval_logistic-std_eif-",
                     c("nsim_300-cvfolds_1-2020-04-17_23:10:16.rds",
                       "nsim_300-cvfolds_2-2020-04-18_23:54:38.rds"))),
      max = c(paste0("0.5-pval_logistic-std_eif-",
                     c("nsim_300-cvfolds_1-2020-04-15_09:55:22.rds",
                       "nsim_300-cvfolds_2-2020-04-16_09:06:37.rds")))
    )
  }

  # concatenate file prefixes and names
  file_prefix <- "sl-true_findings-effects_large-num_probes_"
  sim_files <- lapply(sim_files, function(file_name) {
    paste0(file_prefix, file_name)
  })

  # iteratively process simulation settings
  for (iter in seq_along(sim_files)) {
    # load data
    sim_file <- sim_files[[iter]][1]
    sim_file_cv <- sim_files[[iter]][2]
    no_findings <- nonnull_genes[iter]
    sim_results <- readRDS(here("data", paste0("main_", sim_dir), sim_file))
    sim_results_cv <- readRDS(here("data", paste0("main_", sim_dir),
                                   sim_file_cv))
    file_label <- results_label[iter]

    # properly compute adjusted p-values across methods + sample sizes
    sim_data <- sim_results %>%
      bind_rows(.id = "n_samp") %>%
      mutate(
        n_samp = as.numeric(str_remove(n_samp, "n_")),
      ) %>%
      group_by(method, n_samp, sim_id) %>%
      mutate(
        p_adj_fwer = p.adjust(p_val, method = "bonferroni"),
        p_adj_fdr = p.adjust(p_val, method = "BH")
      )
    no_sim <- length(unique(sim_data$sim_id))

    sim_data_cv <- sim_results_cv %>%
      bind_rows(.id = "n_samp") %>%
      mutate(
        n_samp = as.numeric(str_remove(n_samp, "n_")),
      ) %>%
      group_by(method, n_samp, sim_id) %>%
      mutate(
        p_adj_fwer = p.adjust(p_val, method = "bonferroni"),
        p_adj_fdr = p.adjust(p_val, method = "BH")
      ) %>%
      ungroup() %>%
      mutate(
        method = paste(method, "cv", sep = "_")
      ) %>%
      group_by(method, n_samp, sim_id)

    # summarize using utility functions
    table_summary_fdr <- make_summary_table(sim_data, fdr_cntrl)
    table_summary_fdr_cv <- make_summary_table(sim_data_cv, fdr_cntrl)

    # combine tables from non-sample-split and CV results
    table_fdr_all <- lapply(seq_along(fdr_cntrl), function(idx) {
      # combine full-sample and CV results
      table_fdr_c <- rbind(
        table_summary_fdr$table[[idx]],
        table_summary_fdr_cv$table[[idx]]
      ) %>%
      filter(method != "limma_cv") %>%
      as.data.table() %>%
      melt(id.vars = 1:3, measure.vars = c("fdp", "npv", "tpr", "tnr")) %>%
      mutate(
        n_samp = n_samp,
        method_label = case_when(
          method == "biotmle" ~ "EIF Moderated Variance",
          method == "classic_tmle" ~ "EIF Empirical Variance",
          method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
          method == "classic_tmle_cv" ~ "CV-EIF Empirical Variance",
          method == "limma" ~ "LIMMA Moderated Variance"
        ),
        metric_label = case_when(
          variable == "fdp" ~ "False Discovery Proportions",
          variable == "npv" ~ "Negative Predictive Values",
          variable == "tpr" ~ "True Positive Rates",
          variable == "tnr" ~ "True Negative Rates"
        ),
        est_type = ifelse(method == "limma", "Parametric", "Semiparametric"),
        fdr_label = paste("FDR =", fdr_cntrl[idx]),
        fdr_line = as.numeric(str_remove(fdr_label, "FDR = "))
      )
    }) %>% rbindlist()

    # combine summaries from non-sample-split and CV results
    summary_fdr_all <- lapply(seq_along(fdr_cntrl), function(idx) {
      # combine full-sample and CV results
      summary_fdr_c <- rbind(
        table_summary_fdr$summary[[idx]],
        table_summary_fdr_cv$summary[[idx]]
      ) %>%
      filter(method != "limma_cv") %>%
      as.data.table() %>%
      melt(id.vars = 1:2, measure.vars = "fdr") %>%
      mutate(
        n_samp = n_samp,
        method_label = case_when(
          method == "biotmle" ~ "EIF Moderated Variance",
          method == "classic_tmle" ~ "EIF Empirical Variance",
          method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
          method == "classic_tmle_cv" ~ "CV-EIF Empirical Variance",
          method == "limma" ~ "LIMMA Moderated Variance"
        ),
        est_type = ifelse(method == "limma", "Parametric", "Semiparametric"),
        fdr_label = paste("FDR =", fdr_cntrl[idx]),
        fdr_line = as.numeric(str_remove(fdr_label, "FDR = "))
      )
    }) %>% rbindlist()

    ## make FDR visualization
    p_fdr_main <- summary_fdr_all %>%
      filter(fdr_line == fdr_cutoff) %>%
      ggplot(aes(x = as.factor(n_samp), y = value, group = method_label,
                 shape = method_label, fill = method_label)) +
      geom_line(linetype = "dotted", color = "black", size = 1,
                position = pd_error) +
      geom_point(size = 6, alpha = 0.8, position = pd_error) +
      geom_hline(aes(yintercept = fdr_line), colour = "red",
                 linetype = "dashed") +
      labs(x = "",
           y = "False discovery rate",
           title = "FDR control of all candidate estimators",
           subtitle = paste0("(Nominal FDR = ", fdr_cutoff, ")")
          ) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
      scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
      scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
      theme_bw() +
      theme(legend.background = element_rect(
              fill = "gray90", size = 0.25, linetype = "dotted"
            ),
            legend.position = "none",
            legend.title = element_blank(),
            text = element_text(size = 25),
            axis.text.x = element_text(angle = 25, hjust = 1))
    p_fdr_main_zoom <- summary_fdr_all %>%
      filter(fdr_line == fdr_cutoff, est_type == "Semiparametric") %>%
      ggplot(aes(x = as.factor(n_samp), y = value, group = method_label,
                 shape = method_label, fill = method_label)) +
      geom_line(linetype = "dotted", color = "black", size = 1,
                position = pd_error) +
      geom_point(size = 6, alpha = 0.8, position = pd_error) +
      geom_hline(aes(yintercept = fdr_line), colour = "red",
                 linetype = "dashed") +
      labs(x = "",
           y = "",
           title = "FDR control of efficient estimators",
           subtitle = "(zoomed in from left panel)"
          ) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
      scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
      scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
      theme_bw() +
      theme(legend.background = element_rect(
              fill = "gray90", size = 0.25, linetype = "dotted"
            ),
            legend.position = "none",
            legend.title = element_blank(),
            text = element_text(size = 25),
            axis.text.x = element_text(angle = 25, hjust = 1))

    # plots of empirical distributions of FDP, NPV, TNR, TPR
    p_emp_main <- table_fdr_all %>%
      filter(fdr_line == fdr_cutoff) %>%
      ggplot(aes(x = as.factor(n_samp), y = value,
                 shape = method_label, fill = method_label)) +
      geom_point(size = 4, alpha = 0.01, position = pd_power) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, position = pd_power,
                   show.legend = FALSE) +
      facet_wrap(~ metric_label, nrow = 2, ncol = 2, scales = "fixed") +
      labs(
        x = "Sample size",
        y = "Empirical discovery rates",
        title = ""
      ) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
      scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
      scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
      theme_bw() +
      theme(
        legend.background = element_rect(fill = "gray90", size = 0.25,
                                         linetype = "dotted"),
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(size = 25),
        axis.text.x = element_text(angle = 25, hjust = 1)
      ) +
      guides(fill = guide_legend(
        override.aes = list(alpha = 1, size = 7))
      )

    # combine Type-I error and power plots
    p_main <- (p_fdr_main | p_fdr_main_zoom) / p_emp_main +
      plot_annotation(
        title = paste("Variance moderation of efficient estimators",
                      "enhances control of FDR"),
        theme = theme(plot.title = element_text(size = 36))
      ) +
      plot_layout(heights = c(1, 2))
    ggsave(filename = here("graphs", "manuscript",
                           paste0("fdr_sl_findings_", file_label,
                                  "_logistic_", sim_dir, ".png")),
           plot = p_main, width = 22, height = 16, dpi = 400)
  }
}


##############################################################################
# SIMULATION #3 - true findings, weak effect
##############################################################################

# results files and settings
sim_dir <- "main_g0prob_no"
sim_files <- list(
  min = c(paste0("0.1-pval_logistic-std_eif-",
                 c("nsim_300-cvfolds_1-2020-04-15_09:03:08.rds",
                   "nsim_300-cvfolds_2-2020-04-16_20:45:14.rds"))),
  mod = c(paste0("0.3-pval_logistic-std_eif-",
                 c("nsim_300-cvfolds_1-2020-04-15_16:01:28.rds",
                   "nsim_300-cvfolds_2-2020-04-17_07:37:59.rds"))),
  max = c(paste0("0.5-pval_logistic-std_eif-",
                 c("nsim_300-cvfolds_1-2020-04-15_18:23:47.rds",
                   "nsim_300-cvfolds_2-2020-04-16_20:33:16.rds")))
)
file_prefix <- "sl-true_findings-effects_small-num_probes_"
sim_files <- lapply(sim_files, function(file_name) {
  paste0(file_prefix, file_name)
})

# iteratively process simulation settings
for (iter in seq_along(sim_files)) {
  # load data
  sim_file <- sim_files[[iter]][1]
  sim_file_cv <- sim_files[[iter]][2]
  no_findings <- nonnull_genes[iter]
  sim_results <- readRDS(here("data", sim_dir, sim_file))
  sim_results_cv <- readRDS(here("data", sim_dir, sim_file_cv))
  file_label <- results_label[iter]

  # properly compute adjusted p-values across methods + sample sizes
  sim_data <- sim_results %>%
    bind_rows(.id = "n_samp") %>%
    mutate(
      n_samp = as.numeric(str_remove(n_samp, "n_")),
    ) %>%
    group_by(method, n_samp, sim_id) %>%
    mutate(
      p_adj_fwer = p.adjust(p_val, method = "bonferroni"),
      p_adj_fdr = p.adjust(p_val, method = "BH")
    )
  no_sim <- length(unique(sim_data$sim_id))

  sim_data_cv <- sim_results_cv %>%
    bind_rows(.id = "n_samp") %>%
    mutate(
      n_samp = as.numeric(str_remove(n_samp, "n_")),
    ) %>%
    group_by(method, n_samp, sim_id) %>%
    mutate(
      p_adj_fwer = p.adjust(p_val, method = "bonferroni"),
      p_adj_fdr = p.adjust(p_val, method = "BH")
    ) %>%
    ungroup() %>%
    mutate(
      method = paste(method, "cv", sep = "_")
    ) %>%
    group_by(method, n_samp, sim_id)

  # summarize using utility functions
  table_summary_fdr <- make_summary_table(sim_data, fdr_cntrl)
  table_summary_fdr_cv <- make_summary_table(sim_data_cv, fdr_cntrl)

  # combine tables from non-sample-split and CV results
  table_fdr_all <- lapply(seq_along(fdr_cntrl), function(idx) {
    # combine full-sample and CV results
    table_fdr_c <- rbind(
      table_summary_fdr$table[[idx]],
      table_summary_fdr_cv$table[[idx]]
    ) %>%
    filter(method != "limma_cv") %>%
    as.data.table() %>%
    melt(id.vars = 1:3, measure.vars = c("fdp", "npv", "tpr", "tnr")) %>%
    mutate(
      n_samp = n_samp,
      method_label = case_when(
        method == "biotmle" ~ "EIF Moderated Variance",
        method == "classic_tmle" ~ "EIF Empirical Variance",
        method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
        method == "classic_tmle_cv" ~ "CV-EIF Empirical Variance",
        method == "limma" ~ "LIMMA Moderated Variance"
      ),
      metric_label = case_when(
        variable == "fdp" ~ "False Discovery Proportions",
        variable == "npv" ~ "Negative Predictive Values",
        variable == "tpr" ~ "True Positive Rates",
        variable == "tnr" ~ "True Negative Rates"
      ),
      est_type = ifelse(method == "limma", "Parametric", "Semiparametric"),
      fdr_label = paste("FDR =", fdr_cntrl[idx]),
      fdr_line = as.numeric(str_remove(fdr_label, "FDR = "))
    )
  }) %>% rbindlist()

  # combine summaries from non-sample-split and CV results
  summary_fdr_all <- lapply(seq_along(fdr_cntrl), function(idx) {
    # combine full-sample and CV results
    summary_fdr_c <- rbind(
      table_summary_fdr$summary[[idx]],
      table_summary_fdr_cv$summary[[idx]]
    ) %>%
    filter(method != "limma_cv") %>%
    as.data.table() %>%
    melt(id.vars = 1:2, measure.vars = "fdr") %>%
    mutate(
      n_samp = n_samp,
      method_label = case_when(
        method == "biotmle" ~ "EIF Moderated Variance",
        method == "classic_tmle" ~ "EIF Empirical Variance",
        method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
        method == "classic_tmle_cv" ~ "CV-EIF Empirical Variance",
        method == "limma" ~ "LIMMA Moderated Variance"
      ),
      est_type = ifelse(method == "limma", "Parametric", "Semiparametric"),
      fdr_label = paste("FDR =", fdr_cntrl[idx]),
      fdr_line = as.numeric(str_remove(fdr_label, "FDR = "))
    )
  }) %>% rbindlist()

  ## make FDR visualization
  p_fdr_main <- summary_fdr_all %>%
    filter(fdr_line == fdr_cutoff) %>%
    ggplot(aes(x = as.factor(n_samp), y = value, group = method_label,
               shape = method_label, fill = method_label)) +
    geom_line(linetype = "dotted", color = "black", size = 1,
              position = pd_error) +
    geom_point(size = 6, alpha = 0.8, position = pd_error) +
    geom_hline(aes(yintercept = fdr_line), colour = "red",
               linetype = "dashed") +
    labs(x = "",
         y = "False discovery rate",
         title = "FDR control of all candidate estimators",
         subtitle = paste0("(Nominal FDR = ", fdr_cutoff, ")")
        ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
    scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
    theme_bw() +
    theme(legend.background = element_rect(
            fill = "gray90", size = 0.25, linetype = "dotted"
          ),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 25),
          axis.text.x = element_text(angle = 25, hjust = 1))
  p_fdr_main_zoom <- summary_fdr_all %>%
    filter(fdr_line == fdr_cutoff, est_type == "Semiparametric") %>%
    ggplot(aes(x = as.factor(n_samp), y = value, group = method_label,
               shape = method_label, fill = method_label)) +
    geom_line(linetype = "dotted", color = "black", size = 1,
              position = pd_error) +
    geom_point(size = 6, alpha = 0.8, position = pd_error) +
    geom_hline(aes(yintercept = fdr_line), colour = "red",
               linetype = "dashed") +
    labs(x = "",
         y = "",
         title = "FDR control of efficient estimators",
         subtitle = "(zoomed in from left panel)"
        ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
    scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
    scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
    theme_bw() +
    theme(legend.background = element_rect(
            fill = "gray90", size = 0.25, linetype = "dotted"
          ),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 25),
          axis.text.x = element_text(angle = 25, hjust = 1))

  # plots of empirical distributions of FDP, NPV, TNR, TPR
  p_emp_main <- table_fdr_all %>%
    filter(fdr_line == fdr_cutoff) %>%
    ggplot(aes(x = as.factor(n_samp), y = value,
               shape = method_label, fill = method_label)) +
    geom_point(size = 4, alpha = 0.01, position = pd_power) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, position = pd_power,
                 show.legend = FALSE) +
    facet_wrap(~ metric_label, nrow = 2, ncol = 2, scales = "fixed") +
    labs(
      x = "Sample size",
      y = "Empirical discovery rates",
      title = ""
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_shape_manual(values = c(24, 22, 24, 22, 21)) +
    scale_fill_manual(values = c(pal[1], pal[1], pal[3], pal[3], pal[5])) +
    theme_bw() +
    theme(
      legend.background = element_rect(
        fill = "gray90", size = 0.25, linetype = "dotted"
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      text = element_text(size = 25),
      axis.text.x = element_text(angle = 25, hjust = 1)
    ) +
    guides(fill = guide_legend(
      override.aes = list(alpha = 1, size = 7))
    )

  # combine Type-I error and power plots
  p_main <- (p_fdr_main | p_fdr_main_zoom) / p_emp_main +
    plot_annotation(
      title = paste("Variance moderation of efficient estimators",
                    "enhances control of FDR"),
      theme = theme(plot.title = element_text(size = 36))
    ) +
    plot_layout(heights = c(1, 2))
  ggsave(filename = here("graphs", "manuscript",
                         paste0("fdr_sl_findings_", file_label,
                                "_logistic_small_effect.png")),
         plot = p_main, width = 22, height = 16, dpi = 400)
}
