# packages and programmatic housekeeping
library(here)
library(biotmle)
library(data.table)
library(tidyverse)
library(latex2exp)
library(patchwork)
library(ggsci)
pd_error <- position_dodge(0.2)
pd_power <- position_dodge(0.5)

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

# simulation parameters
fwer_cntrl <- c(0.01, 0.05, 0.1)
fdr_cntrl <- c(0.01, 0.05, 0.1)
nonnull_genes <- 150 * c(0.1, 0.3, 0.5)
results_label <- c("min", "mod", "max")

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

  #############################################################################
  # SUMMARIZE FWER - true findings included
  #############################################################################

  # make FWER table by properly computing rate across simulations
  make_fwer_table <- function(sim_data) {
    table_fwer <- sim_data %>%
      dplyr::filter(!is.na(psi)) %>%
      summarise(
        fwer_count_low = sum(p_adj_fwer < fwer_cntrl[1]),
        fwer_count_med = sum(p_adj_fwer < fwer_cntrl[2]),
        fwer_count_high = sum(p_adj_fwer < fwer_cntrl[3])
      ) %>%
      ungroup(sim_id) %>%
      group_by(method, n_samp) %>%
      summarise(
        fwer_low = mean(fwer_count_low > no_findings),
        fwer_med = mean(fwer_count_med > no_findings),
        fwer_high = mean(fwer_count_high > no_findings)
      )
    return(table_fwer)
  }
  table_fwer <- make_fwer_table(sim_data)
  table_fwer_cv <- make_fwer_table(sim_data_cv)
  table_fwer_all <- rbind(table_fwer, table_fwer_cv) %>%
    dplyr::filter(method != "limma_cv")

  # make FWER visualization
  p_fwer <- table_fwer_all %>%
    reshape2::melt(id.vars = 1:2) %>%
    transmute(
      n_samp = n_samp,
      method = case_when(
        method == "biotmle" ~ "EIF Moderated Variance",
        method == "classic_tmle" ~ "EIF Standard Variance",
        method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
        method == "classic_tmle_cv" ~ "CV-EIF Standard Variance",
        method == "limma" ~ "Moderated Variance w/ LIMMA"
      ),
      fwer_label = paste("FWER =",
                         case_when(variable == "fwer_high" ~ fwer_cntrl[3],
                                   variable == "fwer_med" ~ fwer_cntrl[2],
                                   variable == "fwer_low" ~ fwer_cntrl[1])
                        ),
      fwer_line = as.numeric(str_remove(fwer_label, "FWER = ")),
      fwer = value
    ) %>%
    ggplot(aes(x = as.factor(n_samp), y = fwer, group = method)) +
    geom_point(aes(shape = method, fill = method), size = 8, alpha = 0.6,
               position = pd_error) +
    geom_line(linetype = "dotted", size = 1.5, position = pd_error) +
    geom_hline(aes(yintercept = fwer_line), colour = "red",
               linetype = "dashed") +
    facet_grid(~fwer_label, scales = "fixed") +
    labs(x = "Sample size",
         y = paste("Family-wise error rate (across", no_sim, "simulations)"),
         title = paste("Variance moderation of efficient estimator enhances",
                       "control of family-wise error rate")
        ) +
    scale_shape_manual(values = c(24, 22, 25, 23, 21)) +
    scale_fill_nejm() +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = .25,
                                           linetype = "dotted"),
          legend.position = "bottom",
          legend.title = element_blank(),
          text = element_text(size = 25),
          axis.text.x = element_text(angle = 25, hjust = 1))

  # save plot
  ggsave(filename = here("graphs", "effect_small",
                         paste0("vis_fwer_sl_findings_", file_label,
                                "_logistic.pdf")),
         plot = p_fwer, width = 21, height = 14)

  #############################################################################
  # SUMMARIZE FDR - true findings included
  #############################################################################

  # make FDR table by properly computing rate across simulations
  make_fdr_table <- function(sim_data) {
    table_fdr <- sim_data %>%
      dplyr::filter(!is.na(psi)) %>%
      mutate(
        reject_low = (p_adj_fdr < fdr_cntrl[1]),
        reject_med = (p_adj_fdr < fdr_cntrl[2]),
        reject_high = (p_adj_fdr < fdr_cntrl[3])
      ) %>%
      mutate(
        fdp_low = ifelse(sum(reject_low) > 0,
                         sum(label == "null" & reject_low) /
                             max(1, sum(reject_low)), 0),
        fdp_med = ifelse(sum(reject_med) > 0,
                         sum(label == "null" & reject_med) /
                             max(1, sum(reject_med)), 0),
        fdp_high = ifelse(sum(reject_high) > 0,
                          sum(label == "null" & reject_high) /
                              max(1, sum(reject_high)), 0)
      ) %>%
      ungroup(sim_id) %>%
      group_by(method, n_samp) %>%
      summarise(
        fdr_low = mean(fdp_low),
        fdr_med = mean(fdp_med),
        fdr_high = mean(fdp_high)
      )
    return(table_fdr)
  }
  table_fdr <- make_fdr_table(sim_data)
  table_fdr_cv <- make_fdr_table(sim_data_cv)
  table_fdr_all <- rbind(table_fdr, table_fdr_cv) %>%
    dplyr::filter(method != "limma_cv")

  # make FDR visualization
  p_fdr <- table_fdr_all %>%
    reshape2::melt(id.vars = 1:2) %>%
    transmute(
      n_samp = n_samp,
      method = case_when(
        method == "biotmle" ~ "EIF Moderated Variance",
        method == "classic_tmle" ~ "EIF Standard Variance",
        method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
        method == "classic_tmle_cv" ~ "CV-EIF Standard Variance",
        method == "limma" ~ "Moderated Variance w/ LIMMA"
      ),
      fdr_label = paste("FDR =",
                        case_when(variable == "fdr_high" ~ fdr_cntrl[3],
                                  variable == "fdr_med" ~ fdr_cntrl[2],
                                  variable == "fdr_low" ~ fdr_cntrl[1])
                       ),
      fdr_line = as.numeric(str_remove(fdr_label, "FDR = ")),
      fdr = value
    ) %>%
    ggplot(aes(x = as.factor(n_samp), y = fdr, group = method)) +
    geom_point(aes(shape = method, fill = method), size = 8, alpha = 0.6,
               position = pd_error) +
    geom_line(linetype = "dotted", size = 1.5, position = pd_error) +
    geom_hline(aes(yintercept = fdr_line), colour = "red",
               linetype = "dashed") +
    facet_grid(~fdr_label, scales = "free_y") +
    labs(x = "Sample size",
         y = paste("False discovery rate (across", no_sim, "simulations)"),
         title = paste("Variance moderation of efficient estimator enhances",
                       "control of false discovery rate")
        ) +
    scale_shape_manual(values = c(24, 22, 25, 23, 21)) +
    scale_fill_nejm() +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = .25,
                                           linetype = "dotted"),
          legend.position = "bottom",
          legend.title = element_blank(),
          text = element_text(size = 25),
          axis.text.x = element_text(angle = 25, hjust = 1))

  # save plot
  ggsave(filename = here("graphs", "effect_small",
                         paste0("vis_fdr_sl_findings_", file_label,
                                "_logistic.pdf")),
         plot = p_fdr, width = 21, height = 14)

  #############################################################################
  # SUMMARIZE POWER - true findings only
  #############################################################################

  # make power table by looking at only true findings 
  make_power_table <- function(sim_data, err_type = c("fdr", "fwer"),
                               p_cntrl) {
    table_power <- sim_data %>%
      dplyr::filter(!is.na(psi)) %>%
      dplyr::filter(label == "hit") %>%
      group_by(n_samp, method, gene_id) %>%
      summarise(
        reject_err_low = mean(get(paste0("p_adj_", err_type)) < p_cntrl[1]),
        reject_err_med = mean(get(paste0("p_adj_", err_type)) < p_cntrl[2]),
        reject_err_high = mean(get(paste0("p_adj_", err_type)) < p_cntrl[3]),
      ) %>%
      summarise(
        power_err_low = median(reject_err_low),
        power_err_med = median(reject_err_low),
        power_err_high = median(reject_err_low)
      )
    return(table_power)
  }

  table_power <- make_power_table(sim_data, "fdr", fdr_cntrl)
  table_power_cv <- make_power_table(sim_data_cv, "fdr", fdr_cntrl)
  table_power_c <- rbind(table_power, table_power_cv) %>%
    dplyr::filter(method != "limma_cv")

  p_power <- table_power_c %>%
    reshape2::melt(id.vars = 1:2) %>%
    transmute(
      n_samp = n_samp,
      method = case_when(
        method == "biotmle" ~ "EIF Moderated Variance",
        method == "classic_tmle" ~ "EIF Standard Variance",
        method == "biotmle_cv" ~ "CV-EIF Moderated Variance",
        method == "classic_tmle_cv" ~ "CV-EIF Standard Variance",
        method == "limma" ~ "Moderated Variance w/ LIMMA"
      ),
      fdr_label = paste("FDR =",
                        case_when(variable == "power_err_high" ~ fdr_cntrl[3],
                                  variable == "power_err_med" ~ fdr_cntrl[2],
                                  variable == "power_err_low" ~ fdr_cntrl[1])
                       ),
      fdr_line = as.numeric(str_remove(fdr_label, "FDR = ")),
      fdr = value
    ) %>%
    ggplot(aes(x = as.factor(n_samp), y = fdr, group = method)) +
      geom_point(aes(shape = method, fill = method), size = 8, alpha = 0.6,
                 position = pd_power) +
      geom_line(linetype = "dotted", size = 1.5, position = pd_power) +
      facet_grid(~fdr_label, scales = "free_y") +
      labs(x = "Sample size",
           y = paste("Minimal recovery rate (across", no_sim, "simulations)"),
           title = paste("Variance moderation preserves worst-case power",
                         "while controlling false discovery rate")
          ) +
      scale_shape_manual(values = c(24, 22, 25, 23, 21)) +
      scale_fill_nejm() +
      theme_bw() +
      theme(legend.background = element_rect(fill = "gray90", size = .25,
                                             linetype = "dotted"),
            legend.position = "bottom",
            legend.title = element_blank(),
            text = element_text(size = 25),
            axis.text.x = element_text(angle = 25, hjust = 1))

  # save plot
  ggsave(filename = here("graphs", "effect_small",
                         paste0("vis_power_sl_findings_", file_label,
                                "_logistic.pdf")),
         plot = p_power, width = 21, height = 14)
}
