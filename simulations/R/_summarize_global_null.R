# NOTE: on average, proportion of false positive should be specified rate
library(here)
library(biotmle)
library(data.table)
library(tidyverse)
library(latex2exp)
library(patchwork)
library(ggsci)
pd <- position_dodge(0.2)

###############################################################################
# SIMULATION DATA - global null
###############################################################################

# results files and settings
sim_dir <- "g0prob_no"
sim_files <- paste0("sl-global_null-pval_logistic-std_eif-",
                   c("nsim_300-cvfolds_1-2020-04-15_08:56:56.rds",
                     "nsim_300-cvfolds_2-2020-04-20_18:33:13.rds"))
sim_dir <- "g0prob_yes"
sim_files <- paste0("sl-global_null-pval_logistic-std_eif-",
                   c("nsim_300-cvfolds_1-2020-04-15_20:21:05.rds",
                     "nsim_300-cvfolds_2-2020-04-16_20:31:34.rds"))

# read in results
sim_results <- readRDS(here("data", paste0("main_", sim_dir), sim_files[1]))
sim_results_cv <- readRDS(here("data", paste0("main_", sim_dir), sim_files[2]))

# how many simulations per sample size per estimator?
fwer_cntrl <- c(0.01, 0.05, 0.1)
fdr_cntrl <- c(0.01, 0.05, 0.1)
true_psi <- 0   # A has no effect on Y => ATE = 0

# properly compute adjusted p-values across estimation methods / sample sizes
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

###############################################################################
# SUMMARIZE FWER - global null
###############################################################################

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
      fwer_low = mean(fwer_count_low > 0),
      fwer_med = mean(fwer_count_med > 0),
      fwer_high = mean(fwer_count_high > 0)
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
             position = pd) +
  geom_line(linetype = "dotted", size = 1.5) +
  geom_hline(aes(yintercept = fwer_line), colour = "red", linetype = "dashed") +
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
ggsave(filename = here("graphs", sim_dir, "vis_fwer_sl_null_logistic.pdf"),
       plot = p_fwer, width = 21, height = 14)

###############################################################################
# SUMMARIZE FDR - global null
###############################################################################

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
             position = pd) +
  geom_line(linetype = "dotted", size = 1.5) +
  geom_hline(aes(yintercept = fdr_line), colour = "red", linetype = "dashed") +
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
ggsave(filename = here("graphs", sim_dir, "vis_fdr_sl_null_logistic.pdf"),
       plot = p_fdr, width = 21, height = 14)
