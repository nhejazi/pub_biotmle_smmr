###############################################################################
# Summarize Results
################################################################################

# load packages from renv library
renv::activate(project = here::here())

# load required packages and scripts
library(here)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)
library(GGally)
library(ggpubr)
library(patchwork)
library(MetBrewer)
source(here("R", "parallel_coord_plotting.R"))

# load the results
full_data <- read_rds(here("results", "full-dataset-biotmle.Rds"))
p25_data_ls <- read_rds(here("results", "p25-subsample-ls.Rds"))
p50_data_ls <- read_rds(here("results", "p50-subsample-ls.Rds"))
p75_data_ls <- read_rds(here("results", "p75-subsample-ls.Rds"))

# count the number of CpGs that are significant at the 5% FWER in each subsample
# application
p25_counts <- p25_data_ls %>% map_dbl(function(tbl) sum(tbl$adj_pval <= 0.05))
p50_counts <- p50_data_ls %>% map_dbl(function(tbl) sum(tbl$adj_pval <= 0.05))
p75_counts <- p75_data_ls %>% map_dbl(function(tbl) sum(tbl$adj_pval <= 0.05))

# create plot of distributions of CpG hits across subsamples
plot_tbl <- tibble(
    idx = seq_len(30),
    percentage = c(rep("25%", 10), rep("50%", 10), rep("75%", 10)),
    num_cpgs = c(p25_counts, p50_counts, p75_counts)
  ) %>%
  mutate(percentage = factor(percentage, levels = c("25%", "50%", "75%")))

total_hits_p <- plot_tbl %>%
  ggplot(aes(x = percentage, y = num_cpgs)) +
  geom_violin() +
  geom_point(size = 5, alpha = 0.5) +
  ylab("Number of Significant CpGs") +
  xlab("Subsample Size") +
  ggtitle("Comparison of Significant CpGs") +
  scale_y_log10() +
  theme_classic()

ggsave(
  filename = "subsample_results.jpg",
  path = here("graphs"),
  width = 4,
  height = 4,
  scale = 1.2
)

# create plot of min/median/max p-values of top CpGs across subsamples
top_ten_p <- plot_parallel_coord(
  full_data, p25_data_ls, p50_data_ls, p75_data_ls, rank_range = 1:10,
  title = "Ten Most Significant CpGs in Complete Analysis"
)
top_twenty_p <- plot_parallel_coord(
  full_data, p25_data_ls, p50_data_ls, p75_data_ls, rank_range = 11:20,
  title = "11-20 Most Significant CpGs in Complete Analysis"
)
top_thirty_p <- plot_parallel_coord(
  full_data, p25_data_ls, p50_data_ls, p75_data_ls, rank_range = 21:30,
  title = "21-30 Most Significant CpGs in Complete Analysis"
)
p_coord <- top_ten_p / top_twenty_p / top_thirty_p
ggsave(
  filename = "aggregated_par_coord.jpg",
  plot = p_coord,
  path = here("graphs"),
  width = 16,
  height = 12,
  scale = 1.2
)

# create plot of CpG rankings v median p-values across subsamples
p_cpg_ranks <- plot_cpg_ranks(
  full_data, p25_data_ls, p50_data_ls, p75_data_ls, rank_range = 1:30,
  title = "Top 30 most significant CpGs are stable to sampling perturbations"
)
ggsave(
  filename = "cpg_ranks.jpg",
  plot = p_cpg_ranks,
  path = here("graphs"),
  width = 16,
  height = 12,
  scale = 1.2
)
