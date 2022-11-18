################################################################################
# Parallel Coordiante Plots
################################################################################

plot_parallel_coord <- function(
  full_data, p25_data_ls, p50_data_ls, p_75_data_ls, rank_range, title
) {

  # create the dataframe for plotting
  par_coord_df <- full_data %>%
    arrange(adj_pval) %>%
    mutate(full_adj_pval = -log10(adj_pval)) %>%
    select(ID, full_adj_pval) %>%
    slice(rank_range)
  colnames(par_coord_df) <- c("ID", "median_pval")
  cpg_order <- par_coord_df %>% pull(ID)

  p25_par_coord_df <- p25_data_ls %>%
    map2_dfc(seq_along(p25_data_ls),
             function(df, idx) {
               df <- df %>%
                 select(ID, adj_pval) %>%
                 right_join(par_coord_df, by = "ID") %>%
                 mutate(p25_adj_pval = -log10(
                     if_else(is.na(adj_pval), 1, adj_pval))
                   ) %>%
                 select(p25_adj_pval)
               colnames(df) <- paste("25% Sub. Data", idx)
               df
             }
    )
  p25_median_pval <- p25_par_coord_df %>% apply(1, median)
  p25_min_pval <- p25_par_coord_df %>% apply(1, min)
  p25_max_pval <- p25_par_coord_df %>% apply(1, max)
  p25_df <- data.frame(
    ID = par_coord_df$ID,
    subset = "25% Subsample",
    median_pval = p25_median_pval,
    min_pval = p25_min_pval,
    max_pval = p25_max_pval
  )

  p50_par_coord_df <- p50_data_ls %>%
    map2_dfc(seq_along(p50_data_ls),
             function(df, idx) {
               df <- df %>%
                 select(ID, adj_pval) %>%
                 right_join(par_coord_df, by = "ID") %>%
                 mutate(p50_adj_pval = -log10(
                     if_else(is.na(adj_pval), 1, adj_pval))
                   ) %>%
                 select(p50_adj_pval)
               colnames(df) <- paste("50% Sub. Data", idx)
               df
             }
    )
  p50_median_pval <- p50_par_coord_df %>% apply(1, median)
  p50_min_pval <- p50_par_coord_df %>% apply(1, min)
  p50_max_pval <- p50_par_coord_df %>% apply(1, max)
  p50_df <- data.frame(
    ID = par_coord_df$ID,
    subset = "50% Subsample",
    median_pval = p50_median_pval,
    min_pval = p50_min_pval,
    max_pval = p50_max_pval
  )

  p75_par_coord_df <- p75_data_ls %>%
    map2_dfc(seq_along(p75_data_ls),
             function(df, idx) {
               df <- df %>%
                 select(ID, adj_pval) %>%
                 right_join(par_coord_df, by = "ID") %>%
                 mutate(p75_adj_pval = -log10(
                     if_else(is.na(adj_pval), 1, adj_pval))
                   ) %>%
                 select(p75_adj_pval)
               colnames(df) <- paste("75% Sub. Data", idx)
               df
             }
    )
  p75_median_pval <- p75_par_coord_df %>% apply(1, median)
  p75_min_pval <- p75_par_coord_df %>% apply(1, min)
  p75_max_pval <- p75_par_coord_df %>% apply(1, max)
  p75_df <- data.frame(
    ID = par_coord_df$ID,
    subset = "75% Subsample",
    median_pval = p75_median_pval,
    min_pval = p75_min_pval,
    max_pval = p75_max_pval
  )

  par_coord_df <- par_coord_df %>%
    mutate(
      subset = "Full Data",
      min_pval = median_pval,
      max_pval = median_pval
    )

  par_coord_df <- par_coord_df %>%
    bind_rows(p75_df, p50_df, p25_df) %>%
    mutate(
      subset = factor(subset, levels = c("Full Data", "75% Subsample",
                                         "50% Subsample", "25% Subsample")),
      ID = factor(ID, levels = cpg_order)
    )

  # create the parallel coordinates plot
  par_min <- par_coord_df %>%
    ggplot(aes(x = subset, y = min_pval, group = ID, colour = ID)) +
    geom_line() +
    geom_hline(
      yintercept = -log10(0.05), colour = "black", linetype = 2, alpha = 0.3
    ) +
    ylab("Min -log(adj. p-value)") +
    xlab("") +
    ggtitle("") +
    scale_color_discrete(name = "CpG Site") +
    scale_y_log10(limits = c(0.1, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 16),
          axis.text.y = element_text(size = 16),
          text = element_text(size = 20))
  par_median <- par_coord_df %>%
    ggplot(aes(x = subset, y = median_pval, group = ID, colour = ID)) +
    geom_line() +
    geom_hline(
      yintercept = -log10(0.05), colour = "black", linetype = 2, alpha = 0.3
    ) +
    ylab("Median -log(adj. p-value)") +
    xlab("Datasets") +
    ggtitle(title) +
    scale_color_discrete(name = "CpG Site") +
    scale_y_log10(limits = c(0.1, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 16),
          axis.text.y = element_text(size = 16),
          text = element_text(size = 20))
  par_max <- par_coord_df %>%
    ggplot(aes(x = subset, y = max_pval, group = ID, colour = ID)) +
    geom_line() +
    geom_hline(
      yintercept = -log10(0.05), colour = "black", linetype = 2, alpha = 0.3
    ) +
    ylab("Max -log(adj. p-value)") +
    xlab("") +
    ggtitle("") +
    scale_color_discrete(name = "CpG Site") +
    scale_y_log10(limits = c(0.1, 100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 16),
          axis.text.y = element_text(size = 16),
          text = element_text(size = 20))

  # assemble the panels
  ggarrange(
    par_min, par_median, par_max, nrow = 1, common.legend = TRUE,
    legend = "right"
  )
}

plot_cpg_ranks <- function(
  full_data, p25_data_ls, p50_data_ls, p_75_data_ls, rank_range, title
) {

  # create the dataframe for plotting
  par_coord_df <- full_data %>%
    arrange(adj_pval) %>%
    mutate(full_adj_pval = -log10(adj_pval)) %>%
    select(ID, full_adj_pval) %>%
    slice(rank_range)
  colnames(par_coord_df) <- c("ID", "median_pval")
  cpg_order <- par_coord_df %>% pull(ID)
  
  p25_par_coord_df <- p25_data_ls %>%
    map2_dfc(seq_along(p25_data_ls),
             function(df, idx) {
               df <- df %>%
                 select(ID, adj_pval) %>%
                 right_join(par_coord_df, by = "ID") %>%
                 mutate(p25_adj_pval = -log10(
                   if_else(is.na(adj_pval), 1, adj_pval))
                 ) %>%
                 select(p25_adj_pval)
               colnames(df) <- paste("25% Sub. Data", idx)
               df
             }
    )
  p25_median_pval <- p25_par_coord_df %>% apply(1, median)
  p25_df <- data.frame(
    ID = par_coord_df$ID,
    subset = "25% Subsamples",
    median_pval = p25_median_pval
  )
  
  p50_par_coord_df <- p50_data_ls %>%
    map2_dfc(seq_along(p50_data_ls),
             function(df, idx) {
               df <- df %>%
                 select(ID, adj_pval) %>%
                 right_join(par_coord_df, by = "ID") %>%
                 mutate(p50_adj_pval = -log10(
                   if_else(is.na(adj_pval), 1, adj_pval))
                 ) %>%
                 select(p50_adj_pval)
               colnames(df) <- paste("50% Sub. Data", idx)
               df
             }
    )
  p50_median_pval <- p50_par_coord_df %>% apply(1, median)
  p50_df <- data.frame(
    ID = par_coord_df$ID,
    subset = "50% Subsamples",
    median_pval = p50_median_pval
  )
  
  p75_par_coord_df <- p75_data_ls %>%
    map2_dfc(seq_along(p75_data_ls),
             function(df, idx) {
               df <- df %>%
                 select(ID, adj_pval) %>%
                 right_join(par_coord_df, by = "ID") %>%
                 mutate(p75_adj_pval = -log10(
                   if_else(is.na(adj_pval), 1, adj_pval))
                 ) %>%
                 select(p75_adj_pval)
               colnames(df) <- paste("75% Sub. Data", idx)
               df
             }
    )
  p75_median_pval <- p75_par_coord_df %>% apply(1, median)
  p75_df <- data.frame(
    ID = par_coord_df$ID,
    subset = "75% Subsamples",
    median_pval = p75_median_pval
  )

  par_coord_df <- par_coord_df %>%
    mutate(subset = "Full Data") %>%
    bind_rows(p75_df, p50_df, p25_df) %>%
    mutate(
      subset = factor(subset, levels = c("Full Data", "75% Subsamples",
                                         "50% Subsamples", "25% Subsamples")),
      ID = factor(ID, levels = cpg_order)
    )

  # generate the plot
  par_coord_df %>%
    ggplot(aes(y = median_pval, x = ID, group = subset)) +
      geom_line(linetype = "dotted") +
      geom_point(size = 7, aes(color = subset)) +
      geom_hline(yintercept = -log10(0.05), colour = "black", linetype = 2) +
      labs(
        title = title,
        x = "CpG Site",
        y = "Median -log(Adj. P-Value)"
      ) +
      scale_color_manual(
        name = "Datasets",
        values = met.brewer("Cassatt2", type = "discrete")[c(1, 3, 4, 5)]
      ) +
      theme_bw() +
      theme(
        legend.background = element_rect(
          fill = "gray90", size = 0.25, linetype = "dotted"
        ),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        text = element_text(size = 20)
      )
}
