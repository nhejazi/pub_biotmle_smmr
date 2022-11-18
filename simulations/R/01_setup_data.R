# simulate genomic data and return summarized experiment object
sim_genomic_data <- function(n_obs = 1000,  # sample size
                             p = 1000,      # number of probes
                             findings_type = c("global_null", "true_findings"),
                             positivity_issues = c("no", "yes"),
                             effect_size = c("large", "small"),
                             no_findings = NULL
                            ) {
  # argument munging
  findings_type <- match.arg(findings_type)
  positivity_issues <- match.arg(positivity_issues)
  effect_size <- match.arg(effect_size)

  # set number of truly differentially expressed probes
  if (findings_type != "global_null" & !is.null(no_findings)) {
    n_true_probes <- no_findings
  }

  # set effect size of exposure on differential expression
  if (effect_size == "large") {
    tx_coef <- 5
  } else if (effect_size == "small") {
    tx_coef <- 0.15
  }

  # baseline covariates
  W1 <- runif(n_obs)
  W2 <- runif(n_obs)

  # treatment model with built-in truncation
  if (positivity_issues == "yes") {
    g0 <- plogis(0.5 + 2.5 * W1 - 3 * W2 - 2)
  } else {
    g0 <- plogis(0.5 + 2.5 * W1 - 3 * W2)
  }
  A <- rbinom(n_obs, 1, prob = g0)

  if (findings_type == "global_null") {
    # outcome model for probes -- NO EFFECT OF A ON Y
    y_mean <- 2 + W1 + 0.5 * W2 + W1 * W2

    # generate probe data with specified covariance matrix
    #y_cov <- replicate(n_obs, runif(n = n_obs, min = 0, max = 1))
    #y_cov <- crossprod(y_cov)
    #y_mat <- rmvnorm(n = p, mean = y_mean, sigma = y_cov, method = "chol")

    # generate probe data as independent for each subject
    y_list <- lapply(seq_len(p), function(j) {
      y_mean +  rnorm(n = n_obs, mean = 0, sd = 1)
    })
    y_mat <- t(do.call(cbind, y_list))

  } else if (findings_type == "true_findings") {
    # outcome model for probes -- NO EFFECT OF A ON Y
    y_mean_null <- 2 + W1 + 0.5 * W2 + W1 * W2

    # generate probe data as independent for each subject
    y_list_null <- lapply(seq_len(p - n_true_probes), function(j) {
      y_mean_null + rnorm(n = n_obs, mean = 0, sd = 1)
    })
    y_mat_null <- t(do.call(cbind, y_list_null))

    # outcome model for probes -- TRUE EFFECT OF A ON Y
    y_mean_nonnull <- 2 + W1 + 0.5 * W2 + W1 * W2 + tx_coef * A

    # generate probe data as independent for each subject
    y_list_nonnull <- lapply(seq_len(n_true_probes), function(j) {
      y_mean_nonnull + rnorm(n = n_obs, mean = 0, sd = 0.2)
    })
    y_mat_nonnull <- t(do.call(cbind, y_list_nonnull))

    # combine probes into a single matrix
    y_mat <- rbind(y_mat_null, y_mat_nonnull)
  }


  # full set of probes
  probe_data <- y_mat %>%
    as.data.frame() %>%
    DataFrame()

  # phenotype and intervention variables
  pdata <- cbind(A, W1, W2) %>%
    as.data.frame() %>%
    mutate(
      A = I(A),
      W1 = I(W1),
      W2 = I(W2),
    ) %>%
    DataFrame()

  # setup as summarized experiment for compatibility
  data_se <- SummarizedExperiment(assays = list(probes = probe_data),
                                  colData = pdata)
  return(data_se)
}
