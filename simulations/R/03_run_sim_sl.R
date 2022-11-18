# use custom package library on Savio cluster
if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/nhejazi/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# read in command line arguments and set defaults
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                             defaults = list(# effect type, eg, "true_findings"
                                             sim_type = "global_null",
                                             # should also try out V-fold CV
                                             cv_folds = 1,
                                             # different SL libraries
                                             sl_spec = "main",  # "supp"
                                             # reference dist. for test stat.
                                             stat_ref = "logistic",  # "normal"
                                             # standard v corrupted EIF est.
                                             eif_type = "std_eif",
                                             # number of simulation
                                             n_sim = 300,
                                             # number of probes
                                             p_genes = 150,
                                             # proportion of true findings
                                             prop_genes = 0,
                                             # positivity violations in g0
                                             positivity_issues = "no",
                                             # size of DE effect
                                             effect_size = "large"))

# reference for logging
print(args)

# packages
library(here)
library(tidyverse)
library(mvtnorm)
library(foreach)
library(future)
library(doFuture)
library(doRNG)
library(hal9001)
library(SuperLearner)
library(drtmle)
library(limma)
library(SummarizedExperiment)
library(biotmle)

# parallelization
options(future.globals.maxSize = 1e32)
registerDoFuture()
plan(multiprocess, workers = 24L)
set.seed(7491)

# load scripts, parallelization, PRNG
source(here("R", "01_setup_data.R"))
source(here("R", "02_compare_methods.R"))
seed_int_tmle <- 94153
n_obs <- c(50, 100, 200, 400)    # sample sizes doubling each time

# number of findings when effect type is not the global null
if (args$sim_type == "global_null") {
  findings_total <- NULL
} else {
  findings_total <- round(args$p_genes * args$prop_genes, 1)
}

# use different Super Learner specifications by sample size
sl_lib_main <- list(
  g_lib_small = c("SL.mean", "SL.bayesglm", "SL.glmnet", "SL.speedglm"),
  g_lib = c("SL.mean", "SL.bayesglm", "SL.glmnet", "SL.speedglm", "SL.gam",
            "SL.randomForest"),
  Q_lib_small = c("SL.mean", "SL.bayesglm", "SL.glmnet", "SL.speedglm"),
  Q_lib = c("SL.mean", "SL.bayesglm", "SL.glmnet", "SL.speedglm", "SL.earth",
            "SL.randomForest", "SL.xgboost")
)

# NOTE: the following library works too well for this DGP --- making the effect
#       of variance moderation essentially negligible
# NOTE: As of 18 January, this gives perfect results -- do NOT change this
sl_lib_supp <- list(
  g_lib_small = c("SL.mean", "SL.speedglm", "SL.bayesglm"),
  g_lib = c("SL.mean", "SL.speedglm", "SL.bayesglm", "SL.gam"),
  Q_lib_small = c("SL.mean", "SL.bayesglm", "SL.glm.interaction"),
  Q_lib = c("SL.mean", "SL.bayesglm", "SL.glm.interaction", "SL.earth")
)


# well-specified or mainified simulation setting
if (args$sl_spec == "supp") {
  g_lib <- sl_lib_supp$g_lib
  Q_lib <- sl_lib_supp$Q_lib
  g_lib_small <- sl_lib_supp$g_lib_small
  Q_lib_small <- sl_lib_supp$Q_lib_small
} else if (args$sl_spec == "main") {
  g_lib <- sl_lib_main$g_lib
  Q_lib <- sl_lib_main$Q_lib
  g_lib_small <- sl_lib_main$g_lib_small
  Q_lib_small <- sl_lib_main$Q_lib_small
}

# perform simulation across sample sizes
sim_results <- lapply(n_obs, function(sample_size) {
  # get results in parallel
  results <- foreach(sim_iter = seq_len(args$n_sim),
                     .options.multicore = list(preschedule = FALSE),
                     .errorhandling = "remove") %dorng% {
    # simulate data
    data_sumexp <- sim_genomic_data(n_obs = sample_size,
                                    p = args$p_genes,
                                    findings_type = args$sim_type,
                                    positivity_issues = args$positivity_issues,
                                    effect_size = args$effect_size,
                                    no_findings = findings_total)

    # set library based on sample size
    if (sample_size < 150) {
      g_lib <- g_lib_small
      Q_lib <- Q_lib_small
    }

    # fit each estimator
    out <- fit_estimators(data_se = data_sumexp,
                          stat_ref = args$stat_ref,
                          eif_type = args$eif_type,
                          no_findings = findings_total,
                          seed_int = seed_int_tmle,
                          # followings dot args passed to drtmle and biotmle
                          g_lib = g_lib,
                          Q_lib = Q_lib,
                          cv_folds = args$cv_folds,
                          tolg = 0.00001)

    # output for logging
    print(paste("Finished iteration", sim_iter, "for n =", sample_size))

    # re-organize + add info for output
    out <- out %>%
      mutate(
        sim_id = rep(sim_iter, n())
      )
    return(out)
  }

  # concatenate list into tibble before moving to next sample size
  results <- bind_rows(results)
})

# save results to file
names(sim_results) <- paste("n", n_obs, sep = "_")
timestamp <- str_replace_all(Sys.time(), " ", "_")
if (args$sim_type == "global_null") {
  saveRDS(object = sim_results,
          file = here("data",
                      paste0(args$sl_spec, "_g0prob_", args$positivity_issues),
                      paste0("sl-", args$sim_type,
                             "-pval_", args$stat_ref,
                             "-", args$eif_type,
                             "-nsim_", args$n_sim,
                             "-cvfolds_", args$cv_folds,
                             "-", timestamp, ".rds")))
} else {
  frac_findings <- round(findings_total / args$p_genes, 1)
  saveRDS(object = sim_results,
          file = here("data",
                      paste0(args$sl_spec, "_g0prob_", args$positivity_issues),
                      paste0("sl-", args$sim_type, "-",
                             "effects_", args$effect_size,
                             "-num_probes_", frac_findings,
                             "-pval_", args$stat_ref,
                             "-", args$eif_type,
                             "-nsim_", args$n_sim,
                             "-cvfolds_", args$cv_folds,
                             "-", timestamp, ".rds")))
}
