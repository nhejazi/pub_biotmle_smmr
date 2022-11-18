################################################################################
# function to fit various estimators
################################################################################
fit_estimators <- function(data_se,
                           adjust_method = c("BH", "bonferroni"),
                           g_lib,
                           Q_lib,
                           cv_folds,
                           stat_ref = c("logistic", "normal"),
                           eif_type = c("std_eif", "corrupt_eif"),
                           no_findings = NULL,
                           seed_int,
                           ...) {
  # get utility function for corrupted influence function method
  source(here::here("R", "corrupted_eif.R"))

  # set defaults for simulation
  adjust_method <- match.arg(adjust_method)
  stat_ref <- match.arg(stat_ref)
  eif_type <- match.arg(eif_type)

  # testing biotmle
  set.seed(seed_int)
  suppressWarnings(
    biotmle_fit <- biomarkertmle(se = data_se, varInt = 1, parallel = FALSE,
                                 bppar_debug = TRUE, cv_folds = cv_folds,
                                 g_lib = g_lib, Q_lib = Q_lib, ...)
  )

  # corrupt influence function if so set
  if (eif_type == "corrupt_eif") {
    eif_corrupt <- apply(biotmle_fit@tmleOut, 1, corrupt_eif_sim, delta = 0.1)
    eif_corrupt <- as_tibble(t(eif_corrupt) + rowMeans(biotmle_fit@tmleOut))
    biotmle_fit@tmleOut <- eif_corrupt
  }

  # moderate variance estimate
  biotmle_fit <- modtest_ic(biotmle = biotmle_fit, adjust = adjust_method,
                            confint = 0.95)

  # HOW DOES LIMMA COMPUTE P-VALUES?
  #p.value <- 2*pt(-abs(test$t),df=df.total)
  biotmle_out <- biotmle_fit@topTable %>%
    transmute(
      gene_id = seq_len(n()),
      ci_lwr = CI.L,
      psi = AveExpr,
      ci_upr = CI.R,
      se_psi = sqrt(var_bayes),
      p_val = case_when(stat_ref == "normal" ~ P.Value,
                        stat_ref == "logistic" ~
                          2 * plogis(-abs(t), location = 0,
                                     scale = sqrt(3) / pi)
                       ),
      method = rep("biotmle", n()),
      label =
        if (is.null(no_findings)) {
          rep("null", nrow(data_se))
        } else {
          c(rep("null", nrow(data_se) - no_findings), rep("hit", no_findings))
        }
    )

  # testing standard tmle
  A <- as.numeric(colData(data_se)$A)
  W <- as.data.frame(colData(data_se)[, -1])
  Y <- as.data.frame(t(as.matrix(assay(data_se))))
  a_0 <- sort(unique(A[!is.na(A)]))
  set.seed(seed_int)
  classic_tmle_fit <- lapply(X = seq_along(Y),
                             FUN = function(gene) {
                               suppressWarnings(
                                 tmle_fit <- drtmle(Y = Y[, gene],
                                                    A = A, W = W,
                                                    a_0 = a_0,
                                                    SL_g = g_lib,
                                                    SL_Q = Q_lib,
                                                    cvFolds = cv_folds,
                                                    stratify = TRUE,
                                                    guard = NULL,
                                                    parallel = FALSE,
                                                    use_future = FALSE,
                                                    ...)
                               )
                               psi <- tmle_fit$tmle$est[2] -
                                 tmle_fit$tmle$est[1]
                               eif_delta <- tmle_fit$ic$ic[, 2] -
                                 tmle_fit$ic$ic[, 1]
                               var_psi <- var(eif_delta) / length(A)
                               gene_id <- gene
                               return(cbind(gene_id, psi, var_psi))
                             })
  classic_tmle_fit <- do.call(rbind, classic_tmle_fit)

  # use variance from corrupted influence function
  if (eif_type == "corrupt_eif") {
    n_samp <- length(colData(data_se)$A)
    var_eif_corrupt <- apply(eif_corrupt, 1, var) / n_samp
    classic_tmle_fit[, "var_psi"] <- var_eif_corrupt
  }

  # HOW DOES TMLE COMPUTE P-VALUES?
  #pvalue <- 2 * pnorm(-abs(psi / sqrt(var_psi)))
  classic_tmle_out <- classic_tmle_fit %>%
    as_tibble() %>%
    transmute(
      gene_id = gene_id,
      ci_lwr = psi - sqrt(var_psi) * abs(stats::qnorm(p = (1 - 0.95) / 2)),
      psi = psi,
      ci_upr = psi + sqrt(var_psi) * abs(stats::qnorm(p = (1 - 0.95) / 2)),
      se_psi = sqrt(var_psi),
      p_val = case_when(stat_ref == "normal" ~ 2 *
                          pnorm(-abs(psi / sqrt(var_psi))),
                        stat_ref == "logistic" ~
                          2 * plogis(-abs(psi) / sqrt(var_psi),
                                     location = 0, scale = sqrt(3) / pi)
                       ),
      method = rep("classic_tmle", n()),
      label =
        if (is.null(no_findings)) {
          rep("null", nrow(data_se))
        } else {
          c(rep("null", nrow(data_se) - no_findings), rep("hit", no_findings))
        }
    )

  # testing limma
  limma_fit <- lmFit(assay(data_se), as.matrix(colData(data_se)))
  limma_fit <- eBayes(limma_fit)
  limma_test <- topTable(limma_fit, coef = 1, number = Inf, sort.by = "none",
                         adjust.method = adjust_method, confint = 0.95)
  limma_out <- limma_test %>%
    as_tibble() %>%
    transmute(
      gene_id = seq_len(n()),
      ci_lwr = CI.L,
      psi = logFC,
      ci_upr = CI.R,
      # how limma computes its multiplier for CI construction
      se_psi = sqrt(limma_fit$s2.post) * limma_fit$stdev.unscaled[, "A"],
      p_val = P.Value,
      method = rep("limma", n()),
      label =
        if (is.null(no_findings)) {
          rep("null", nrow(data_se))
        } else {
          c(rep("null", nrow(data_se) - no_findings), rep("hit", no_findings))
        }
    )

  # output
  sim_results <- rbind(biotmle_out, classic_tmle_out, limma_out)
  return(sim_results)
}
