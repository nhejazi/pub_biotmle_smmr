################################################################################
# biotmle and related helper functions
################################################################################

# FWER Correction a la Tuglus et van der Laan ##################################
# Inspired by: https://github.com/nhejazi/methyvim/blob/master/R/fdr_msa.R

#' FWER-MSA correction
#'
#' Modified FWER Controlling Procedure for Multi-Stage Analyses (orignally for
#' FDR control, from MJ van der Laan and C Tuglus, 2009,
#' <doi:10.2202/1544-6115.1397>)
#'
#' @param pvals Numeric vector containing the p-values that result from any
#'  chosen statistical hypothesis testing procedure.
#' @param total_cpgs Numeric indicating the total number of observations that
#'  would have been available for testing prior to the selection procedure
#'  employed in the multi-stage analysis performed.
#'
#' @return A \code{numeric} vector of corrected p-values, controlling the False
#'  Discovery Rate, using the method of Tuglus and van der Laan.
#'
#' @examples
#' g <- 1e4
#' n <- 1e2
#' p <- abs(rnorm(n, mean = 1e-8, sd = 1e-2))
#' # treating the vector p as one of p-values, FDR-MSA may be applied
#' fwer_p <- fwer_msa(pvals = p, total_obs = g)
fwer_msa <- function(pvals, total_cpgs) {
  pvals_not_tested <- rep(1, total_cpgs - length(pvals))
  pvals_all <- c(pvals, pvals_not_tested)
  fwer_adj <- stats::p.adjust(pvals_all, method = "holm")
  fwer_out <- fwer_adj[seq_along(pvals)]
  return(fwer_out)
}


# Apply limma for initial filtering ############################################

#' Apply limma's variance-moderated robust linear regression methodology to the
#' CpG sites of a GenomicRationSet. The p-values resulting from the differential
#' methylation tests are then used to filter the CpGs for further testing.
#' 
#' @param grs A GenomicRatioSet object.
#' @param design_mat A design matrix. Make sure that the coefficient of interest
#'   is located in the second column immediately following the intercept.
#' @param pval_cutoff A numeric indicating the FWER-adjusted p-value cutoff used
#'   to filter the CpG sites.
#'
#' @return A character vector of CpGs whose moderated t-tests' nominal p-values
#'   are greater than the p-value cutoff.
limma_filter <- function(grs, design_mat, pval_cutoff) {
  # transform the design_mat into an actual design matrix
  design_mat <- model.matrix(
    ~ smoking_status + sex + age + race + CD8T + CD4T + NK + Bcell + Mono + Gran,
    data = design_mat
  )

  # extract M values
  mvals <- minfi::getM(grs)

  # fit the linear model
  lm_fit <- limma::lmFit(object = mvals, design = design_mat)

  # fit the empirical Bayes model
  e_fit <- limma::eBayes(lm_fit, robust = TRUE)

  # create table of results
  p_top <- limma::topTable(
    e_fit, adjust.method = "holm", coef = "smoking_status",
    number = Inf, p.value = pval_cutoff, sort.by = "p"
  )

  # return the CpGs satisfying the cutoff
  cpgs <- rownames(p_top)
  return(cpgs)
}


# Apply biotmle with basic library of learners #################################

#' A wrapper function for the biotmle methodology. Simply provide a
#' GenomicRatioSet, a design matrix, and a p-value cutoff in order to identify
#' CpG sites that predicted to be differentially methylated.
#' 
#' @param grs A GenomicRatioSet object.
#' @param design_mat A design matrix. Make sure that the coefficient of interest
#'   is located in the second column immediately following the intercept.
#' @param pval_cutoff A numeric indicating the p-value cutoff used to filter the
#'   CpG sites.
#' @param dist_choice A character indicating the choice of p value type to use.
#' @param known_cpgs A character vector of CpGs believed a priori to be
#'   associated with the exposure of interest. These CpGs will not be filtered
#'   during the intial filtering step. Defaults to NULL.
#'
#' @return The biotmle results table, augmented with the specially computed
#'   FWER adjusted p-values using the method of Tuglus and van der Laan.
apply_biotmle <- function(
  grs,
  design_mat,
  pval_cutoff,
  dist_choice = "normal",
  known_cpgs = NULL) {
  # get the number of CpGs
  num_cpgs <- nrow(grs)

  # filter the GenomicRatioSet using the limma filtering strategy
  keep_cpgs <- limma_filter(grs, design_mat, pval_cutoff)
  keep_cpgs <- dplyr::union(keep_cpgs, known_cpgs)
  grs <- grs[keep_cpgs, ]

  # turn GenomicRatioSet into a SummarizedExperiment
  se <- SummarizedExperiment(
    assays = list(counts = minfi::getM(grs)),
    colData = design_mat
  )

  # estimate ATEs
  biotmle_res <- biomarkertmle(
    se = se,
    varInt = 1,
    g_lib = c("SL.mean", "SL.glm", "SL.glmnet", "SL.earth"),
    Q_lib = c("SL.mean", "SL.glm", "SL.glmnet", "SL.earth", "SL.ranger"),
    cv_folds = 2
  )

  # apply variance moderation, and generate the topTable object
  # NOTE: `sort.by = None` and `number = Inf` are hardcoded in `modtest_ic`
  # TODO: should we check `pval_type="logistic"`? we recommend this in the
  #       manuscript, but empirically (from running a subset of CpGs in the
  #       data analysis), it makes an order of magnitude difference in the raw
  #       p-value, so it may be too conservative a correction when used in
  #       combination with filtering and FDR-MSA
  biotmle_res <- modtest_ic(
    biotmle = biotmle_res, adjust = "none", pval_type = dist_choice
  )
  p_top <- biotmle_res@topTable
  p_top$adj.P.Val <- NULL # just the nominal p-value

  # compute BH adjusted p-values using method of Tuglus and van der Laan
  p_top$adj_pval <- fwer_msa(p_top$P.Value, num_cpgs)
  return(p_top)
}
