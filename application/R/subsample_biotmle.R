################################################################################
# Subsample from full dataset, and apply biotmle
################################################################################

#' Subsample the GenomicRatioSet and Apply the biotmle Pipeline
#' 
#' This function subsamples the full dataset according to the desired subsample
#' size, and then applies the biotmle pipeline to the reduced data.
#' 
#' @param grs The full dataset GenomicRatioSet.
#' @param sample_prop A numeric between 0 and 1 dictating the size of the
#'   subsample.
#' @param dist_choice A character indicating the distribution to be used when
#'   computing p-values. Defaults to "normal".
#'
#' @return The biotmle results table of the subsampled data, augmented with the
#'   specially computed FWER adjusted p-values using the method of Tuglus and
#'   van der Laan.
subsample_biotmle <- function(grs, sample_prop, dist_choice = "normal") {
  
  # sample without replacement from the full dataset
  n_samp <- ncol(grs)
  sample_idx <- seq_len(n_samp)
  grs_sub <- grs[, sample(sample_idx, size = round(sample_prop * n_samp))]
  
  # prepare the design matrix
  design_mat <- pData(grs_sub) %>%
    as.data.frame() %>%
    select(smoking_status, sex, age, race, CD8T, CD4T, NK, Bcell, Mono, Gran) %>%
    mutate(
      smoking_status = if_else(smoking_status == "never", 0, 1),
      sex = if_else(sex == "male", 0, 1),
      race = if_else(race == "white", 0, 1)
    )
  
  # apply biotmle
  biotmle_results <- apply_biotmle(
    grs = grs_sub,
    design_mat = design_mat,
    pval_cutoff = 0.05,
    dist_choice = dist_choice
  )
  
  # return the topTable output
  return(biotmle_results)
  
}