################################################################################
# Apply biotmle to Subsampled Datasets
################################################################################

# load packages from renv library
renv::activate(here::here())

# load required packages and scripts
library(here)
library(dplyr)
library(readr)
library(limma)
library(minfi)
library(biotmle)
library(SuperLearner)
library(RhpcBLASctl)
source(here("R", "apply_biotmle.R"))
source(here("R", "subsample_biotmle.R"))
set.seed(7612452)

# explicitly load packages from SL library for renv awareness
library(earth)
library(glmnet)
library(ranger)

# disable multithreading
blas_set_num_threads(1L)
omp_set_num_threads(1L)

# load the prepped GenomicRatioSet
scratch_path <- "/global/scratch/users/philippe_boileau/su_et_al_data"
grs <- read_rds(here(scratch_path, "prepared-grs.Rds"))

# number of subsamples per sample size
n_subsamples <- 10 

# 25% (63 observations)
p25_subsample_ls <- lapply(
  seq_len(n_subsamples),
  function(idx) {subsample_biotmle(grs, sample_prop = 0.25)}
)

# 50% (127 observations)
p50_subsample_ls <- lapply(
  seq_len(n_subsamples),
  function(idx) {subsample_biotmle(grs, sample_prop = 0.50)}
)

# 75% (190 observations)
p75_subsample_ls <- lapply(
  seq_len(n_subsamples),
  function(idx) {subsample_biotmle(grs, sample_prop = 0.75)}
)

# save the results
write_rds(p25_subsample_ls, file = here("results", "p25-subsample-ls.Rds"))
write_rds(p50_subsample_ls, file = here("results", "p50-subsample-ls.Rds"))
write_rds(p75_subsample_ls, file = here("results", "p75-subsample-ls.Rds"))