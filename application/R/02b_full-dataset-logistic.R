################################################################################
# Apply biotmle to Full Dataset ################################################
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
set.seed(4512316)

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

# create the design matrix
design_mat <- pData(grs) %>%
  as.data.frame() %>%
  select(smoking_status, sex, age, race, CD8T, CD4T, NK, Bcell, Mono, Gran) %>%
  mutate(
    smoking_status = if_else(smoking_status == "never", 0, 1),
    sex = if_else(sex == "male", 0, 1),
    race = if_else(race == "white", 0, 1)
  )

# apply biotmle
biotmle_results <- apply_biotmle(
  grs = grs,
  design_mat = design_mat,
  pval_cutoff = 0.05,
  dist_choice = "logistic"
)

# save the results
saveRDS(biotmle_results, file = here("results", "full-dataset-biotmle-logistic.Rds"))
