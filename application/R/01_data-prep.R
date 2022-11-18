################################################################################
# Prepare Su et al. Data #######################################################
################################################################################

# Data pertaining to the epigenetic effects of tobacco smoking on whole blood
# from Su et al. is prepared for analysis in this script. See
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166486 for
# more background. The GEO series number is GSE85210.

# load packages from renv library
renv::activate(project = here::here())

# load the required packages
library(here)
library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GEOquery)

# path to scratch directory
scratch_path <- "/global/scratch/users/philippe_boileau/su_et_al_data/"

# download the normalized beta values
# scp'ed it from local computer to scratch due to download timeout
#getGEOSuppFiles(
  #"GSE85210",
  #filter_regex = "*Matrix_processed*",
  #baseDir = scratch_path
#)
betas <- gzfile(
  description = here(
    scratch_path, "GSE85210", "GSE85210_Matrix_processed.txt.gz"
  )
) %>% read_table(col_names = FALSE, skip = 2)

# extract the colnames (they are poorly formatted, causing issues)
colnames_betas <- gzfile(
  description = here(
    scratch_path, "GSE85210", "GSE85210_Matrix_processed.txt.gz"
  )
) %>%
  read_tsv() %>%
  colnames()
colnames(betas) <- colnames_betas
betas <- betas[, -507]
betas <- betas %>% dplyr::select(!contains("Detection Pval"))
rownames_betas <- betas$ID_REF
betas <- betas[, -1] %>% as.matrix()
rownames(betas) <- rownames_betas

# load the separately provided metadata (cell type composition, covariates)
su_metadata <- read_xlsx(
  path = here("data", "se-smokers-metadata-for-phillipe.xlsx")
) %>%
  mutate(
    sample_name = Sample_Name1,
    smoking_status = if_else(
      `Sample_Group Never=0 Ever=1` == 0, "never", "ever"
    ),
    smoking_status = factor(smoking_status, levels = c("never", "ever")),
    sex = if_else(`sex M=1 F=2` == 1, "male", "female"),
    sex = factor(sex, levels = c("male", "female")),
    race = if_else(`race W=1, B=2` == 1, "white", "black"),
    race = factor(race, levels = c("white", "black")),
    pack_years = `pack-years`
  ) %>%
  select(
    sample_name, sex, race, age, smoking_status, pack_years,
    CD8T, CD4T, NK, Bcell, Mono, Gran
  ) %>%
  as.data.frame()
rownames(su_metadata) <- su_metadata$sample_name

# extract the metadata from GEO
# commented out on slurm due to download timeout -- only usedto check for
# sample data missmatching, which is confirmed not to be an issue when run
# locally

# Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
# geo_metadata <- getGEO("GSE85210")$GSE85210_series_matrix.txt.gz

# check for sample missmatching
# geo_pheno_data <- pData(geo_metadata) %>%
#   as.data.frame() %>%
#   mutate(sample_name = str_extract(title, "SE(.*)"))
# meta_sample_order <- order(geo_pheno_data$sample_name, su_metadata$sample_name)
# geo_pheno_data <- geo_pheno_data[meta_sample_order, ]
# stopifnot(
#   identical(if_else(
#     geo_pheno_data$`subject status:ch1` == "smoker", "ever", "never"),
#     as.character(su_metadata$smoking_status)
#   )
# )

# create a single genomic ratio set
betas_order <- order(colnames(betas), su_metadata$sample_name)
betas <- betas[, betas_order]
grs <- makeGenomicRatioSetFromMatrix(
  betas, rownames = rownames(betas), pData = su_metadata,
  array = "IlluminaHumanMethylation450k", what = "Beta"
)

# save the genomic ratio set
write_rds(grs, file = here(scratch_path, "prepared-grs.Rds"))
