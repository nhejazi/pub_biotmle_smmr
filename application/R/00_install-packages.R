################################################################################
# Install Required Packages ####################################################
################################################################################

# In this script, the required packages are installed from CRAN and
# Bioconductor.

# set user-specific package library
.libPaths("/global/scratch/users/philippe_boileau/R")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
lib_loc <- "/global/scratch/users/philippe_boileau/R"

# set CRAN mirror
options(repos = structure(c(CRAN = "https://cran.rstudio.com/")))

# CRAN packages:
install.packages(
  c("BiocManager", "here", "dplyr", "readr", "readxl", "stringr", "speedglm",
    "arm", "ranger", "gam", "glmnet", "earth", "SuperLearner", "remotes"),
  lib = lib_loc
)

# Bioconductor packages:
BiocManager::install(
  c("minfi", "limma", "IlluminaHumanMethylation450kanno.ilmn12.hg19"),
  lib = lib_loc
)

# GitHub packages
remotes::install_github("nhejazi/biotmle")

# update all packages
update.packages(ask = FALSE, lib.loc = lib_loc)
