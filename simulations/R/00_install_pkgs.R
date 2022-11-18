if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/nhejazi/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# set CRAN mirror
options(repos = structure(c(CRAN = "https://cran.rstudio.com/")))
unixtools::set.tempdir("/global/scratch/nhejazi/rtmp")

# lie to pkgbuild, as per Jeremy
pkgbuild:::cache_set("has_compiler", TRUE)

# from CRAN
install.packages(c("R.utils", "Rcpp", "here", "remotes", "devtools",
                   "BiocManager", "tidyverse", "data.table", "Rsolnp", "nnls",
                   "foreach", "doRNG", "future", "future.apply", "doFuture",
                   "mvtnorm", "glmnet", "randomForest", "earth", "arm", "gam",
                   "nnet", "xgboost", "speedglm", "SuperLearner", "drtmle"),
                   lib = "/global/scratch/nhejazi/R")

# from Bioconductor
BiocManager::install(c("S4Vectors", "BiocGenerics", "BiocParallel",
                       "SummarizedExperiment", "limma", "biotmleData"),
                     lib = "/global/scratch/nhejazi/R",
                     update = FALSE)

# use remotes to install from GitHub
remotes::install_github(c("nhejazi/biotmle@master",
                          "tlverse/hal9001@master"),
                        lib = "/global/scratch/nhejazi/R")

# update all packages
#update.packages(ask = FALSE, lib.loc = "/global/scratch/nhejazi/R")
