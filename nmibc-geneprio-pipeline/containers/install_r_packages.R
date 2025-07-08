# install_r_packages.R
libdir <- "~/tools/r_env/lib/R"
dir.create(libdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(libdir)

install.packages(c(
  "remotes", "R.utils", "coloc", "susieR", "ggplot2", "viridis",
  "optparse", "data.table", "Rcpp", "RcppEigen", "dplyr"
), repos = "https://cloud.r-project.org")

remotes::install_github("gabraham/plink2R")
