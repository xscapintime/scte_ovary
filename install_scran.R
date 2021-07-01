#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

##=== may solve some error ===##
Sys.setenv(R_INSTALL_STAGED = FALSE)

#BiocManager::install("Biobase")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("scran")

## Bioconductor package
# set Bioc to R 3.6 version
BiocManager::install(verison = "3.10")

# install
#BiocManager::install("Rhdf5lib", verison = "3.10", INSTALL_opts = c('--no-lock')) # failed everytime, and finally installed by conda


# reinstalling RcppAnnoy v0.0.14 keeps BiocNeighbors working
# https://github.com/LTLA/BiocNeighbors/issues/10
utils::install.packages(
    pkgs = paste(
        "https://cran.r-project.org",
        "src",
        "contrib",
        "Archive",
        "RcppAnnoy",
        "RcppAnnoy_0.0.14.tar.gz",
        sep = "/"
    ),
    repos = NULL,
    type = "source"
)

BiocManager::install("DelayedMatrixStats", verison = "3.10", INSTALL_opts = c('--no-lock'))
BiocManager::install("HDF5Array", verison = "3.10", INSTALL_opts = c('--no-lock'))
BiocManager::install("BiocNeighbors", verison = "3.10", INSTALL_opts = c('--no-lock'))
### FINALLY ###
BiocManager::install("scran", verison = "3.10", INSTALL_opts = c('--no-lock'))


# CRAN package
# packages not available (for R version 3.6.1)

install.packages("https://cran.r-project.org/src/contrib/Archive/rsvd/rsvd_1.0.1.tar.gz") # older version


#BiocManager::install(c("rsvd", "Rhdf5lib", "BiocNeighbors", "BiocSingular", "rhdf5", "HDF5Array", "DelayedMatrixStats", "scater", "scran"))
