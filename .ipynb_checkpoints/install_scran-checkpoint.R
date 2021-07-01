#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

##=== may solve some error ===##
Sys.setenv(R_INSTALL_STAGED = FALSE)

#BiocManager::install("Biobase")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("scran")


BiocManager::install(c("rsvd", "Rhdf5lib", "BiocNeighbors", "BiocSingular", "rhdf5", "HDF5Array", "DelayedMatrixStats", "scater", "scran"))