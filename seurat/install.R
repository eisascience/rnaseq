install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()

toInstall <- c('DESeq2', 'destiny', 'MAST', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'dplyr', 'dtplyr', 'enrichR')
BiocManager::install(toInstall, update = TRUE, ask = FALSE)

devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0', dependencies = T, upgrade = 'always')