install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()

toInstall <- c('DESeq2', 'destiny', 'MAST', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'dplyr', 'dtplyr', 'enrichR', 'rtracklayer', 'data.table', 'naturalsort', 'Rlabkey')
BiocManager::install(toInstall, update = FALSE, ask = FALSE)

suppressWarnings(BiocManager::install(update=TRUE, ask=FALSE))

devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0', dependencies = T, upgrade = 'always')