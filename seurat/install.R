install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')
devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0', dependencies = T, upgrade = 'always')