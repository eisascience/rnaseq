list.of.packages <- c('digest', 'plyr', 'perm', 'reshape', 'reshape2', 'knitr', 'ggplot2', 'Rcpp', 'lattice', 'tidyverse', 'stringr', 'gplots', 'pheatmap', 'Rtsne', 'Rlabkey', 'devtools', 'codetools', 'RcppEigen', 'Cairo', 'rgl', 'limma', 'fry', 'goana', 'kegga', 'plotMD', 'plotWithHighlights', 'DESeq2', 'edgeR', 'glmQLFit', 'statmod', 'biomaRt', 'sva', 'scater', 'scran', 'destiny', 'WGCNA')

source('http://bioconductor.org/biocLite.R')
biocLite(list.of.packages, ask=FALSE, dependencies=TRUE)

# see: 
# https://github.com/hms-dbmi/scde/issues/48
# https://groups.google.com/forum/#!topic/singlecellstats/rbFUTOQ9wu4
require(devtools)
remove.packages('flexmix')
install_version('flexmix', version = '2.3-13', repos = 'http://cran.us.r-project.org')
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE, force = TRUE)

devtools::install_github("satijalab/seurat", ref = "v2.1.0")