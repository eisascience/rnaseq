list.of.packages <- c('digest', 'plyr', 'perm', 'reshape', 'reshape2', 'knitr', 'ggplot2', 'Rcpp', 'lattice', 'tidyverse', 'stringr', 'edgeR', 'gplots', 'pheatmap', 'Rtsne', 'Rlabkey', 'devtools')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]
if(length(new.packages)){
	install.packages(new.packages, dependencies=TRUE, repos='http://cran.rstudio.com')
} else {
	print('no updates needed')
}

source('http://bioconductor.org/biocLite.R')
biocLite(c('RcppEigen', 'Cairo', 'rgl', 'limma', 'fry', 'goana', 'kegga', 'plotMD', 'plotWithHighlights', 'DESeq2', 'edgeR', 'glmQLFit', 'statmod', 'biomaRt', 'sva', 'scater', 'scran', 'destiny', 'WGCNA', ask=FALSE, dependencies=TRUE)
#biocLite('scde', ask=FALSE, dependencies=TRUE)

# see: 
# https://github.com/hms-dbmi/scde/issues/48
# https://groups.google.com/forum/#!topic/singlecellstats/rbFUTOQ9wu4
require(devtools)
uninstall.packages('flexmix')
install_version('flexmix', version = '2.3-13', repos = 'http://cran.us.r-project.org')
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE, force = TRUE)
