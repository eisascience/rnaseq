list.of.packages <- c('digest', 'plyr', 'perm', 'reshape', 'reshape2', 'knitr', 'ggplot2', 'Rcpp', 'lattice', 'tidyverse', 'stringr', 'edgeR', 'gplots', 'pheatmap', 'Rtsne', 'Rlabkey', 'devtools')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
	install.packages(new.packages, dependencies=TRUE, repos='http://cran.rstudio.com')
} else {
	print('no updates needed')
}

# see: https://groups.google.com/forum/#!topic/singlecellstats/rbFUTOQ9wu4
require(devtools)
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")

source("http://bioconductor.org/biocLite.R")
biocLite("RcppEigen", ask=FALSE, dependencies=TRUE)
biocLite("Cairo", ask=FALSE, dependencies=TRUE)
biocLite("rgl", ask=FALSE, dependencies=TRUE)
biocLite("limma", ask=FALSE, dependencies=TRUE)
biocLite(c('fry', 'goana', 'kegga', 'plotMD', 'plotWithHighlights'), ask=FALSE, dependencies=TRUE)
biocLite("DESeq2", ask=FALSE, dependencies=TRUE)
biocLite("edgeR", ask=FALSE, dependencies=TRUE)
biocLite("glmQLFit", ask=FALSE, dependencies=TRUE)
biocLite("statmod", ask=FALSE, dependencies=TRUE)
biocLite("biomaRt", ask=FALSE, dependencies=TRUE)
biocLite("sva", ask=FALSE, dependencies=TRUE)
biocLite("scater", ask=FALSE, dependencies=TRUE)
biocLite("scran", ask=FALSE, dependencies=TRUE)
biocLite('destiny', ask=FALSE, dependencies=TRUE)
biocLite('scde', ask=FALSE, dependencies=TRUE)
biocLite('WGCNA', ask=FALSE, dependencies=TRUE)

# see: https://github.com/hms-dbmi/scde/issues/48
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)
