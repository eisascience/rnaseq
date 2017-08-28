list.of.packages <- c('digest', 'plyr', 'perm', 'reshape', 'reshape2', 'knitr', 'ggplot2', 'Rcpp', 'lattice', 'tidyverse', 'stringr', 'edgeR', 'gplots', 'pheatmap', 'Rtsne', 'Rlabkey')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
	install.packages(new.packages, dependencies=TRUE, repos='http://cran.rstudio.com')
} else {
	print('no updates needed')
}

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
