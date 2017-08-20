list.of.packages <- c('digest', 'plyr', 'perm', 'reshape', 'reshape2', 'knitr', 'ggplot2', 'Rcpp', 'lattice', 'tidyverse', 'stringr', 'edgeR', 'gplots', 'pheatmap', 'WGCNA')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
	install.packages(new.packages, dependencies=TRUE, repos='http://cran.rstudio.com')
} else {
	print('no updates needed')
}

source("http://bioconductor.org/biocLite.R")
biocLite("rgl", ask=FALSE)
biocLite("limma", ask=FALSE)
biocLite(c('fry', 'goana', 'kegga', 'plotMD', 'plotWithHighlights'), ask=FALSE)
biocLite("DESeq2", ask=FALSE)
biocLite("edgeR", ask=FALSE)
biocLite("glmQLFit", ask=FALSE)
biocLite("statmod", ask=FALSE)
biocLite("biomaRt", ask=FALSE)
biocLite("sva", ask=FALSE)
biocLite("scater", ask=FALSE)
biocLite("scran", ask=FALSE)
