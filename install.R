list.of.packages <- c('digest', 'plyr', 'perm', 'reshape', 'reshape2', 'knitr', 'ggplot2', 'Rcpp', 'lattice', 'tidyverse', 'stringr', 'edgeR', 'gplots')
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
biocLite("DESeq2", ask=TRUE)
biocLite("edgeR", ask=TRUE)
biocLite("glmQLFit", ask=TRUE)
biocLite("statmod", ask=TRUE)
