require(devtools)

# see: 
# https://github.com/hms-dbmi/scde/issues/48
# https://groups.google.com/forum/#!topic/singlecellstats/rbFUTOQ9wu4
remove.packages('flexmix')
install_version('flexmix', version = '2.3-13', repos = 'http://cran.us.r-project.org')
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE, force = TRUE)

# Seurat
devtools::install_github('satijalab/seurat', ref = 'v2.1.0')

# Monocle
source('http://bioconductor.org/biocLite.R')
biocLite()
biocLite('monocle', 'scater')