source('https://raw.github.com/bbimber/rnaseq/master/scripts/main.R')

args = commandArgs(trailingOnly=TRUE)

setwd('/work')
register(MulticoreParam(12))
Sys.setenv('ALLOW_WGCNA_THREADS' = 12)

geneCountTableFile=args[1]
groupColName = args[2]

allowableIdsTable <- read.table('Metadata.txt', header = TRUE, sep = '\t')

l <- prepareTables(allowableIdsTable$ReadsetId, geneCountTableFile, groupColName, minLibrarySize=500000)
meta <- l$meta
geneCountMatrix <- l$geneCounts

if(!exists('designF')){
	designF = ' ~ AnimalId + GroupCol'
}

designFormula <- as.formula(designF)
coefs <- length(colnames(design))
contrast <- c(groupColName, levels(meta$GroupCol)[nlevels(meta$GroupCol)], levels(meta$GroupCol)[1])

rmarkdown::render('RNASeq.rmd', clean=TRUE)
