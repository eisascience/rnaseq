source('https://raw.github.com/bbimber/rnaseq/master/scripts/main.R')

args = commandArgs(trailingOnly=TRUE)

setwd('/work')
register(MulticoreParam(12))

geneCountTableFile=args[1]
groupColName = args[2]

allowableIdsTable <- read.table('Metadata.txt', header = TRUE, sep = '\t')

metaUnfilter <- prepareMetadataTable(fullMetadata, allowableIdsTable$OutputFileId, groupColName)
geneCountMatrix <- parseGeneCountsMatrix(geneCountTableFile, as.character(metaUnfilter$OutputFileId))
metaUnfilter <- prepareMetadataTable2(metaUnfilter, geneCountMatrix)

if(!exists('designF'){
	designF = ' ~ AnimalId + GroupCol'
}

designFormula <- as.formula(designF)
coefs <- length(colnames(design))
contrast <- c(groupColName, levels(metaUnfilter$GroupCol)[nlevels(metaUnfilter$GroupCol)], levels(metaUnfilter$GroupCol)[1])

rmarkdown::render('RNASeq.rmd', clean=TRUE)