generateQCPlots <- function(counts, metaUnfilter){
	counts <- counts[rowSums(counts)>0,]
	nGenes <- length(counts[,1])
	coverage <- colSums(counts)/nGenes
	write.table(data.frame(sampleId=colnames(counts), group=metaUnfilter$GroupCol, coverage=coverage), file='coveragePerFeature.txt', sep='\t', row.names = FALSE, quote=FALSE)
	
	barplot(coverage,xaxt='n',col=metaUnfilter$GroupCol, ylab="Counts per gene")

	nonZero <- colSums(counts > 0)
	barplot(nonZero,xaxt='n',col=metaUnfilter$GroupCol, ylab="Non-zero genes per sample")

	nonZero <- colSums(counts > 10)
	barplot(nonZero,xaxt='n',col=metaUnfilter$GroupCol, ylab="Genes >10 per sample")
}

