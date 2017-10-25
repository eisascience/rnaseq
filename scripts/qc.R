generateQCPlots <- function(counts, meta){
	counts <- counts[rowSums(counts)>0,]
	nGenes <- length(counts[,1])
	coverage <- colSums(counts)/nGenes
	write.table(data.frame(sampleId=colnames(counts), group=meta$GroupCol, coverage=coverage), file='coveragePerFeature.txt', sep='\t', row.names = FALSE, quote=FALSE)
	
	barplot(coverage,xaxt='n',col=meta$GroupCol, ylab="Counts per gene")

	nonZero <- colSums(counts > 0)
	barplot(nonZero,xaxt='n',col=meta$GroupCol, ylab="Non-zero genes per sample")

	nonZero <- colSums(counts > 10)
	barplot(nonZero,xaxt='n',col=meta$GroupCol, ylab="Genes >10 per sample")
	
	qplot(EstimatedLibrarySize, TotalNonZeroFeatures, data=meta, color=Peptide)
	
	qplot(EstimatedLibrarySize, TotalNonZeroFeatures, data=meta, color=AnimalId)
}

