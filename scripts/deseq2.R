prepareDDS <- function(geneCountMatrix, sizeFactorsFromEdgeR=FALSE){
	countdata=data.frame(geneCountMatrix)
	
	#not necessary if columns have proper names
	#names(countdata) <- gsub('X', '', names(countdata))

	coldata=metaUnfilter 
	rownames(coldata)=names(countdata)

	dds <- DESeqDataSetFromMatrix(countData = countdata,colData = coldata,design = designFormula)
	sum(rowSums(counts(dds)) ==0)  #this dataset does have all 0 count genes
	dds <- dds[ rowSums(counts(dds)) > 1, ]

	if (sizeFactorsFromEdgeR == TRUE){
		#Note: because these samples errored due to all genes having at least one zero value, calculate norm factors w/ edgeR
		sizeFactors(dds) <- calcNormFactors(counts(dds))

		#or: https://support.bioconductor.org/p/63229/
		#ddsCounts<-counts(dds)
		#geoMeans = apply(ddsCounts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
		#dds<-estimateSizeFactors(dds, geoMeans=geoMeans)
	}
	
	dds <- DESeq(dds, parallel=TRUE)
	
	return(dds)
}

doGenerateDeseq2Summary <- function(dds){
	plotSparsity(dds)

	counts <- counts(dds, normalized=TRUE)
	counts <- counts[rowSums(counts)>0,]
	nGenes <- length(counts[,1])
	plot(sizeFactors(dds),colSums(counts)/nGenes)

	DESeq2result <- transform(DESeq2result, Ensembl = colsplit(GeneID, split = "\\|", names = c('Id', 'Name')))
	DESeq2result$Ensembl.Name[as.character(DESeq2result$Ensembl.Name) == as.character(DESeq2result$Ensembl.Id)] <- c(NA)
	write.table(DESeq2result,'DESeq2result.txt', sep='\t', quote=FALSE, row.names=FALSE)

	topGenes <- DESeq2result[DESeq2result$padj < 0.05,]
	topGeneIds <- topGenes$GeneID
	topGenes <- arrange(topGenes, padj)
	if (nrow(topGenes) > 0){
		topGenes <- addEnsembl(topGenes)
	}
	write.table(topGenes ,'DESeq2result_top.txt', sep='\t', quote=FALSE, row.names=FALSE)

	#histogram of p-values
	hist(DESeq2result$pvalue)

	return(topGenes)
}