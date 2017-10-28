runEdgeRQL <- function(geneCountMatrix, geneAnnotations=NULL){
	y = DGEList(geneCountMatrix,genes=geneAnnotations) 
	keep = rowSums(cpm(y)>1)>=3 #filtering

	y = y[keep,,keep.lib.sizes=FALSE]
	y = calcNormFactors(y) ##TMM normalization

	#this initial dispersion estimation is need for glmQLFit dispersion squeeze
	#this function appends the original object
	y_QL = estimateDisp(y, design, robust=TRUE) 

	return(y_QL)
}

writeEdgeRSummary <- function(qlf2, y_QL, outputFile){
	QLresult = topTags(qlf2,n=nrow(y_QL)) 
	QLresult <- as.data.frame(QLresult)
	#write.table(QLresult, file='edgeR_QLresult.txt')
	
	QLresult$GeneID <- rownames(QLresult)
	QLresult = QLresult %>% dplyr::select(GeneID,logFC,logCPM,PValue,FDR) %>% mutate(dir=sign(logFC)*(FDR<0.05))%>% arrange(GeneID)
	row.names(QLresult)=QLresult$GeneID

	QLresult$Ensembl.Id <- QLresult$GeneID
	write.table(QLresult,outputFile, quote=FALSE, sep='\t', row.names=FALSE)
	
	return(QLresult)
}

writeEdgeRTopGenes <- function(qlf2, y_QL, pval = 0.05, outputFile){
	QLresultTop = topTags(qlf2, p.value=pval) 
	QLresultTop <- as.data.frame(QLresultTop)
	if (nrow(QLresultTop) == 0){
		return()
	}

	QLresultTop$GeneID <- rownames(QLresultTop)
	QLresultTop = QLresultTop %>% dplyr::select(GeneID,logFC,logCPM,PValue,FDR) %>% mutate(dir=sign(logFC)*(FDR<0.05)) %>% arrange(PValue)

	QLresultTop$Ensembl.Id <- QLresultTop$GeneID
	QLresultTop <- addEnsembl(QLresultTop)
	write.table(QLresultTop, outputFile, quote=FALSE, sep='\t', row.names=FALSE)
	
	return(QLresultTop)
}

doPlotMDS <- function(y_QL, meta){
	#for mds:
	for (col in c('AnimalId', 'Peptide', 'Population', 'Activated', 'NumCDR3s', 'CellClass', 'Treatment', 'DistinctLoci')){
		plotMDS.DGEList(y_QL, labels = meta[[col]])
	}
}

generateEdgeRSummary <- function(y_QL, qlf2, meta, suffix){
	#histogram of p-values
	hist(qlf2$table$PValue)

	plotBCV(y_QL)

	#generate outputs
	QLresult <- writeEdgeRSummary(qlf2, y_QL, outputFile=paste0('edgeR_', suffix, '_allgenes.txt'))

	#top genes
	QLresultTop <- writeEdgeRTopGenes(qlf2, y_QL, outputFile=paste0('edgeR_',suffix,'_topGenes.txt'))

	doPlotMDS(y_QL, meta)

	summary(de <- decideTestsDGE(qlf2, p.value=0.0000001, lfc=3))
	detags <- rownames(y_QL)[as.logical(de)]
	print(length(detags))

	plotSmear(qlf2, de.tags=detags)
	abline(h=c(-1, 1), col="blue")

	if (length(detags) > 1){
	  logcpm <- cpm(y_QL[detags,], prior.count=0.5, log=TRUE)
	  my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
	  heatmap.2(logcpm , labCol=meta$GroupCol, col=my_palette)
	}
	
	return(QLresultTop)
}
