library(limma)
library(lattice)
library(edgeR)
library(tidyverse)
library(DESeq2)
library(stringr)
library(gplots)
library(ggplot2)
library(reshape)
library(biomaRt)
library(dplyr)
library("BiocParallel")
library(Rlabkey)
library(Matrix)

parseGeneCountsMatrix <- function(fileName, allowableColumnIDs = NA, filterZeroCounts=TRUE){
  geneCounts <- read.table(fileName, sep = '\t', header=TRUE, row.names = 1)
  geneCounts <- geneCounts[, !(names(geneCounts) %in% c('GeneId', 'GeneDescription', 'SamplesWithReads', 'GeneName')) ] 
  geneCounts <- geneCounts[!(rownames(geneCounts) %in% c('N_ambiguous', 'N_multimapping', 'N_noFeature', 'N_unmapped')),] 
  

  origCols <- ncol(geneCounts)
  
  colnames(geneCounts) <- gsub('X', '', colnames(geneCounts))
  geneCounts <- geneCounts[, sort(colnames(geneCounts))]
  
  if (!is.na(allowableColumnIDs)){
    toSelect <- intersect(allowableColumnIDs, colnames(geneCounts))
    geneCounts <- geneCounts[ toSelect ] 
    print(paste0('columns filtered from original gene count matrix because they were not in metadata: ', (origCols - ncol(geneCounts))))
  }
  
  geneCountMatrix <- data.matrix(geneCounts)
  
  if (filterZeroCounts){
    origRow <- nrow(geneCountMatrix)
    geneCountMatrix <- geneCountMatrix[ ,colSums(geneCountMatrix) > 0 ] 
    print(paste0('rows filtered from original gene count matrix because of zero gene feature counts: ', (origRow - nrow(geneCountMatrix))))
  }
  
  print(paste0('total samples: ', ncol(geneCountMatrix)))
               
  return(geneCountMatrix)
}

addClusterDifferentiation <- function(df, cdFile=NA, yField='GeneSymbol', xField='external_gene_name'){
  ref <- read.table(cdFile, sep='\t', quote='"', header=TRUE)
  
  return(merge(df, ref, by.x=c(xField), by.y=c(yField), all.x=TRUE))
}

addEnsembl <- function(df, geneField = 'Ensembl.Id'){
  ensembl = useEnsembl(biomart="ensembl", dataset="mmulatta_gene_ensembl")
  ensemblIds = unique(df[[geneField]])
  attrs <- unique(getBM(attributes=c('ensembl_gene_id','external_gene_name','description','name_1006', 'hgnc_symbol'), filters ='ensembl_gene_id', values=ensemblIds, mart = ensembl))
  
  attrs <- attrs %>%
    group_by(ensembl_gene_id, external_gene_name, hgnc_symbol, description) %>%
    summarise(goAnnotations = toString(sort(unique(name_1006))))
  attrs <- data.frame(attrs)
  
  return(merge(df, attrs, by.x=geneField, by.y='ensembl_gene_id'))
}

prepareTables <- function(allowableIds, geneCountTableFile, grouping, minLibrarySize=500000, minReads=500000, minFeatures=1000){
  metaUnfilter <- prepareMetadataTable(allowableIds, grouping, 'ReadsetId')
  geneCountMatrix <- parseGeneCountsMatrix(geneCountTableFile, as.character(metaUnfilter$ReadsetId))
  
  return(prepareMetadataTable2(metaUnfilter, geneCountMatrix, minLibrarySize = minLibrarySize, minReads = minReads, minFeatures = minFeatures))
}

prepareMetadataTable <- function(allowableIds, grouping, metadataFieldName='ReadsetId', groupingDataframe=NULL){
	#pull full metadata from LK
	metadataFull <- pullTCRMetaFromLabKey()
	metaUnfilter <- metadataFull[metadataFull[[metadataFieldName]] %in% allowableIds,]

	rownames(metaUnfilter) <- as.character(metaUnfilter[[metadataFieldName]])
	metaUnfilter$AnimalId <- as.factor(metaUnfilter$AnimalId)
	metaUnfilter$ReadsetId <- as.factor(metaUnfilter$ReadsetId)
	metaUnfilter$HasCDR3s <- c(FALSE)
	metaUnfilter$HasCDR3s[metaUnfilter$NumCDR3s > 0] <- c(TRUE)

	#metaUnfilter$EstimatedLibrarySize <- as.numeric(gsub(',','',metaUnfilter$EstimatedLibrarySize))
	if (is.character(grouping)){
	  metaUnfilter$GroupCol <- as.factor(metaUnfilter[[grouping]])  
	} else if(is.data.frame(grouping)) {
	  toJoin <- data.frame(ReadsetId=grouping$ReadsetId, GroupCol=grouping$GroupCol)
	  metaUnfilter <- merge(metaUnfilter, toJoin, by.x=c('ReadsetId'))
	}

	return(metaUnfilter)
}

prepareMetadataTable2 <- function(metaUnfilter, geneCountMatrix, minLibrarySize = 2000000, minFeatures = 1000, metadataFieldName='ReadsetId', minReads=500000){
  #first, limit to only those in the gene table, in case we dropped metadata rows
  metaUnfilter <- metaUnfilter[metaUnfilter[[metadataFieldName]] %in% colnames(geneCountMatrix),]
  
  metaUnfilter$TotalNonZeroFeatures <- colSums(geneCountMatrix != 0)
  
  #origRows <- nrow(metaUnfilter)
  #metaUnfilter <- metaUnfilter[metaUnfilter$EstimatedLibrarySize > minLibrarySize,]
  #print(paste0('rows dropped due to library size: ', (origRows - nrow(metaUnfilter))))
  
  origRows <- nrow(metaUnfilter)
  metaUnfilter <- metaUnfilter[is.na(metaUnfilter$Status),]	
  print(paste0('rows dropped due to non-blank status: ', (origRows - nrow(metaUnfilter))))
  
	origRows <- nrow(metaUnfilter)
	metaUnfilter <- metaUnfilter[metaUnfilter$TotalNonZeroFeatures > minFeatures,]
	print(paste0('rows dropped due to low non-zero features: ', (origRows - nrow(metaUnfilter))))
	
	origRows <- nrow(metaUnfilter)
	metaUnfilter <- metaUnfilter[metaUnfilter$TotalForwardReads >= minReads,]	
	print(paste0('rows dropped due to low reads: ', (origRows - nrow(metaUnfilter))))

	metaUnfilter$GroupCol <- as.factor(as.character(metaUnfilter$GroupCol))
	
	#also update gene table for those dropped rows
	origCols <- ncol(geneCountMatrix)
	toSelect <- intersect(metaUnfilter[[metadataFieldName]], colnames(geneCountMatrix))
	geneCountMatrix <- geneCountMatrix[ ,toSelect ] 
	
	if (ncol(geneCountMatrix) != nrow(metaUnfilter)){
	  stop('Rowcount of metadata does not match gene count!')
	}
	
	print(paste0('gene count cols dropped due to metadata drops: ', (origCols - ncol(geneCountMatrix))))
	
	return(list(meta=metaUnfilter, geneCounts=geneCountMatrix))
}

pullTCRMetaFromLabKey <- function(requireOuputFileId = FALSE, replicateAsSuffix = FALSE){
	df <- labkey.selectRows(
		baseUrl="https://prime-seq.ohsu.edu", 
		folderPath="/Labs/Bimber/", 
		schemaName="tcrdb", 
		queryName="cdnas", 
		viewName="", 
		colSelect=c(
		  'readsetId','readsetId/status','readsetId/workbook','readsetId/totalForwardReads','readsetId/numCDR3s','readsetId/distinctLoci','readsetId/numTcrRuns',
		  'sortId', 'sortId/stimId','sortId/stimId/animalId','sortId/stimId/date','sortId/stimId/stim','sortId/stimId/treatment','sortId/stimId/activated','sortId/stimId/background',
		  'sortId/population','sortId/replicate','sortId/cells', "rowid", "sortid/stimid/stim/category", "sortid/stimid/stim/type"
		),
		containerFilter=NULL,
		colNameOpt='rname'
	)

	#print(str(df))
	
	names(df)[names(df)=="readsetid"] <- "ReadsetId"
	df$ReadsetId <- as.integer(df$ReadsetId)
	
	names(df)[names(df)=="readsetid_totalforwardreads"] <- "TotalForwardReads"
	df$TotalForwardReads <- as.integer(df$TotalForwardReads)
	
	names(df)[names(df)=="readsetid_numcdr3s"] <- "NumCDR3s"
	df$NumCDR3s <- as.integer(df$NumCDR3s)
	
	names(df)[names(df)=="readsetid_distinctloci"] <- "DistinctLoci"
	df$DistinctLoci <- as.factor(df$DistinctLoci)
	
	names(df)[names(df)=="readsetid_status"] <- "Status"
	df$Status <- as.factor(df$Status)
	
	names(df)[names(df)=="sortid_stimid_animalid"] <- "AnimalId"
	df$AnimalId <- as.factor(df$AnimalId)
	
	names(df)[names(df)=="sortid_stimid_date"] <- "SampleDate"
	names(df)[names(df)=="rowid"] <- "DatasetId"
	
	names(df)[names(df)=="sortid_cells"] <- "Cells"
	df$Cells <- as.integer(df$Cells)
	
	names(df)[names(df)=="sortid_replicate"] <- "Replicate"
	df$Replicate <- as.factor(df$Replicate)
	
	names(df)[names(df)=="sortid"] <- "SortId"
	df$SortId <- as.factor(df$SortId)
	
	df$IsSingleCell <- df['Cells'] == 1
	
	#TODO: consider if this is the best approach
	df$Replicate[df$IsSingleCell] <- c(NA)
	
	names(df)[names(df)=="sortid_stimid_treatment"] <- "Treatment"
	df$Treatment <- as.factor(df$Treatment)
	
	names(df)[names(df)=="sortid_stimid_stim"] <- "Peptide"
	df$Peptide <- as.factor(df$Peptide)
	
	names(df)[names(df)=="sortid_stimid_stim_category"] <- "StimCategory"
	df$StimCategory <- as.factor(df$StimCategory)
	
	names(df)[names(df)=="sortid_stimid_stim_type"] <- "StimType"
	df$StimType <- as.factor(df$StimType)

		names(df)[names(df)=="sortid_population"] <- "Population"
	df$Population <- as.factor(df$Population)
	
	names(df)[names(df)=="sortid_stimid"] <- "StimId"
	df$StimId <- as.factor(df$StimId)
	
	df$Label <- as.character(df$Peptide)
	
	if (replicateAsSuffix){
	  df$Label[!is.na(df$Replicate)] <- paste0(df$Label[!is.na(df$Replicate)], '_', df$Replicate[!is.na(df$Replicate)])  
	}
	
	df$Label[df$IsSingleCell] <- paste0(df$Label[df$IsSingleCell], '**')
	
	#TODO: restore this
	#df$EstimatedLibrarySize <- df$metadata_estimatedlibrarysize

	df$Activated <- c(FALSE)
	df$Activated[df$Treatment == 'TAPI-0' & df$Population == 'TNF-Pos'] <- c(TRUE)

	return(df)
}

downloadGenotypes <- function(readsetIds, genomeName, targetField = 'lineages', requireAll = TRUE){
	readsetIds <- unique(readsetIds)  

	df2 <- labkey.selectRows(
  	baseUrl="https://prime-seq.ohsu.edu", 
  	folderPath="/Labs/Bimber", 
  	schemaName="sequenceanalysis", 
  	queryName="sequence_analyses", 
  	colFilter=makeFilter(c("readset", "IN", paste0(readsetIds, collapse=';')), c("library_id/name", "EQUALS", genomeName), c('totalSbtReads', 'GT', 0)),
  	colSelect=c('rowid', 'readset'), 
  	containerFilter=NULL,
  	colNameOpt='rname'
	)
	foundReadsetIds <- unique(df2$readset)  
	analysisIds <- unique(df2$rowid)  

	if (requireAll & length(foundReadsetIds) != length(readsetIds)){
  	stop('Did not find a matching alignment with SBT reads for all readsets')
	}

	df3 <- labkey.selectRows(
  	baseUrl="https://prime-seq.ohsu.edu", 
  	folderPath="/Labs/Bimber", 
  	schemaName="sequenceanalysis", 
  	queryName="alignment_summary_by_lineage", 
  	colFilter=makeFilter(c("analysis_id", "IN", paste0(analysisIds, collapse=';')), c('totalLineages', 'EQUALS', 1)),
  	colSelect=c('analysis_id', 'analysis_id/readset', 'lineages', 'total'), 
  	containerFilter=NULL,
  	colNameOpt='rname'
	)

	merged <- merge(df1, df3, by.x='readset', by.y='analysis_id_readset')
	merged <- merged[,c('rowid', 'lineages', 'total')]
	ret <- cast(merged, lineages ~ rowid, value='total')
	rownames(ret) <- ret$lineages
	ret <- ret[,!(names(ret) %in% c('lineages')), drop = FALSE]

	return(ret)
}

downloadGeneAnnotations <- function(geneIds){
  df  <- data.frame(geneId=geneIds) 
  
  genes.table = NULL
  if (!file.exists("cache.genes.table")) {
    ensembl = useEnsembl(biomart="ensembl", dataset="mmulatta_gene_ensembl")
    attrs <- unique(getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name", "description", "name_1006"), values= df$geneId, mart= ensembl))
    
    attrs <- attrs %>%
      group_by(ensembl_gene_id, external_gene_name, description) %>%
      summarise(goAnnotations = toString(sort(unique(name_1006))))
    attrs <- data.frame(attrs)
    
    save(attrs, file= "cache.genes.table")
  } else {
    load("cache.genes.table")
  }
  
  ret <- merge(x=df, y=attrs, by.x="geneId", by.y="ensembl_gene_id",all.x=T, all.y=F)
  rownames(ret) <- ret$geneId
  
  return(ret)
}
