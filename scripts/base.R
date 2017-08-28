parseGeneCountsMatrix <- function(fileName, allowableColumnIDs){
  geneCounts <- read.table(fileName, sep = '\t', header=TRUE, row.names = 1)
  geneCounts <- geneCounts[, !(names(geneCounts) %in% c('TranscriptId', 'GeneDescription', 'SamplesWithReads', 'GeneName')) ] 
  geneCounts <- geneCounts[!(rownames(geneCounts) %in% c('N_ambiguous', 'N_multimapping', 'N_noFeature', 'N_unmapped')),] 
  
  #debugging
  #geneCounts <- geneCounts[c(1:500),c(2:8)]
  
  colnames(geneCounts) <- gsub('X', '', colnames(geneCounts))
  toSelect <- intersect(allowableColumnIDs, colnames(geneCounts))
  geneCounts <- geneCounts[ toSelect ] 
  
  geneCountMatrix <- data.matrix(geneCounts, rownames.force = NA)
  geneCountMatrix <- geneCountMatrix[ ,colSums(geneCountMatrix) > 0 ] 
  
  return(geneCountMatrix)
}
	
addEnsembl <- function(df){
  ensembl = useEnsembl(biomart="ensembl", dataset="mmulatta_gene_ensembl")
  ensemblIds = unique(df$Ensembl.Id)
  attrs <- unique(getBM(attributes=c('ensembl_gene_id','external_gene_name','description','name_1006'), filters ='ensembl_gene_id', values=ensemblIds, mart = ensembl))
  
  attrs <- attrs %>%
    group_by(ensembl_gene_id, external_gene_name, description) %>%
    summarise(goAnnotations = toString(sort(unique(name_1006))))
  attrs <- data.frame(attrs)
  
  return(merge(df, attrs, by.x='Ensembl.Id', by.y='ensembl_gene_id'))
}

prepareMetadataTable <- function(allowableIds, groupColName, minLibrarySize = 2000000){
	#pull full metadata from LK
	metadataFull <- pullTCRMetaFromLabKey()

	metaUnfilter <- metadataFull[metadataFull$OutputFileId %in% allowableIds,]
	
	rownames(metaUnfilter) <- as.character(metaUnfilter$OutputFileId)
	metaUnfilter$AnimalId <- as.factor(metaUnfilter$AnimalId)
	metaUnfilter$ReadsetId <- as.factor(metaUnfilter$ReadsetId)
	metaUnfilter$HasCDR3s <- c(FALSE)
	metaUnfilter$HasCDR3s[metaUnfilter$NumCDR3s > 0] <- c(TRUE)

	metaUnfilter$EstimatedLibrarySize <- as.numeric(gsub(',','',metaUnfilter$EstimatedLibrarySize))
	metaUnfilter <- metaUnfilter[metaUnfilter$EstimatedLibrarySize > minLibrarySize,]
	metaUnfilter <- metaUnfilter[metaUnfilter$AnimalId != '1247',]
	metaUnfilter <- metaUnfilter[is.na(metaUnfilter$Status),]

	metaUnfilter$GroupCol <- as.factor(metaUnfilter[[groupColName]])

	return(metaUnfilter)
}

prepareMetadataTable2 <- function(metaUnfilter, geneCountMatrix){
	#limit to only those in the gene table, in case we dropped 0 count libraries
	metaUnfilter <- metaUnfilter[metaUnfilter$OutputFileId %in% colnames(geneCountMatrix),]
	metaUnfilter$GroupCol <- as.factor(as.character(metaUnfilter[[groupColName]]))

}

pullTCRMetaFromLabKey <- function(){
	df <- labkey.selectRows(
		baseUrl="https://prime-seq.ohsu.edu", 
		folderPath="/Internal/Bimber/145", 
		schemaName="lists", 
		queryName="TCR_Datasets", 
		viewName="", 
		colSelect=c('ReadsetId/rowid','ReadsetId','ReadsetId/application','ReadsetId/status','ReadsetId/workbook','StimId','StimId/AnimalId','StimId/Date','StimId/Peptide','StimId/Treatment','Population','Replicate','Cells','StimId/ActivatedFreq','StimId/CellClass','StimId/Background','ReadsetId/numCDR3s','ReadsetId/distinctLoci','ReadsetId/numTcrRuns','SequenceComments','GroupId','GeneTable','Activated','Metadata/geneCountFiles','Metadata/estimatedLibrarySize'), 
		containerFilter=NULL,
		colNameOpt='rname'
	)

	print(str(df))
	
	df$AnimalId <- df$stimid_animalid
	df$ReadsetId <- df$readsetid
	df$EstimatedLibrarySize <- df$metadata_estimatedlibrarysize
	df$Peptide <- df$stimid_peptide
	df$CellClass <- df$stimid_cellclass
	df$NumCDR3s <- df$readsetid_numcdr3s
	df$Treatment <- df$stimid_treatment
	df$DistinctLoci <- df$readsetid_distinctloci
	df$OutputFileId <- as.integer(df$metadata_genecountfiles)
	
	df$Activated <- c(FALSE)
	df$Activated[df$Treatment == 'TAPI-0' & df$Population == 'TNF-Pos'] <- c(TRUE)

	return(df)
}
