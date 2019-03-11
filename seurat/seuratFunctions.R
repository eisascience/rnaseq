library(Seurat)

createSeuratObj <- function(seuratData = NA, project = NA, minFeatures = 25){
  seuratObj <- CreateSeuratObject(counts = seuratData, min.cells = 0, min.features = minFeatures, project = project)

  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObj), value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts'))
  seuratObj[['percent.mito']] <- percent.mito
  
  return(seuratObj)
}

appendCellHashing <- function(seuratObj, barcodeFile) {
  barcodeData <- read.table(barcodeFile, sep = ',', header = T, row.names = 1)
  barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% c('no_match', 'total_reads'))),]
  print(paste0('Initial barcodes in HTO data: ', ncol(barcodeData)))

  toDrop <- sum(colSums(barcodeData) == 0)
  if (toDrop > 0){
    print(paste0('cells dropped due to zero HTO counts: ', toDrop))
    barcodeData <- barcodeData[,which(colSums(barcodeData) > 0)]
    print(paste0('Final barcodes in HTO data: ', ncol(barcodeData)))
  }
  
  print(paste0('Initial barcodes in GEX data: ', ncol(seuratObj)))
  
  joint_bcs <- intersect(colnames(barcodeData),colnames(seuratObj))
  print(paste0('Total barcodes shared between HTO and GEX data: ', length(joint_bcs)))
  
  seuratObj <- subset(x = seuratObj, cells = joint_bcs)
  barcodeData <- as.matrix(barcodeData[,joint_bcs])
  
  print(paste0('Initials HTOs: ', length(rownames(barcodeData))))
  toDrop <- sum(rowSums(barcodeData) == 0)        
  if (toDrop > 0){
    print(paste0('HTOs dropped due to zero counts: ', toDrop))
    print(names(which(rowSums(barcodeData) > 0)))
    barcodeData <- barcodeData[which(rowSums(barcodeData) > 0),]
    print(paste0('Final HTOs: ', nrow(barcodeData)))
  }
  
  seuratObj[['HTO']] <- CreateAssayObject(counts = barcodeData)
  
  #doHtoDemux(seuratObj)
  
  return(seuratObj)
}

doHtoDemux <- function(seuratObj) {
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR", display.progress = FALSE)
  seuratObj <- HTODemux(seuratObj, assay = "HTO", positive.quantile =  0.99)
  
  #report outcome
  print(table(seuratObj$HTO_classification.global))
  print(table(seuratObj$hash.ID))
  
  # Group cells based on the max HTO signal
  seuratObj_hashtag <- seuratObj
  Idents(seuratObj_hashtag) <- "hash.ID"
  htos <- rownames(GetAssayData(seuratObj_hashtag,assay = "HTO"))
  for (hto in htos){
    print(RidgePlot(seuratObj_hashtag, features = c(hto), assay = 'HTO', ncol = 1))
  }
  
  print(HTOHeatmap(seuratObj, assay = "HTO", classification = "HTO_classification", global.classification = "HTO_classification.global", ncells = min(3000, ncol(seuratObj)), singlet.names = NULL))
  
  return(seuratObj)
}

printQcPlots1 <- function(seuratObj) {
  print(VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  
  print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  
  #10x-like plot
  nUMI <- Matrix::colSums(GetAssayData(object = seuratObj, slot = "counts"))
  nUMI <- sort(nUMI)
  
  countAbove <-unlist(lapply(nUMI, function(x){
    sum(nUMI >= x)
  }))
  
  #TODO: labeling?
  plot(log(countAbove), log(nUMI), pch=20)  
}

hasStepRun <- function(seuratObj, name) {
  return(!is.null(seuratObj@misc[[paste0(name, 'Run')]]))
}

markStepRun <- function(seuratObj, name, saveFile = NULL) {
  seuratObj@misc[paste0(name, 'Run')] <- T
  if (!is.null(saveFile)){
    saveRDS(seuratObj, file = saveFile)  
  }
  
  return(seuratObj)
}

mergeSeuratObjs <- function(seuratObjs, data){
  for (exptNum in names(data)) {
    print(paste0('adding expt: ', exptNum))
    seuratObjs[[exptNum]] <- RenameCells(object = seuratObjs[[exptNum]], add.cell.id = exptNum)
    seuratObjs[[exptNum]]
    
    if (is.null(seuratObj)) {
      seuratObj <- seuratObjs[[exptNum]]
    } else {
      seuratObj <- merge(x = seuratObj,
                         y = seuratObjs[[exptNum]],
                         project = outPrefix,
                         do.normalize = F
      )
      
    }
    
    print('after merge')
    print(seuratObj)
  }

  return(seuratObj)  
}

processSeurat1 <- function(seuratObj){
  if (!hasStepRun(seuratObj, 'NormalizeData')) {
    seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", do.par = doPar, num.cores = numCores)
    seuratObj <- markStepRun(seuratObj, 'NormalizeData', saveFile)
  }
  
  
  if (!hasStepRun(seuratObj, 'FindVariableFeatures')) {
    seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
    seuratObj <- markStepRun(seuratObj, 'FindVariableFeatures', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'ScaleData')) {
    seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), vars.to.regress = c("nCount_RNA", "percent.mito"), do.par = doPar, num.cores = numCores)
    seuratObj <- markStepRun(seuratObj, 'ScaleData', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'RunPCA')) {
    seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(object = seuratObj), verbose = FALSE)
    seuratObj <- markStepRun(seuratObj, 'RunPCA', saveFile)
  }
  
  
  if (!hasStepRun(seuratObj, 'ProjectDim')) {
    seuratObj <- ProjectDim(object = seuratObj)
    seuratObj <- markStepRun(seuratObj, 'ProjectDim', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'JackStraw')) {
    seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, do.par = T, num.cores = numCores)
    seuratObj <- markStepRun(seuratObj, 'JackStraw', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'ScoreJackStraw')) {
    seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
    seuratObj <- markStepRun(seuratObj, 'ScoreJackStraw')
  }
  
  return(seuratObj)
}

appendTcrClonotypes <- function(seuratData = NA, clonotypeFile = NA){

}