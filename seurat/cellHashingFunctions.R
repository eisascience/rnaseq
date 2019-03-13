library(Seurat)
library(data.table)

generateCellHashCalls <- function(barcodeFile) {
  barcodeData <- fread(barcodeFile, sep = ',', header = T, row.names = 1)
  barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% c('no_match', 'total_reads'))),]
  print(paste0('Initial barcodes in HTO data: ', ncol(barcodeData)))
  
  toDrop <- sum(colSums(barcodeData) == 0)
  if (toDrop > 0){
    print(paste0('cells dropped due to zero HTO counts: ', toDrop))
    barcodeData <- barcodeData[,which(colSums(barcodeData) > 0)]
    print(paste0('Final barcodes in HTO data: ', ncol(barcodeData)))
  }
  
  seuratObj <- CreateSeuratObject(barcodeData, assay = 'HTO')   
  
  tryCatch({
    doHtoDemux(seuratObj)
  }, warning = function(w){
    print(w)
  }, error = function(e){
    print(e)
  })
}

appendCellHashing <- function(seuratObj, barcodeFile) {
  barcodeData <- fread(barcodeFile, sep = ',', header = T, row.names = 1)
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
  
  tryCatch({
    doHtoDemux(seuratObj)
  }, warning = function(w){
    print(w)
  }, error = function(e){
    print(e)
  })
  
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