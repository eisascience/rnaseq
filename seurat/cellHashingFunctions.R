library(Seurat)
library(data.table)
library(reshape2)
library(ggplot2)
library(knitr)
library(KernSmooth)

processCiteSeqCount <- function(barcodeFile) {
  barcodeData <- read.table(barcodeFile, sep = ',', header = T, row.names = 1)
  barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% c('no_match', 'total_reads'))),]
  print(paste0('Initial barcodes in HTO data: ', ncol(barcodeData)))
  
  toDrop <- sum(colSums(barcodeData) == 0)
  if (toDrop > 0){
    print(paste0('cells dropped due to zero HTO counts: ', toDrop))
    barcodeData <- barcodeData[,which(colSums(barcodeData) > 0)]
    print(paste0('Final barcodes in HTO data: ', ncol(barcodeData)))
  }

  toDrop <- sum(rowSums(barcodeData) == 0)
  if (toDrop > 0){
    print(paste0('HTOs dropped due to zero cells with counts: ', toDrop))
    print(paste(rownames(barcodeData)[which(rowSums(barcodeData) == 0)], collapse = ', '))
    barcodeData <- barcodeData[which(rowSums(barcodeData) > 0),]
    print(paste0('Final HTOs: ', nrow(barcodeData)))
  }
  
  return(barcodeData)  
}

generateQcPlots <- function(barcodeData){
  #Plot counts/cell:
  countsPerCell <- Matrix::colSums(barcodeData)
  countsPerCell <- sort(countsPerCell)
  
  countAbove <-unlist(lapply(countsPerCell, function(x){
    sum(countsPerCell >= x)
  }))
  
  plot(log(countAbove), log(countsPerCell), pch=20, ylab = "log(Reads/Cell)", xlab = "log(Total Cells)")  
  
  print(kable(sort(tail(countsPerCell, n = 20), decreasing = T)))
  barcodeMatrix <- as.matrix(barcodeData)
  
  #TODO: smart plot removing outlier???
  #boxplot per HTO:
  melted <- setNames(melt(barcodeMatrix), c('HTO', 'CellBarcode', 'Count'))
  print(ggplot(melted, aes(x = HTO, y = Count)) +
      geom_boxplot() +
      xlab('HTO') +
      ylab('Count') +
      ggtitle('Counts By HTO')
  )
  
  print(ggplot(melted, aes(x = HTO, y = Count)) +
          geom_boxplot() +
          xlab('HTO') +
          #TODO: maybe we should add 0.5 to counts to avoid zeros?
          scale_y_continuous(trans='log10') +
          ylab('Count') +
          ggtitle('Counts By HTO (log)')
  )
  
  
  #normalize columns, print top 2 per cell:
  normalizedBarcodes <- sweep(barcodeMatrix,2,colSums(barcodeMatrix),`/`)
  topValue <- apply(normalizedBarcodes,2,function(x){
    max(x)
  })
  secondValue <- apply(normalizedBarcodes,2,function(x){
    n <- length(x)
    sort(x)[n-1]
  })
  
  #TODO: does this make sense?
  df <- data.frame(Barcode1 = topValue, Barcode2 = secondValue)
  print(ggplot(df, aes(x = Barcode1, y = Barcode2)) +
          geom_point() +
          xlab('Top Barcode Fraction') +
          ylab('Second Barcoe Fraction')
  )
  
}

generateCellHashCallsSeurat <- function(barcodeData) {
  seuratObj <- CreateSeuratObject(barcodeData, assay = 'HTO')   
  
  tryCatch({
    seuratObj <- doHtoDemux(seuratObj)
  }, warning = function(w){
    print(w)
  }, error = function(e){
    print(e)
  })
  
  #TODO: we probably should return a simple dataframe, not a seurat object
  return(seuratObj)
}

appendCellHashing <- function(seuratObj, barcodeData) {
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

generateCellHashCallsMultiSeq <- function(barcodeData) {
  
}