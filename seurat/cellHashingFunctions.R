library(data.table)
library(dplyr)
library(dtplyr)

library(Seurat)
library(reshape2)
library(ggplot2)
library(knitr)
library(KernSmooth)

source('https://raw.githubusercontent.com/chris-mcginnis-ucsf/MULTI-seq/master/R/MULTIseq.Classification.Suite.R')

#TODO: can we make these thresholds adaptive based on our data?
processCiteSeqCount <- function(barcodeFile, minRowSum = 5, minColSum = 5, minRowMax = 40, minColMax = 5, minMeanNonZeroCount = 5) {
  barcodeData <- read.table(barcodeFile, sep = ',', header = T, row.names = 1)
  barcodeData <- barcodeData[which(!(rownames(barcodeData) %in% c('no_match', 'total_reads'))),]
  print(paste0('Initial barcodes in HTO data: ', ncol(barcodeData)))
  
  barcodeData <- doCellFiltering(barcodeData, minColSum, minColMax)

  #rowsum
  toDrop <- sum(rowSums(barcodeData) < minRowSum)
  if (toDrop > 0){
    print(paste0('HTOs dropped due to low total counts: ', toDrop))
    print(paste(rownames(barcodeData)[which(rowSums(barcodeData) < minRowSum)], collapse = ', '))
    barcodeData <- barcodeData[which(rowSums(barcodeData) >= minRowSum),]
    print(paste0('HTOs after filter: ', nrow(barcodeData)))
  }
  
  #summarize
  rowSummary <- generateByRowSummary(barcodeData)
  print(kable(rowSummary, caption = 'HTO Summary', row.names = F))
  
  #Drop HTOs with zero strong cells:
  toDrop <- rowSummary$max < minRowMax
  if (sum(toDrop) > 0){
    print(paste0('HTOs dropped due to low max counts: ', sum(toDrop)))
    print(paste(rownames(barcodeData)[toDrop], collapse = ', '))
    barcodeData <- barcodeData[!toDrop,]
    print(paste0('HTOs after filter: ', nrow(barcodeData)))
  }
  
  #Now filter HTO by mean count among non-zero cells:
  barcodeMatrix <- as.matrix(barcodeData)
  meanNonZero <- (rowSums(barcodeMatrix) / rowSums(!!barcodeMatrix))
  print(kable(data.frame(HTO = rownames(barcodeMatrix), MeanCountOfNonZeroCells = meanNonZero), row.names = F))

  toDrop <- meanNonZero < minMeanNonZeroCount
  if (sum(toDrop) > 0){
    print(paste0('HTOs dropped due to insufficient mean non-zero count: ', sum(toDrop)))
    print(paste(rownames(barcodeData)[toDrop], collapse = ', '))
    barcodeData <- barcodeData[!toDrop,]
    print(paste0('HTOs after filter: ', nrow(barcodeData)))
    print(paste(rownames(barcodeData), collapse = ', '))
  }
  
  # repeat colsum filter.  because we potentially dropped rows, some cells might be filtered out:
  barcodeData <- doCellFiltering(barcodeData, minColSum, minColMax)
  
  #repeat summary:
  if (nrow(barcodeData) == 0) {
    print('No HTOs remaining')
  } else {
    rowSummary <- generateByRowSummary(barcodeData)
    print(kable(rowSummary, caption = 'HTO Summary After Filter', row.names = F))
  }
  
  return(barcodeData)  
}

generateByRowSummary <- function(barcodeData) {
  barcodeMatrix <- as.matrix(barcodeData)
  df <- data.frame(HTO = rownames(barcodeData), min = apply(barcodeData, 1, min), max = apply(barcodeData, 1, max), mean = apply(barcodeData, 1, mean), nonzero = apply(barcodeData, 1, function(x){
    sum(x > 0)
  }), mean_nonzero = (rowSums(barcodeMatrix) / rowSums(!!barcodeMatrix)), total_gt1 = apply(barcodeMatrix, 1, function(x){
    sum(x > 1)  
  }), mean_gt1 = apply(barcodeMatrix, 1, function(x){
    mean(sapply(x, function(y){if (y > 1) y else NA}), na.rm = T)  
  }))
  
  return(df)
}

doCellFiltering <- function(barcodeData, minColSum, minColMax){
  #colsum
  toDrop <- sum(colSums(barcodeData) < minColSum)
  if (toDrop > 0){
    print(paste0('cells dropped due to low total counts per column: ', toDrop))
    barcodeData <- barcodeData[,which(colSums(barcodeData) >= minColSum)]
    print(paste0('Final cell barcodes: ', ncol(barcodeData)))
  }
  
  #colmax
  barcodeMatrix <- as.matrix(barcodeData)
  cm <- apply(barcodeMatrix, 2, max)
  toDrop <- sum(cm < minColMax)
  if (toDrop > 0){
    print(paste0('cells dropped due to low max counts per column: ', toDrop))
    barcodeData <- barcodeData[,(cm >= minColMax)]
    print(paste0('Final cell barcodes: ', ncol(barcodeData)))
  }
  
  return(barcodeData)
}

generateQcPlots <- function(barcodeData){
  print('Generating QC Plots')
  
  #Plot counts/cell:
  countsPerCell <- Matrix::colSums(barcodeData)
  countsPerCell <- sort(countsPerCell)
  countAbove <-unlist(lapply(countsPerCell, function(x){
    sum(countsPerCell >= x)
  }))
  plot(log(countAbove), log(countsPerCell), pch=20, ylab = "log(Reads/Cell)", xlab = "log(Total Cells)")  
  
  topBarcodes <- sort(tail(countsPerCell, n = 20), decreasing = T)
  
  print(kable(data.frame(CellBarcode = names(topBarcodes), Count = topBarcodes), row.names = F))

  #boxplot per HTO:
  barcodeMatrix <- as.matrix(barcodeData)
  melted <- setNames(melt(barcodeMatrix), c('HTO', 'CellBarcode', 'Count'))
  print(ggplot(melted, aes(x = HTO, y = Count)) +
      geom_boxplot() +
      xlab('HTO') +
      ylab('Count') +
      ggtitle('Counts By HTO') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
  
  melted$Count <- melted$Count + 0.5
  print(ggplot(melted, aes(x = HTO, y = Count)) +
          geom_boxplot() +
          xlab('HTO') +
          scale_y_continuous(trans='log10') +
          ylab('Count') +
          ggtitle('Counts By HTO (log)') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
  
  
  #normalize columns, print top barcode fraction:
  normalizedBarcodes <- sweep(barcodeMatrix,2,colSums(barcodeMatrix),`/`)
  topValue <- apply(normalizedBarcodes,2,function(x){
    max(x)
  })

  df <- data.frame(Barcode1 = topValue)
  print(ggplot(df, aes(x = Barcode1)) +
          geom_histogram(binwidth = 0.05) +
          xlab('Top Barcode Fraction') +
          ylab('Count')
  )
  
  print(paste0('Total cells where top barcode is >0.75 of counts: ', length(topValue > 0.75)))
  
}

generateCellHashCallsSeurat <- function(barcodeData) {
  seuratObj <- CreateSeuratObject(barcodeData, assay = 'HTO')   
  
  tryCatch({
    seuratObj <- doHtoDemux(seuratObj)
    
    return(data.table(Barcode = as.factor(colnames(seuratObj)), HTO_classification = seuratObj$hash.ID, HTO_classification.all = seuratObj$HTO_classification, HTO_classification.global = seuratObj$HTO_classification.global, key = c('Barcode')))
  }, error = function(e){
    print(e)
    print('Error generating seurat calls, aborting')
    return(NA)
  })
}

appendCellHashing <- function(seuratObj, barcodeCallTable) {
  print(paste0('Initial called barcodes in HTO data: ', nrow(barcodeCallTable)))
  print(paste0('Initial barcodes in GEX data: ', ncol(seuratObj)))
  
  joint_bcs <- intersect(barcodeCallTable$CellBarcode,colnames(seuratObj))
  print(paste0('Total barcodes shared between HTO and GEX data: ', length(joint_bcs)))
  
  seuratObj <- subset(x = seuratObj, cells = joint_bcs)
  barcodeCallTable <- barcodeCallTable[colnames(seuratObj),]
  
  seuratObj[['HTO']] <- barcodeCallTable$HTO
  seuratObj[['HTO_Classification']] <- barcodeCallTable$HTO_Classification
  
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

generateCellHashCallsMultiSeq <- function(barcodeData, showTSNE = T) {
  #transpose CITE-seq count input:
  bar.table.full <- data.frame(t(barcodeData))
  
  if (showTSNE) {
    bar.tsne <- barTSNE(bar.table.full)
  
    for (i in 3:ncol(bar.tsne)) {
      g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
        geom_point() +
        scale_color_gradient(low = "black", high = "red") +
        ggtitle(colnames(bar.tsne)[i]) +
        theme(legend.position = "right")
      print(g)
    }
  }
  
  final.calls <- NA
  neg.cells <- c()
  r1 <- performRoundOfMultiSeqCalling(bar.table.full, 1)
  if (is.list(r1)){
    neg.cells <- c(neg.cells, r1$neg.cells)
    final.calls <- r1$final.calls
    
    if (length(r1$neg.cells > 0)) {
      r2 <- performRoundOfMultiSeqCalling(r1$bar.table, 2)
      if (is.list(r2)){
        neg.cells <- c(neg.cells, r2$neg.cells)
        final.calls <- r2$final.calls
    
        if (length(r2$neg.cells > 0)) {    
          r3 <- performRoundOfMultiSeqCalling(r2$bar.table, 3)
          if (is.list(r3)){
            neg.cells <- c(neg.cells, r3$neg.cells)
            final.calls <- r3$final.calls
          }
        }
      }
    }
  }
  
  if (all(is.na(final.calls))) {
    print('No calls, aborting')
    return(NA)
  }
  
  neg.cells <- unique(neg.cells)
  final.names <- c(names(final.calls),neg.cells)
  
  #MultiSeq replaces hyphens in names
  final.calls <- gsub(x = final.calls, pattern = '\\.', '-')
  final.calls <- c(final.calls, rep("Negative",length(neg.cells)))
  names(final.calls) <- final.names
  
  print(table(final.calls))

  global <- as.character(final.calls)
  global[!(global %in% c('Doublet', 'Negative'))] <- 'Singlet'
  global <- as.factor(global)
  
  if (length(final.calls) == 0){
    return(NA)
  }
  
  return(data.table(
    Barcode = as.factor(names(final.calls)), 
    HTO_classification = as.factor(final.calls), 
    HTO_classification.all = as.factor(final.calls), 
    HTO_classification.global = global,
    key = c('Barcode')))
}

performRoundOfMultiSeqCalling <- function(bar.table, roundNum) {
  ## Perform Quantile Sweep
  print(paste0("Round ", roundNum ," calling..."))
  print(paste0('Initial cells: ', nrow(bar.table)))
  bar.table_sweep.list <- list()
  n <- 0
  for (q in seq(0.01, 0.99, by=0.02)) {
    #print(q)
    n <- n + 1
    bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
  }
  
  ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
  threshold1 <- findThresh(call.list=bar.table_sweep.list)
  
  print(ggplot(data=threshold1$res, aes(x=q, y=Proportion, color=Subset)) + 
    geom_line() + 
    theme(legend.position = "right") +
    geom_vline(xintercept=threshold1$extrema, lty=2) + 
    scale_color_manual(values=c("red","black","blue")) +
    ggtitle(paste0("Round ", roundNum)                         
    )
  )
  
  if (length(threshold1$extrema) == 0){
    print("Unable to find extrema, attempting to default to 0.9")
    threshold1$extrema <- c(0.15)
  }
  
  ## Finalize round 1 classifications, remove negative cells
  extrema <- threshold1$extrema[length(threshold1$extrema)]  #assume we use max value
  print(paste0('Round ', roundNum ,' Threshold: ', extrema))
  round1.calls <- classifyCells(bar.table, q=findQ(threshold1$res, extrema))
  neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
  print(paste0('Negative cells dropped: ', length(neg.cells)))
  if (length(neg.cells) > 0) {
    bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
    print(paste0('Remaining: ', nrow(bar.table)))
  }
  
  return(list(
    'bar.table' = bar.table,
    'neg.cells' = neg.cells,
    'final.calls' = round1.calls
  ))  
}

reclassifyByMultiSeq <- function(bar.table.full, final.calls){
  reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
  reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)
  
  ## Visualize Results
  print(ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
    geom_point() + xlim(c(nrow(pool.reclass.res)-1,1)) + 
    ylim(c(0,1.05)) +
    geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)
  )
  
  ## Finalize negative cell rescue results
  final.calls.rescued <- final.calls
  rescue.ind <- which(reclass.cells$ClassStability >= 16) ## Note: Value will be dataset-specific
  final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
}

processEnsemblHtoCalls <- function(mc, sc, barcodeData, outFile = 'combinedHtoCalls.txt', allCallsOutFile = NA) {
  if (all(is.na(sc)) && all(is.na(mc))){
    print('MULTI-seq and Seurat failed to produce calls, aborting')
    return()
  }
  
  if (all(is.na(sc))){
    print('No calls for seurat found')  
    dt <- data.table(CellBarcode = mc$Barcode, HTO = mc$HTO_classification, HTO_Classification = mc$HTO_classification.global, key = 'CellBarcode', Seurat = c(F), MultiSeq = c(T))
    dt <- printFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)
    
    return(dt)
  }
  
  if (all(is.na(mc))){
    print('No calls for MULTI-seq found')  
    dt <- data.table(CellBarcode = sc$Barcode, HTO = sc$HTO_classification, HTO_Classification = sc$HTO_classification.global, key = 'CellBarcode', Seurat = c(T), MultiSeq = c(F))
    dt <- printFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)
    
    return(dt)
  }
  
  mc$Barcode <- as.character(mc$Barcode)
  sc$Barcode <- as.character(sc$Barcode)
  merged <- merge(mc, sc, all = T, by = 'Barcode', suffixes = c('.MultiSeq', '.Seurat'))
  
  merged$Concordant <- as.character(merged$HTO_classification.MultiSeq) == as.character(merged$HTO_classification.Seurat)
  merged$ConcordantNoNeg <- !(!merged$Concordant & merged$HTO_classification.MultiSeq != 'Negative' & merged$HTO_classification.Seurat != 'Negative')
  merged$GlobalConcordant <- as.character(merged$HTO_classification.global.MultiSeq) == as.character(merged$HTO_classification.global.Seurat)
  merged$HasSeuratCall <- !is.na(merged$HTO_classification.Seurat) & merged$HTO_classification.Seurat != 'Negative'
  merged$HasMultiSeqCall <- !is.na(merged$HTO_classification.MultiSeq) & merged$HTO_classification.MultiSeq != 'Negative'
  
  print(paste0('Total concordant: ', nrow(merged[merged$Concordant])))
  print(paste0('Total discordant: ', nrow(merged[!merged$Concordant])))
  print(paste0('Total discordant, ignoring negatives: ', nrow(merged[!merged$ConcordantNoNeg])))
  print(paste0('Total discordant global calls: ', nrow(merged[!merged$GlobalConcordant])))
  
  discord <- merged[!merged$GlobalConcordant]
  discord <- discord %>% group_by(HTO_classification.global.MultiSeq, HTO_classification.global.Seurat) %>% summarise(Count = n())
  
  print(qplot(x=HTO_classification.global.MultiSeq, y=HTO_classification.global.Seurat, data=discord, fill=Count, geom="tile") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle('Discordance By Global Call') + ylab('Seurat') + xlab('MULTI-seq')
  )
  
  discord <- merged[!merged$Concordant]
  discord <- discord %>% group_by(HTO_classification.MultiSeq, HTO_classification.Seurat) %>% summarise(Count = n())
  
  print(qplot(x=HTO_classification.MultiSeq, y=HTO_classification.Seurat, data=discord, fill=Count, geom="tile") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle('Discordance By HTO Call') + ylab('Seurat') + xlab('MULTI-seq')
  )
  
  discord <- merged[!merged$ConcordantNoNeg]
  discord <- discord %>% group_by(HTO_classification.MultiSeq, HTO_classification.Seurat) %>% summarise(Count = n())
  
  print(qplot(x=HTO_classification.MultiSeq, y=HTO_classification.Seurat, data=discord, fill=Count, geom="tile") + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggtitle('Discordance By HTO Call, Ignoring Negatives') + ylab('Seurat') + xlab('MULTI-seq')
  )
  ret <-merged[merged$ConcordantNoNeg,]
  
  # These calls should be identical, except for possibly negatives from one method that are non-negative in the other
  # For the time being, accept those as correct.
  ret$FinalCall <- ret$HTO_classification.MultiSeq
  ret$FinalCall[ret$HTO_classification.MultiSeq == 'Negative'] <- ret$HTO_classification.Seurat[ret$HTO_classification.MultiSeq == 'Negative']
  
  ret$FinalClassification <- ret$HTO_classification.global.MultiSeq
  ret$FinalClassification[ret$HTO_classification.global.MultiSeq == 'Negative'] <- ret$HTO_classification.global.Seurat[ret$HTO_classification.global.MultiSeq == 'Negative']

  if (!is.na(allCallsOutFile) && nrow(merged) > 0) {
    write.table(merged, file = allCallsOutFile, row.names = F, sep = '\t', quote = F)
  }
  
  if (nrow(ret) > 0){
    dt <- data.table(CellBarcode = ret$Barcode, HTO = ret$FinalCall, HTO_Classification = ret$FinalClassification, key = 'CellBarcode', Seurat = ret$HasSeuratCall, MultiSeq = ret$HasMultiSeqCall)
    dt <- printFinalSummary(dt, barcodeData)
    write.table(dt, file = outFile, row.names = F, sep = '\t', quote = F)
    
    return(dt)
    
  } else {
    print('No rows, not saving ')
  }
}

printFinalSummary <- function(dt, barcodeData){
  #Append raw counts:
  bc <- t(barcodeData)
  x <- melt(bc)
  names(x) <- c('CellBarcode', 'HTO', 'Count')
  
  merged <- merge(dt, x, by = c('CellBarcode', 'HTO'), all.x = T, all.y = F)
  
  bc <- as.data.frame(bc)
  bc$CellBarcode <- rownames(bc)
  merged <- merge(merged, bc, by = c('CellBarcode'), all.x = T, all.y = F)
 
  t <- table(SeuratCall = merged$Seurat, MultiSeqCall = merged$MultiSeq)
  
  colnames(t)[colnames(t) == T] <- c('MultiSeq Call')
  colnames(t)[colnames(t) == F] <- c('MultiSeq No Call')
  
  rownames(t)[rownames(t) == T] <- c('Seurat Call')
  rownames(t)[rownames(t) == F] <- c('Seurat No Call')
  
  print(kable(t))
  
  print(ggplot(merged, aes(x = HTO)) +
          geom_histogram(stat = 'count') +
          xlab('HTO') +
          ylab('Count') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
  
  print(kable(table(Classification = merged$HTO)))
  
  print(ggplot(merged, aes(x = HTO_Classification)) +
          geom_histogram(stat = 'count') +
          xlab('Classification') +
          ylab('Count') +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
  
  print(kable(table(Classification = merged$HTO_Classification)))
  
    
  
  return(merged)
}