library('DESeq2')
library('ggplot2')
library('data.table')
library(gridExtra)
library(grid)

options(echo=TRUE);
setwd('/work')

load('/work/step4.rdata')

groupColName <- 'Peptide'

deSeqCounts <- function(geneTable, dds, ct, geneId, geneName){
  d1 <- plotCounts(dds, geneId, intgroup = 'OutputFileId', normalized = FALSE, main = 'Title', xlab='OutputFileId', transform = FALSE, returnData = TRUE, replaced = FALSE)  
  d2 <- plotCounts(dds, geneId, intgroup = 'OutputFileId', normalized = TRUE, main = 'Title', xlab='OutputFileId', transform = FALSE, returnData = TRUE, replaced = FALSE)  
  d3 <- plotCounts(dds, geneId, intgroup = 'OutputFileId', normalized = TRUE, main = 'Title', xlab='OutputFileId', transform = TRUE, returnData = TRUE, replaced = FALSE)  
  d4 <- plotCounts(dds, geneId, intgroup = 'OutputFileId', normalized = FALSE, main = 'Title', xlab='OutputFileId', transform = TRUE, returnData = TRUE, replaced = FALSE)  
  
  toAdd <- data.frame(OutputFileId=d2$OutputFileId, deseq.norm=d2$count)
  #colnames(toAdd)[colnames(toAdd) == 'deseq'] <- paste0('deseq.', geneName)
  colnames(toAdd)[colnames(toAdd) == 'deseq.norm'] <- paste0('deseq.norm.', geneName)
  #colnames(toAdd)[colnames(toAdd) == 'deseq.transform'] <- paste0('deseq.transform.', geneName)
  #colnames(toAdd)[colnames(toAdd) == 'deseq.transform.nonorm'] <- paste0('deseq.transform.nonorm.', geneName)
  
  geneTable <- merge(geneTable, toAdd, by='OutputFileId', suffixes = c('', geneName))

  return(geneTable)
}

joinGene <- function(withGenes, ct, geneName, colName){
  toJoin <- data.frame(counts=ct[geneName,], OutputFileId=as.integer(colnames(ct)))
  withGenes$OutputFileId <- as.integer(withGenes$OutputFileId)
  withGenes <- merge(withGenes, toJoin, by.x='OutputFileId', by.y='OutputFileId')
  names(withGenes)[names(withGenes) == 'counts'] <- colName
  
  return(withGenes)
}

findGeneId <- function(name){
  allGenes <- rownames(geneCountMatrix)
  ret <- allGenes[grepl(paste0('\\|',name,'$'), allGenes)]
  if (length(ret) == 1){
    return(ret[1])
  } else if (length(ret) > 1) {
    print(ret)
  }
  
  return(NA)
}

makePlot2 <- function(withGenes, geneName){
  png(paste0('/work/', geneName, '.png'), width=1200, height=600)
  toPlot <- data.frame(withGenes)
  toPlot$Counts = toPlot[,colnames(toPlot)[colnames(toPlot) == geneName]]
  toPlot$CountsN = toPlot[,colnames(toPlot)[colnames(toPlot) == paste0('cpmN.',geneName)]]
  #toPlot$DeSeq = toPlot[,colnames(toPlot)[colnames(toPlot) == paste0('deseq.',geneName)]]
  toPlot$DeSeqNorm = toPlot[,colnames(toPlot)[colnames(toPlot) == paste0('deseq.norm.', geneName)]]
  #toPlot$DeSeqTrans = toPlot[,colnames(toPlot)[colnames(toPlot) == paste0('deseq.transform.', geneName)]]
  #toPlot$DeSeqTransNoNorm = toPlot[,colnames(toPlot)[colnames(toPlot) == paste0('deseq.transform.nonorm.', geneName)]]
  #toPlot$mask <- c(0)
  #toPlot$mask[toPlot$DeSeqNorm < 500] <- c(1)
  
  #P1 <- ggplot(toPlot, aes(x=Peptide, y=Counts, col=AnimalId, label=SubjectId)) + geom_text(fontface = "bold") + geom_jitter(width = 0.2) + labs(title = 'Counts Per Million, Not Normalized') + theme(legend.position="none")
  #P2 <- ggplot(toPlot, aes(x=Peptide, y=CountsN, col=AnimalId, label=SubjectId)) + geom_text(fontface = "bold") + geom_jitter(width = 0.2) + labs(title = 'Counts Per Million, Normalized') + theme(legend.position="none")
  #P1 <- P1 + coord_trans(y = "log10") 
  #P2 <- P2 + coord_trans(y = "log10") 
  
  cutoff <- data.frame(yintercept=1, cutoff=factor(0))
  #DE1 <- ggplot(toPlot, aes(x=Peptide, y=DeSeq, col=HasCDR3s, label=SubjectId)) + geom_text(fontface = "bold") + geom_jitter(width = 0.2) + labs(title = 'DESeq2 Counts: No Norm/Transform') + theme(legend.position="none")
  DE1 <- ggplot(toPlot, aes(x=AnimalId, y=DeSeqNorm, col=AnimalId, label=AnimalId)) + 
    #geom_text(fontface = "bold",position=position_jitter(width=0.2,height=0.2)) + 
    geom_point(position=position_jitter(width=0.1,height=0.1)) + 
    labs(title = paste0('DESeq2 Normalized Counts: ', geneName)) + 
    geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff, show_guide=TRUE) +
    ylim(0.5,NA) +
    theme(legend.position="none") 
  
  DE2 <- DE1 + coord_trans(y = "log10")
  DE2 <- DE2 + labs(title = paste0('DESeq2 Normalized Counts, Log Scale: ', geneName)) 
  
    #facet_grid(mask ~ ., scales="free", space="fixed") + 
    #theme(strip.background = element_blank(), strip.text.y = element_blank())
  
  #DE3 <- ggplot(toPlot, aes(x=Peptide, y=DeSeqTrans, col=HasCDR3s, label=SubjectId)) + geom_text(fontface = "bold",position=position_jitter(width=0.2,height=0.2)) + labs(title = paste0('DESeq2 Transformed/Normalized Counts, : ', geneName)) + theme(legend.position="none")
  #DE4 <- ggplot(toPlot, aes(x=Peptide, y=DeSeqTransNoNorm, col=HasCDR3s, label=SubjectId)) + geom_text(fontface = "bold",position=position_jitter(width=0.2,height=0.2)) + labs(title = paste0('DESeq2 Transformed Counts, No Normalization: ', geneName)) + theme(legend.position="none")
  
  grid.arrange(DE1, DE2, ncol=2)
  
  dev.off()
}

ct <- fpm(dds, robust=TRUE)
ctNorm <- counts(dds, normalized=TRUE)

genes <- c('CCL2', 'CCL19', 'CCL21', 'DPP4', 'GLTSCR2', 'TLR4', 'IL12A', 'IL12B', 'CSF2', 'CCR3', 'CXCR3', 'CD52', 'SIGLEC6', 'SIGLEC10', 'PRF1', 'FCER2', 'CXCL13', 'PRKCE', 'TLR5', 'TCL1A', 'BTLA', 'CXCL3', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCR4', 'CCR4', 'EBF1', 'BLK', 'BTG2', 'TCF3', 'NCAM1', 'CD69','TNF', 'FCGR3', 'CCR8', 'CD3G','FAS', 'CD28', 'CD27', 'ABCB1','RORC', 'PTPRCAP', 'CD8A','CD8B','B3GAT1','CCL1','CCL2', 'CCL4','CCL4L2', 'CCL5','CCL8', 'CCL17', 'CCR5', 'CCR7', 'CD160','CD4','CD44','FOXP3','GATA3','HAVCR2','ICOS','IFNG','IL2','IL2RB','IL4','IL4R','IL5','IL6','IL9','IL10','IL13','IL15','IL17A', 'IL17F', 'IL21', 'IL22', 'KLRB1','KLRC1','KLRC2','KLRC3','KLRD1','NKG2D','PD1','PTPRC','TBX21','CD1C','CD19','GZMA','GZMB','GZMK', 'TGFB3')
withGenes <- data.table(metaUnfilter)
withGenes$CellClass <- as.character(withGenes$CellClass)
withGenes$CellClass[withGenes$Population == 'TNF-Neg'] <- 'TNF-Neg'
withGenes$CellClass <- as.factor(withGenes$CellClass)
for (geneName in genes){
  id <- findGeneId(geneName)
  if (is.na(id)){
    print(paste0('not found: ', geneName))
    next
  }
  
  if (id %in% rownames(ct)){
    withGenes <- joinGene(withGenes, ct, id, geneName)  
    withGenes <- joinGene(withGenes, ctNorm, id, paste0('cpmN.', geneName))  
    withGenes <- deSeqCounts(withGenes, dds, ct, id, geneName)
    makePlot2(withGenes, geneName)
  } else {
    print(paste0('no counts: ', id))
  }
}

write.table(withGenes, 'metadataWithGenes.txt', sep='\t', row.names = FALSE, quote = FALSE)