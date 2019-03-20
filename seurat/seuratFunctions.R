library(Seurat)
library(Rlabkey)
library(data.table)
library(dplyr)
library(naturalsort)
library(DropletUtils)

labkey.setDefaults(baseUrl = "https://prime-seq.ohsu.edu")

createSeuratObj <- function(seuratData = NA, project = NA, minFeatures = 25){
  seuratObj <- CreateSeuratObject(counts = seuratData, min.cells = 0, min.features = minFeatures, project = project)
  
  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObj), value = TRUE)
  p.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts'))
  seuratObj[['p.mito']] <- p.mito
  
  return(seuratObj)
}

printQcPlots1 <- function(seuratObj) {
  print(VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA", "p.mito"), ncol = 3))
  
  print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "p.mito"))
  print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  
  #10x-like plot
  nUMI <- Matrix::colSums(GetAssayData(object = seuratObj, slot = "counts"))
  nUMI <- sort(nUMI)
  
  countAbove <-unlist(lapply(nUMI, function(x){
    sum(nUMI >= x)
  }))
  
  plot(log(countAbove), log(nUMI), pch=20, ylab = "UMI/Cell", xlab = "# Cells")  
}

performEmptyDropletFiltering <- function(seuratRawData, fdrThreshold=0.01) {
  br.out <- barcodeRanks(seuratRawData)
  
  # Making a plot.
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=br.out$knee, col="dodgerblue", lty=2)
  abline(h=br.out$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  
  e.out <- emptyDrops(seuratRawData)
  e.out$is.cell <- e.out$FDR <= fdrThreshold
  print(paste0('Passing cells: ', sum(e.out$is.cell, na.rm=TRUE)))
  print(paste0('Failing cells: ', sum(!e.out$is.cell, na.rm=TRUE)))
  print(table(Limited=e.out$Limited, Significant=e.out$is.cell))
  
  toPlot <- e.out[e.out$LogProb != -Inf,]
  if (nrow(toPlot) > 0) {
    print(plot(toPlot$Total, -toPlot$LogProb, col=ifelse(toPlot$is.cell, "red", "black"), xlab="Total UMI count", ylab="-Log Probability"))
  } else {
    print('Probabilities all -Inf, unable to plot')  
  }
  
  passingCells <- rownames(e.out)[e.out$is.cell]
  return(seuratRawData[,passingCells])
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

processSeurat1 <- function(seuratObj, saveFile = NULL, doCellCycle = T){
  if (!hasStepRun(seuratObj, 'NormalizeData')) {
    seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize")
    seuratObj <- markStepRun(seuratObj, 'NormalizeData', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'FindVariableFeatures')) {
    seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
    seuratObj <- markStepRun(seuratObj, 'FindVariableFeatures', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'ScaleData')) {
    seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), vars.to.regress = c("nCount_RNA", "percent.mito"))
    seuratObj <- markStepRun(seuratObj, 'ScaleData')
  }
  
  if (!hasStepRun(seuratObj, 'RunPCA')) {
    seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(object = seuratObj), verbose = FALSE)
    seuratObj <- markStepRun(seuratObj, 'RunPCA')
  }
  
  if (doCellCycle & !hasStepRun(seuratObj, 'CellCycle')) {
    seuratObj <- removeCellCycle(seuratObj)
    seuratObj <- markStepRun(seuratObj, 'CellCycle', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'ProjectDim')) {
    seuratObj <- ProjectDim(object = seuratObj)
    seuratObj <- markStepRun(seuratObj, 'ProjectDim')
  }
  
  if (!hasStepRun(seuratObj, 'FindNeighbors')) {
      seuratObj <- FindNeighbors(object = seuratObj)
      seuratObj <- markStepRun(seuratObj, 'FindNeighbors', saveFile)
   }
  
  if (!hasStepRun(seuratObj, 'JackStraw')) {
    seuratObj <- JackStraw(object = seuratObj, num.replicate = 100)
    seuratObj <- markStepRun(seuratObj, 'JackStraw', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'ScoreJackStraw')) {
    seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
    seuratObj <- markStepRun(seuratObj, 'ScoreJackStraw')
  }
  
  return(seuratObj)
}

appendTcrClonotypes <- function(seuratObject = NA, clonotypeFile = NA){
  tcr <- processAndAggregateTcrClonotypes(clonotypeFile)
  
  origRows <- nrow(tcr)
  tcr <- tcr[tcr$barcode %in% rownames(seuratObject),]
  
  pct <- nrow(tcr) / origRows * 100
  
  print(paste0('Barcodes with clonotypes: ', origRows, ', intersecting with GEX data: ', nrow(tcr), " (", pct, "%)"))
  
  merged <- merge(data.frame(barcode = colnames(seuratObj)), tcr, by = c('barcode'), all.x = T)
  rownames(merged) <- merged$barcode
  
  for (colName in colnames(tcr)[colnames(tcr) != 'barcode']) {
    toAdd <- merged[[colName]]
    names(toAdd) <- merged[[barcode]]
    seuratObject[[colName]] <- toAdd
  }
  
  return(seuratObject)
}

# Expects all_contig_annotations.csv from cellranger vdj    
processAndAggregateTcrClonotypes <- function(clonotypeFile){  
  tcr <- read.table(clonotypeFile, header=T, sep = ',')
  tcr <- tcr[tcr$cdr3 != 'None',]
  
  # drop cellrange '-1' suffix
  tcr$barcode <- gsub("-1", "", tcr$barcode)
  
  #Download named clonotypes and merge:
  # Add clone names:
  labelDf <- labkey.selectRows(
    folderPath="/Labs/Bimber/", 
    schemaName="tcrdb", 
    queryName="clones", 
    showHidden=TRUE,
    colSelect=c('clonename','chain','cdr3','animals'),
    containerFilter=NULL,
    colNameOpt='rname'
  )
  
  labelDf$LabelCol <- paste0(labelDf$clonename)
  
  labelDf <- labelDf %>% 
    group_by(chain, cdr3) %>% 
    summarize(CloneName = paste(sort(unique(LabelCol)), collapse = ","))
  
  tcr <- merge(tcr, labelDf, by.x = c('chain', 'cdr3'), by.y = c('chain', 'cdr3'), all.x = TRUE, all.y = FALSE)
  
  # Many TRDV genes can be used as either alpha or delta TCRs.  10x classifies and TRDV/TRAJ/TRAC clones as 'Multi'.  Re-classify these:
  tcr$chain[tcr$chain == 'Multi' & grepl(pattern = 'TRD', x = tcr$v_gene) & grepl(pattern = 'TRAJ', x = tcr$j_gene) & grepl(pattern = 'TRAC', x = tcr$c_gene)] <- c('TRA')
  
  # Add chain-specific columns:
  tcr$ChainCDR3s <- paste0(tcr$chain, ':', tcr$cdr3)
  for (l in c('TRA', 'TRB', 'TRD', 'TRG')){
    tcr[[l]] <- c(NA)
    tcr[[l]][tcr$chain == l] <- tcr$cdr3[tcr$chain == l]
  }
  
  # Summarize, grouping by barcode
  tcr <- tcr %>% group_by(barcode) %>% summarize(
    ChainCDR3s = paste(sort(unique(ChainCDR3s)), collapse = ","),
    CDR3s = paste(sort(unique(cdr3)), collapse = ","),
    TRA = paste(sort(unique(as.character(TRA))), collapse = ","),
    TRB = paste(sort(unique(as.character(TRB))), collapse = ","),
    TRD = paste(sort(unique(as.character(TRD))), collapse = ","),
    TRG = paste(sort(unique(as.character(TRG))), collapse = ","),
    CloneNames = paste(sort(unique(CloneName)), collapse = ",")  #this is imprecise b/c we count a hit if we match any chain, but this is probably what we often want
  )
  
  tcr$CloneNames <- sapply(strsplit(as.character(tcr$CloneNames), ",", fixed = TRUE), function(x) paste(naturalsort(unique(x)), collapse = ","))
  tcr$CloneNames[tcr$CloneNames == 'NA'] <- NA
  
  tcr$barcode <- as.factor(tcr$barcode)
  for (colName in colnames(tcr)[colnames(tcr) != 'barcode']) {
    v <- tcr[[colName]]
    v <- as.character(v)
    v[v == ''] <- NA
    
    tcr[[colName]] <- as.factor(v)
  }
  
  return(tcr)
}

removeCellCycle <- function(seuratObj) {
  print("Performing cell cycle cleaning ...")
  
  con <- url("https://raw.githubusercontent.com/bbimber/rnaseq/master/data/cellCycle/regev_lab_cell_cycle_genes.txt")
  cc.genes <- readLines(con = con, warn = F)
  close(con)
  if (length(cc.genes) != 97) {
    stop('Something went wrong downloading gene list')
  }
  
  con <- url("https://raw.githubusercontent.com/bbimber/rnaseq/master/data/cellCycle/G2M.txt")
  g2m.genes <- readLines(con =  con)
  close(con)
  
  # We can segregate this list into markers of G2/M phase and markers of S-phase
  s.genes <- cc.genes[1:43]
  g2m.genes <- unique(c(g2m.genes, cc.genes[44:97]))
  
  s.genes <- s.genes[which(s.genes %in% rownames(seuratObj))]
  g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(seuratObj))]
  
  print("Running PCA with cell cycle genes")
  seuratObj <- RunPCA(object = seuratObj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
  print(PCAPlot(seuratObj))
  
  #store values to append later
  SeuratObjsCCPCA <- as.data.frame(SeuratObj@dr$pca@cell.embeddings)
  colnames(SeuratObjsCCPCA) <- paste(colnames(SeuratObjsCCPCA), "CellCycle", sep="_")
  
  seuratObj <- CellCycleScoring(object = seuratObj, 
                                s.genes = s.genes, 
                                g2m.genes = g2m.genes, 
                                set.ident = TRUE)
  
  
  print(PCAPlot(object = SeuratObjs, dim.1 = 1, dim.2 = 2))
  
  print("Regressing out S and G2M score ...")
  seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("S.Score", "G2M.Score"), display.progress = T)
  
  print("Running PCA with variable genes ...")
  seuratObj <- RunPCA(object = seuratObj, pc.genes = VariableFeatures(object = seuratObj), do.print = F)
  
  for (colName in colnames(SeuratObjsCCPCA)) {
    seuratObj[colName] <- SeuratObjsCCPCA[rownames(seuratObj)]  
  }
  
  return(seuratObj)
}