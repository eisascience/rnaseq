library(Seurat)
library(Rlabkey)
library(data.table)
library(dplyr)
library(naturalsort)
library(DropletUtils)
library(Matrix)

labkey.setDefaults(baseUrl = "https://prime-seq.ohsu.edu")

readAndFilter10xData <- function(dataDir, datasetName) {
  seuratRawData <- Read10X(data.dir = dataDir)
  seuratRawData <- performEmptyDropletFiltering(seuratRawData)
  
  seuratObj <- createSeuratObj(seuratRawData, project = datasetName)
  printQcPlots1(seuratObj)

  return(seuratObj)  
}

createSeuratObj <- function(seuratData = NA, project = NA, minFeatures = 25, minCells = 0){
  seuratObj <- CreateSeuratObject(counts = seuratData, min.cells = minCells, min.features = minFeatures, project = project)
  
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
  print(paste0('Input cells: ', nrow(e.out)))
  e.out <- e.out[!is.na(e.out$LogProb),]
  e.out$is.cell <- e.out$FDR <= fdrThreshold
  print(paste0('Passing cells: ', sum(e.out$is.cell, na.rm=TRUE)))
  print(paste0('Failing cells: ', sum(!e.out$is.cell, na.rm=TRUE)))
  
  #If there are any entries with FDR above the desired threshold and Limited==TRUE, it indicates that npts should be increased in the emptyDrops call.
  print(table(Limited=e.out$Limited, Significant=e.out$is.cell))
  
  toPlot <- e.out[is.finite(e.out$LogProb),]
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
  seuratObj <- NULL
  for (exptNum in names(data)) {
    print(paste0('adding expt: ', exptNum))
    prefix <- paste0(exptNum)
    seuratObjs[[exptNum]] <- RenameCells(object = seuratObjs[[exptNum]], add.cell.id = prefix)
    seuratObjs[[exptNum]][['BarcodePrefix']] <- c(prefix)
    
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

processSeurat1 <- function(seuratObj, saveFile = NULL, doCellCycle = T, doCellFilter = F,
                           nUMI.high = 20000, nGene.high = 3000, pMito.high = 0.15,
                           nUMI.low = 0.99, nGene.low = 200, pMito.low = -Inf){
  
  if(doCellFilter & !hasStepRun(seuratObj, 'FilterCells')) {
    print("Filtering Cells...")
    seuratObj@misc$OriginalCells <- length(colnames(x = seuratObj))
    seuratObj <- subset(x = seuratObj, subset = nCount_RNA > nGene.low & nCount_RNA < nGene.high)
    seuratObj <- subset(x = seuratObj, subset = nFeature_RNA > nUMI.low & nFeature_RNA < nUMI.high)
    seuratObj <- subset(x = seuratObj, subset = p.mito > pMito.low & p.mito < pMito.high)
    
    print(paste0('Initial cells: ', seuratObj@misc$OriginalCells, ', after filter: ', length(colnames(x = seuratObj))))
    
    seuratObj <- markStepRun(seuratObj, 'FilterCells')
  }
  
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
  
  if (!hasStepRun(seuratObj, 'JackStraw')) {
    seuratObj <- JackStraw(object = seuratObj, num.replicate = 100)
    seuratObj <- markStepRun(seuratObj, 'JackStraw', saveFile)
  }
  
  if (!hasStepRun(seuratObj, 'ScoreJackStraw')) {
    seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
    seuratObj <- markStepRun(seuratObj, 'ScoreJackStraw')
  }
  
  print(paste0('Variable genes: ', length(x = VariableFeatures(object = seuratObj))))
  
  print(VizDimLoadings(object = seuratObj, dims = 1:2))
  print(DimPlot(object = seuratObj))
  
  print(DimHeatmap(object = seuratObj, dims = 1, cells = 500, balanced = TRUE))
  print(DimHeatmap(object = seuratObj, dims = 1:20, cells = 500, balanced = TRUE))
  
  print(JackStrawPlot(object = seuratObj, dims = 1:20))
  print(ElbowPlot(object = seuratObj))
  
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
  g2m.genes <- readLines(con =  con, warn = F)
  close(con)
  
  # We can segregate this list into markers of G2/M phase and markers of S-phase
  s.genes <- cc.genes[1:43]
  g2m.genes <- unique(c(g2m.genes, cc.genes[44:97]))
  
  s.genes <- s.genes[which(s.genes %in% rownames(seuratObj))]
  g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(seuratObj))]
  
  print("Running PCA with cell cycle genes")
  seuratObj <- RunPCA(object = seuratObj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
  print(DimPlot(object = seuratObj, reduction = "pca"))
  
  #store values to append later
  SeuratObjsCCPCA <- as.data.frame(seuratObj@reductions$pca@cell.embeddings)
  colnames(SeuratObjsCCPCA) <- paste(colnames(SeuratObjsCCPCA), "CellCycle", sep="_")
  
  seuratObj <- CellCycleScoring(object = seuratObj, 
                                s.features = s.genes, 
                                g2m.features = g2m.genes, 
                                set.ident = TRUE)
  
  print(DimPlot(object = seuratObj, reduction = "pca"))
  
  print("Regressing out S and G2M score ...")
  seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("S.Score", "G2M.Score"), display.progress = T)
  
  print("Running PCA with variable genes ...")
  seuratObj <- RunPCA(object = seuratObj, pc.genes = VariableFeatures(object = seuratObj), do.print = F)
  
  for (colName in colnames(SeuratObjsCCPCA)) {
    seuratObj[[colName]] <- SeuratObjsCCPCA[colnames(seuratObj),colName]  
  }
  
  return(seuratObj)
}

findClustersAndDimRedux <- function(seuratObj, dimsToUse, saveFile = NULL) {
  if (!hasStepRun(seuratObj, 'FindNeighbors')) {
    seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
    seuratObj <- markStepRun(seuratObj, 'FindNeighbors')
  }
  
  if (!hasStepRun(seuratObj, 'FindClusters')) {
    for (resolution in c(0.2, 0.4, 0.8, 1.2, 0.6)){
      seuratObj <- FindClusters(object = seuratObj, resolution = resolution)
      seuratObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = seuratObj)
      seuratObj <- markStepRun(seuratObj, 'FindClusters', saveFile)
    }
  }
  
  if (!hasStepRun(seuratObj, 'RunTSNE')) {
    #See: https://github.com/satijalab/seurat/issues/167
    seuratObj <- RunTSNE(object = seuratObj, dims.use = dimsToUse, check_duplicates = FALSE)
    seuratObj <- markStepRun(seuratObj, 'RunTSNE', saveFile)
  }

  if (!hasStepRun(seuratObj, 'RunUMAP')) {
    seuratObj <- RunUMAP(seuratObj,
                         dims = dimsToUse,
                         n.neighbors = 40L,
                         min.dist = 0.2,
                         metric = "correlation",
                         seed.use = 1234)
    seuratObj <- markStepRun(seuratObj, 'RunUMAP', saveFile)
  }
  
  for (reduction in c('tsne', 'umap')){
    plot1 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_0.2", label = TRUE) + ggtitle('Resolution: 0.2')
    plot2 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_0.4", label = TRUE) + ggtitle('Resolution: 0.4')
    plot3 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_0.6", label = TRUE) + ggtitle('Resolution: 0.6')
    plot4 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_0.8", label = TRUE) + ggtitle('Resolution: 0.8')
    plot5 <- DimPlot(object = seuratObj, reduction = reduction, group.by = "ClusterNames_1.2", label = TRUE) + ggtitle('Resolution: 1.2')
    
    print(CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5), legend = 'none'))
    
    #print(DimPlot(object = seuratObj, reduction = reduction, group.by = "BarcodePrefix", label = TRUE) + ggtitle('Dataset'))
  }
  

  return(seuratObj)
}

findMarkers <- function(seuratObj, resolutionToUse, outFile, saveFileMarkers = NULL) {
  Idents(seuratObj) <- seuratObj[[paste0('ClusterNames_',resolutionToUse)]]
  
  if (file.exists(saveFileMarkers)) {
    print('resuming from file')
    seuratObj.markers <- readRDS(saveFileMarkers)
  } else {
    seuratObj.markers <- FindAllMarkers(object = seuratObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    saveRDS(seuratObj.markers, file = saveFileMarkers)
  }
  
  toPlot <- seuratObj.markers %>% filter(p_val_adj < 0.001) %>% group_by(cluster)  %>% filter(avg_logFC > 0.5) %>% top_n(9, avg_logFC) %>% select(gene)
  
  write.table(toPlot, file = outFile, sep = '\t', row.names = F, quote = F)
  
  print(DimPlot(object = seuratObj, reduction = 'tsne'))
  
  top10 <- seuratObj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
  print(DoHeatmap(object = seuratObj, features = top10$gene))
  
  return(toPlot)
}

createExampleData <- function(nRow = 100, nCol = 10){
  my.counts <- matrix(rpois(1000, lambda=5), ncol=nCol, nrow=nRow)
  my.counts <- as(my.counts, "dgCMatrix")
  cell.ids <- paste0("BARCODE-", seq_len(ncol(my.counts)))
  
  ngenes <- nrow(my.counts)
  gene.ids <- paste0("ENSG0000", seq_len(ngenes))
  gene.symb <- paste0("GENE", seq_len(ngenes))
  
  # Writing this to file:
  tmpdir <- tempfile()
  write10xCounts(tmpdir, my.counts, gene.id=gene.ids, 
                 gene.symbol=gene.symb, barcodes=cell.ids)
  return(tmpdir)
}
                           
                           
                           
AddModuleScore_SERIII <- function(
  #using Seurat v2, AddModuleScore and converting to v3
  #this worked in v2 but their new version 3 breaks so this replaces it
  
  object,
  genes.list = NULL,
  genes.pool = NULL,
  n.bin = 25,
  seed.use = 1,
  ctrl.size = 100,
  use.k = FALSE,
  enrich.name = "Cluster",
  random.seed = 1
) {
  set.seed(seed = random.seed)
  genes.old <- genes.list
  if (use.k) {
    genes.list <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = genes.list)
  } else {
    if (is.null(x = genes.list)) {
      stop("Missing input gene list")
    }
    genes.list <- lapply(
      X = genes.list,
      FUN = function(x) {
        #return(intersect(x = x, y = rownames(x = object@data)))
        return(intersect(x = x, y = rownames(object)))
      }
    )
    cluster.length <- length(x = genes.list)
  }
  if (!all(LengthCheck(values = genes.list))) {
    warning(paste(
      'Could not find enough genes in the object from the following gene lists:',
      paste(names(x = which(x = ! LengthCheck(values = genes.list)))),
      'Attempting to match case...'
    ))
    genes.list <- lapply(
      X = genes.old,
      FUN = CaseMatch, match = rownames(x = object@data)
    )
  }
  if (!all(LengthCheck(values = genes.list))) {
    stop(paste(
      'The following gene lists do not have enough genes present in the object:',
      paste(names(x = which(x = ! LengthCheck(values = genes.list)))),
      'exiting...'
    ))
  }
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(object)
  }
  data.avg <- Matrix::rowMeans(object@assays$RNA@data[genes.pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(
    x = data.avg,
    m = round(x = length(x = data.avg) / n.bin)
  ))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(x = genes.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[genes.use[j]])],
          size = ctrl.size,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object@assays$RNA@data)
  )

  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = object@assays$RNA@data[genes.use, ])
  }
  genes.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object@assays$RNA@data)
  )
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    data.use <- object@assays$RNA@data[genes.use, , drop = FALSE]
    genes.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
  rownames(x = genes.scores.use) <- colnames(x = object@assays$RNA@data)
  object <- AddMetaData(
    object = object,
    metadata = genes.scores.use,
    col.name = colnames(x = genes.scores.use)
  )
  gc(verbose = FALSE)
  return(object)
}
