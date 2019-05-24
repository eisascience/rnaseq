
### Read in Seurat objects, pre-processed (or not, just note to set scale = T/F for FindIntegrationAnchors())


SerProcFiles <- list.files(paste0(counts.path, "SerProc"), full.names = T)

SerProcFiles
SerProcFiles <- SerProcFiles[c(2,3,4,5,6,7,10,11,12, 15, 16, 17, 18,19,20,23,24,25)] #remove bulk
SerProcFilesNames <- unlist(lapply(SerProcFiles, function(SerXf){
  # SerXf = SerProcFiles [1]
  tempName <- gsub("_", "", gsub("-", "_", gsub("\\.", "", gsub(".seurat.rds", "", basename(SerXf)))))
  
}))

if(!file.exists(paste0(save.path, "/SerProcLS.rds"))){
  SerProcLS <- unlist(lapply(1:length(SerProcFiles), function(SerXfN){
  # SerXfN = 1
  print(SerXfN)
  tempSeur <- readRDS(SerProcFiles[SerXfN])
  tempSeur@meta.data$RowID <- rownames(tempSeur@meta.data)
  tempSeur <- downloadAndAppendTcrClonotypes(tempSeur, outPath = './')
  tempSeur$CondX <- SerProcFilesNames[SerXfN]
  tempSeur
}))

names(SerProcLS) <- paste0("ID", SerProcFilesNames)
  saveRDS(SerProcLS, paste0(save.path, "/SerProcLS.rds"))
} else {
  SerProcLS <- readRDS(paste0(save.path, "/SerProcLS.rds"))
}


## Find Anchors


if(!file.exists(paste0(save.path, "/immune.anchors.rds"))){
  immune.anchors <- FindIntegrationAnchors(object.list = 
                                             SerProcLS, 
                                           dims = 1:20, scale = F)
  saveRDS(immune.anchors, paste0(save.path, "/immune.anchors.rds"))
} else {
  immune.anchors <- readRDS(paste0(save.path, "/immune.anchors.rds"))
}


if(!file.exists(paste0(save.path, "/immune.combined_base.rds"))){
  immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
  
  DefaultAssay(immune.combined) <- "integrated"

  # Run the standard workflow for visualization and clustering
  immune.combined <- ScaleData(immune.combined, verbose = T)
  
 
  
   immune.combined <- FindVariableFeatures(object = immune.combined, selection.method="vst")

  
  immune.combined <- RunPCA(immune.combined, npcs = 35, verbose = T, features = )
  
  ElbowPlot(immune.combined)
  

  # t-SNE and Clustering
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
  immune.combined <- FindClusters(immune.combined, resolution = 0.5)

  saveRDS(immune.combined, paste0(save.path, "/immune.combined_base.rds"))
  
} else {
  immune.combined <- readRDS(paste0(save.path, "/immune.combined_base.rds"))
}
