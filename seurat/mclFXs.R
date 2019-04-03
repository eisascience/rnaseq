# Functions related to Seurat but not directly in pipeline or in development



## Make Combo Seur Obj From Several Pre-processed SerObj.rds 

QuickSerCombObjs <- function(save.fig.path="./Figs", 
                            working.serobjs.path="./data/10X/ser/proc", returnComboObj=F){
  
  require(Matrix)
  require(Seurat)
  
  CleaningLS <- list()
  CleaningLS$SeurObjs <- list()
  getwd()
  

  
  if(!dir.exists(save.fig.path)) dir.create(save.fig.path, recursive = T)
  
  CombObj.path <- paste(working.serobjs.path, "/SeurComboObj.rds", sep="")
  
  if(!file.exists( CombObj.path)){
    
    CleaningLS$all_Seurat_files  <- list.files(paste(getwd(), working.serobjs.path, sep=""), 
                                               full.names = T, 
                                               pattern = "SeuratObj.rds")
    
    CleaningLS$SeurObj.rds_files <-  CleaningLS$all_Seurat_files[grep("SeuratObj.rds", CleaningLS$all_Seurat_files)]
    
    
    
    SerObj.files <- CleaningLS$SeurObj.rds_files#[1:4]
    
    for(fN in 1:length(SerObj.files)){
      # fN=1
      print(fN)
      print("reading in Ser obj")
      exptNum <- gsub("-", "_", gsub("_SeuratObj.rds", "", basename(SerObj.files[fN])))

      seuratObjs <- readRDS(SerObj.files[fN])
      
      seuratObjs[['OrigFileName']] <- c(exptNum)
      
      CleaningLS$SeurObjs[[fN]] <- seuratObjs
    }; remove(seuratObjs)
    
    print("all loaded in now preping for merging..")
    length(CleaningLS$SeurObjs)
    
    names(CleaningLS$SeurObjs) <- paste("ID", as.character(unlist(lapply(SerObj.files, function(xN){
      
      gsub("-", "_", gsub("_SeuratObj.rds", "", basename(xN)))
    }))), sep="_")
    
    
    TempA <- CleaningLS$SeurObjs[[names(CleaningLS$SeurObjs)[1]]]
    TempB <- c(CleaningLS$SeurObjs[[names(CleaningLS$SeurObjs)[1]]])
    
    for(ij in 3:length(CleaningLS$SeurObjs)){
      #ij = 2
      TempB <- append(TempB, CleaningLS$SeurObjs[[ij]])
      
      
    }
    
    
    print("merging ser objs")
    SeurComboObj <- merge(x = TempA,
                          y = TempB, 
                          add.cell.ids = gsub("_", "", gsub("-", "_", gsub("\\.", "-", names(CleaningLS$SeurObjs)))), 
                          do.normalize = F, project = "214")
    
    print("completed merging ser objs, now basic pre-processing")
    SeurComboObj <- NormalizeData(object = SeurComboObj)
    SeurComboObj <- FindVariableFeatures(object = SeurComboObj)
    SeurComboObj <- ScaleData(object = SeurComboObj)
    SeurComboObj <- RunPCA(object = SeurComboObj)
    SeurComboObj <- FindNeighbors(object = SeurComboObj)
    SeurComboObj <- FindClusters(object = SeurComboObj)
    #SeurComboObj <- RunTSNE(object = SeurComboObj, check_duplicates = F, verbose=T)
    #SeurComboObj <- RunUMAP(object = SeurComboObj)
    
    print("saving combo obj")
    saveRDS(SeurComboObj, CombObj.path) 
  } else {
    print("Already Made.... delete or rename to remake")
    if(returnComboObj) SeurComboObj <- readRDS(CombObj.path)
  }
  if(returnComboObj) return(SeurComboObj)
}

