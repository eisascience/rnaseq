library(Rlabkey)
library(data.table)

labkey.setDefaults(baseUrl = "https://prime-seq.ohsu.edu")

getCdnaId <- function(readsetId, type = 'GEX'){
  df <- getCdnaRecord(readsetId, type)
  
  if (nrow(df) > 1) {
    
  }
}


getCdnaRecords <- function(readsetId, type = 'GEX') {
  fieldName <- switch(type, 
    'GEX' = 'readsetId',
    'TCR' = 'enrichedReadsetId',
    'HTO' = 'hashingReadsetId'
  )
  
  df <- labkey.selectRows(
    folderPath="/Labs/Bimber/", 
    schemaName="tcrdb", 
    queryName="clones", 
    showHidden=TRUE,
    #colSelect=c(''),
    colFilter=makeFilter(c(fieldName, "EQUAL", readsetId)), 
    containerFilter=NULL,
    colNameOpt='rname'
  )
  
  return(df)
}