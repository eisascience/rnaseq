
datasets <- list(
  '267-1' = 'cellHashing/267-1-citeSeqCounts.txt',
  '267-2' = 'cellHashing/267-2-citeSeqCounts.txt',
  '253-1' = 'cellHashing/253-1-citeSeqCounts.txt',
  '260-4' = 'cellHashing/260-4-HTO_cellHashingRawCounts.txt',
  '278-1' = 'cellHashing/278-1-HTO_cellHashingRawCounts.txt',
  '282-1' = 'cellHashing/282-1-HTO_cellHashingRawCounts.txt'
)

expectations <- list(
  '267-1' = list(htos = c(1:4), gexBarcodeFile = NULL),
  '267-2' = list(htos = c(1:4), gexBarcodeFile = NULL),
  '253-1' = list(htos = c(1:7), gexBarcodeFile = NULL),
  '260-4' = list(htos = c(5:8), gexBarcodeFile = NULL),
  '278-1' = list(htos = c(6:9), gexBarcodeFile = NULL),
  '282-1' = list(htos = c(1:3, 8, 10, 12), gexBarcodeFile = 'cellHashing/282-1-whitelist.txt')
)

#relative to ./tests (assumed to be current working dir)
outDir <- './outs/'
if (!dir.exists(outDir)){
  dir.create(outDir)
}

#now adjust for relative location to RMD file:
outDirRmd <- paste0('../tests/', outDir)

outcomes <- data.frame(file = character(), expected = character(), matrixActual = character(), callActual = character(), diff = character())
for (dataset in names(datasets)) {
  print(paste0('Processing dataset: ', dataset))
  
  barcodeFile <- datasets[[dataset]]
  print(barcodeFile)
  finalCallFile <- paste0(outDirRmd, dataset, '.calls.txt')
  outputFile <- paste0(outDirRmd, dataset, '.html')
  barcodeFile <- paste0('../exampleData/', barcodeFile)
  
  datasetExpectations <- expectations[[dataset]]
  
  #this will result in a concordance report being created
  if (!is.null(datasetExpectations[['gexBarcodeFile']])){
    whitelistFile <- datasetExpectations[['gexBarcodeFile']]   
    summaryFile <- paste0(outDirRmd, dataset, '.summary.txt')
  }
  
  #Note: paths are relative to the RMD file
  rmarkdown::render('../exampleData/htoPipeline.rmd', output_file = outputFile, output_format = 'html_document')
  
  if (exists('whitelistFile')) {
    rm(whitelistFile)
    rm(summaryFile)
  }
  
  print('Expected HTOs')
  expectedHtos <- sort(paste0('HTO-', datasetExpectations$htos))
  print(expectedHtos)
  
  print('Actual HTOs (barcode matrix)')
  actualHtosMatrix <- sort(unname(simplifyHtoNames(rownames(barcodeData))))

  print('Actual HTOs (final calls)')
  actualHtos <- as.character(unique(dt$HTO))
  actualHtos <- sort(actualHtos[!(actualHtos %in% c('Negative', 'Doublet'))])
  print(actualHtos)
  
  differing <- c(setdiff(expectedHtos, actualHtosMatrix), setdiff(actualHtosMatrix, expectedHtos))
  if (length(differing) > 0 ) {
    print('HTO set did not match!')
    print(differing)
  }
  
  print(paste0('Total cells after filter: ', ncol(barcodeData)))
  
  outcomes <- rbind(outcomes, data.frame(file = barcodeFile, 
                                         expected = paste0(expectedHtos, collapse = ','), 
                                         matrixActual = paste0(actualHtosMatrix, collapse = ','), 
                                         callActual = paste0(actualHtos, collapse = ','),
                                         diff = paste0(differing, collapse = ',')))
}

write.table(outcomes, file = paste0(outDir,'testResults.txt'), quote = F, sep = '\t', row.names = F)