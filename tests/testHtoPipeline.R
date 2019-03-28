
datasets <- list(
  '267-1' = 'cellHashing/267-1-citeSeqCounts.txt',
  '267-2' = 'cellHashing/267-2-citeSeqCounts.txt',
  '253-1' = 'cellHashing/253-1-citeSeqCounts.txt',
  '260-4' = 'cellHashing/260-4-HTO_cellHashingRawCounts.txt',
  '278-1' = 'cellHashing/278-1-HTO_cellHashingRawCounts.txt',
  '282-1' = 'cellHashing/282-1-HTO_cellHashingRawCounts.txt'
)

expectations <- list(
  '267-1' = list(htos = c(1:4), gexBarcodeFile = ''),
  '267-2' = list(htos = c(1:4), gexBarcodeFile = ''),
  '253-1' = list(htos = c(1:7), gexBarcodeFile = ''),
  '260-4' = list(htos = c(5:8), gexBarcodeFile = ''),
  '278-1' = list(htos = c(6:9), gexBarcodeFile = ''),
  '282-1' = list(htos = c(1:3, 8, 10, 12), gexBarcodeFile = '')
)

#relative to ./tests (assumed to be current working dir)
outDir <- './outs/'
if (!dir.exists(outDir)){
  dir.create(outDir)
}

#now adjust for relative location to RMD file:
outDirRmd <- '../tests/' + outDir

outcomes <- data.frame(file = character(), expected = character(), matrixActual = character(), callActual = character(), diff = character())
for (dataset in names(datasets)) {
  print(paste0('Processing dataset: ', dataset))
  
  barcodeFile <- datasets[[dataset]]
  print(barcodeFile)
  finalCallFile <- paste0(outDirRmd, dataset, '.calls.txt')
  outputFile <- paste0(outDirRmd, dataset, '.html')
  barcodeFile <- paste0('../exampleData/', barcodeFile)
  
  #Note: paths are relative to the RMD file
  rmarkdown::render('../exampleData/htoPipeline.rmd', output_file = outputFile, output_format = 'html_document')
  
  datasetExpectations <- expectations[[dataset]]
  
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

write.table(outcomes, file = 'testResults.txt', quote = F, sep = '\t', row.names = F)