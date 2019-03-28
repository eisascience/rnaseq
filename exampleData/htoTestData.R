
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

outDir <- './outs/'
if (!dir.exists(outDir)){
  dir.create(outDir)
}

for (dataset in names(datasets)) {
  print(paste0('Processing dataset: ', dataset))
  
  barcodeFile <- datasets[[dataset]]
  print(barcodeFile)
  finalCallFile <- paste0(outDir, dataset, '.calls.txt')
  outputFile <- paste0(outDir, dataset, '.html')
  barcodeFile <- paste0('../exampleData/', barcodeFile)
  
  
  rmarkdown::render('../exampleData/htoPipeline.rmd', 
                    output_file = outputFile, 
                    output_format = 'html_document')
  
  datasetExpectations <- expectations[[dataset]]
  
  print('Expected HTOs')
  expectedHtos <- sort(paste0('HTO-', datasetExpectations$htos))
  print(expectedHtos)
  
  print('Actual HTOs')
  actualHtos <- as.character(unique(dt$HTO))
  actualHtos <- sort(actualHtos[!(actualHtos %in% c('Negative', 'Doublet'))])
  print(actualHtos)
  
  if (!all(expectedHtos == actualHtos)) {
    print('HTO set did not match!')
  }
  
  print(paste0('Total cells after filter: ', ncol(barcodeData)))
}
