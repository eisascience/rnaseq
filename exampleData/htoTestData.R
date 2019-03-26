
datasets <- list(
  #'267-1' = 'cellHashing/267-1-citeSeqCounts.txt',
  #'267-2' = 'cellHashing/267-2-citeSeqCounts.txt',
  #'253-1' = 'cellHashing/253-1-citeSeqCounts.txt',
  #'260-4' = 'cellHashing/260-4-HTO_cellHashingRawCounts.txt'
  '278-1' = 'cellHashing/278-1-HTO_cellHashingRawCounts.txt',
  '282-1' = 'cellHashing/282-1-HTO_cellHashingRawCounts.txt'
)

expectedHtos <- list(
  '267-1' = c(1:4),
  '267-2' = c(1:4),
  '253-1' = c(1:7),
  '260-4' = c(5:8),
  '278-1' = c(6:9),
  '282-1' = c(1:3, 8, 10, 12)
)

getwd()
outDir <- './outs'
if (!dir.exists(outDir)){
  dir.create(outDir)
}
for (dataset in names(datasets)) {
  # dataset = names(datasets)[1]
  print(paste0('Processing dataset: ', dataset))
  
  barcodeFile <- datasets[[dataset]]
  print(barcodeFile)
  finalCallFile <- paste0('./outs/', dataset, '.calls.txt')
  outputFile <- paste0('../outs/', dataset, '.html')
  barcodeFile <- paste0('../exampleData/', barcodeFile)
  
  
  rmarkdown::render('./exampleData/htoPipeline.rmd', 
                    output_file = outputFile, 
                    output_format = 'html_document')
  
  print('Expected HTOs')
  print(paste0('HTO-', expectedHtos[[dataset]]))
  
  print('Actual HTOs')
  actualHtos <- sapply(rownames(barcodeData), function(x){
    paste0(unlist(strsplit(x, '-'))[1:2], collapse = "-")
  })
  names(actualHtos) <- NULL
  print(actualHtos)
  
  print(paste0('Total cells: ', ncol(barcodeData)))
}
