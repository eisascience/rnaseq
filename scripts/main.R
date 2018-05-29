library(limma)
library(lattice)
library(edgeR)
library(tidyverse)
library(DESeq2)
library(stringr)
library(gplots)
library(ggplot2)
library(reshape)
library(biomaRt)
library(dplyr)
library(WGCNA)
library("BiocParallel")
library(Rlabkey)
library(Matrix)
library(RCurl)

script <- getURL("https://raw.github.com/bbimber/rnaseq/master/scripts/base.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

script <- getURL("https://raw.github.com/bbimber/rnaseq/master/scripts/qc.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

script <- getURL("https://raw.github.com/bbimber/rnaseq/master/scripts/edgeR.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

script <- getURL("https://raw.github.com/bbimber/rnaseq/master/scripts/deseq2.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

script <- getURL("https://raw.github.com/bbimber/rnaseq/master/scripts/wgcna.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#source('https://raw.github.com/bbimber/rnaseq/master/scripts/scde.R')
