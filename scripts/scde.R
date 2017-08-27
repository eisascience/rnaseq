library(DESeq2)
library("BiocParallel")
register(MulticoreParam(8))

counts <- geneCountMatrix

nGenes <- length(counts[,1])
coverage <- colSums(counts)/nGenes

counts <- counts[,coverage>1]
cell.labels <- colnames(counts)
coverage <- coverage[coverage>1]
nCells <- length(cell.labels)

counts.norm <- t(apply(counts,1,function(x) x/coverage)) # simple normalization method
top.genes <- tail(order(rowSums(counts.norm)),10)
expression <- log2(counts.norm[top.genes,]+1) # add a pseudocount of 1



means <- apply(counts.norm,1,mean)
excess.var <- apply(counts,1,var)-means
excess.var[excess.var < 0] <- NA
overdispersion <- excess.var / means^2

png('hist.png')
hist(log2(overdispersion),main="Variance of read counts is higher than Poisson")
dev.off()

groups <- factor(metaUnfilter$Peptide[coverage>1])
ord <- order(groups)

dds <- DESeqDataSetFromMatrix(counts,DataFrame(groups), ~groups)

dds <- DESeq(dds,parallel=TRUE)
res <- results(dds)

find.significant.genes <- function(de.result,alpha=0.05) {
  # filter out significant genes based on FDR adjusted p-values
  filtered <- de.result[(de.result$padj < alpha) & !is.infinite(de.result$log2FoldChange) & !is.nan(de.result$log2FoldChange) & !is.na(de.result$padj),]
  # order by p-value, and print out only the gene name, mean count, and log2 fold change
  sorted <- filtered[order(filtered$padj),c(1,2,6)]
}

de2.genes <- find.significant.genes(res)



library(scde)
n.cores <- 2

scde.fitted.model <- scde.error.models(counts=counts,groups=groups,n.cores=n.cores,save.model.plots=F)
save(scde.fitted.model,file="scde_fit.RData")

scde.prior <- scde.expression.prior(models=scde.fitted.model,counts=counts)

ediff <- scde.expression.difference(scde.fitted.model,counts,scde.prior,groups=groups,n.cores=n.cores)
p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
significant.genes <- which(p.values.adj<0.05)
length(significant.genes)

ord <- order(p.values.adj[significant.genes]) # order by p-value
de <- cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
colnames(de) <- c("Lower bound","log2 fold change","Upper bound","p-value")

de[1:15,]


#plot:
#scde.test.gene.expression.difference("Tdh",models=scde.fitted.model,counts=counts,prior=scde.prior)

save.image('checkpoint.rdata')