#### BiocManager packages install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#### Load libraries
library(BiocManager)
library(DESeq2)
library(ggplot2)

#### Load data
data = read.csv('inputfile.csv', sep = ',', header = T, row.names = 1)

treatments = DataFrame(condition=factor(c('t1', 't1','t2', 't2', 't3', 't3', 't4', 't4')))  
dds = DESeqDataSetFromMatrix(data, treatments, formula(~condition))

#### Filtering out genes with less than 10 reads in at least 8 samples
filter = rowSums(counts(dds)>=10)>=8
ddsCountFlt = dds[filter,]

#### Differential expression analysis
ddsCountFlt = DESeq(ddsCountFlt)

normCounts = counts(ddsCountFlt, normalized = TRUE)
write.csv(normCounts, file = 'norm_counts.csv')
readsPerSamp_Norm = colSums(normCounts)

res = results(ddsCountFlt, contrast = c('condition', 't1', 't2'))
write.csv(as.data.frame(res), file = 'resultst1_t2.csv')
mcols(res, use.names = TRUE)

resOrdered = res[order(res$padj),]
sigLog2p005 = resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.05 & abs(resOrdered$log2FoldChange) >= 2,]
write.csv(sigLog2p005, file = 'significant_FC2_p0.05.csv')
