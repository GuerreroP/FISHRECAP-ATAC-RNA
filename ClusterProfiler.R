##### Packages installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

##### Load libraries
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library(enrichplot)

##### Load input data and 
gene <- read.csv("input.csv", header=TRUE, stringsAsFactors = FALSE, sep=",")
TERM2GENE_BP <- read.csv("TERM2GENE_BP.csv", sep=";")
TERM2NAME_BP <- read.csv("TERM2NAME_BP.csv", sep=";")

#### Gene ontology enrichment
result <- enricher(
  gene$X,
  pvalueCutoff = 5,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 5,
  TERM2GENE=TERM2GENE_BP,
  TERM2NAME=TERM2NAME_BP
)

#### Save and plot results
write.csv(result@result, "Final_result.csv")

pdf("file.pdf", paper = "a4r", width = 1920)
plot <- barplot(result, title="GO (BP)", showCategory=10, x="GeneRatio", font.size = 12) 
plot +
  theme(
    panel.grid = element_blank()
  ) + scale_color_brewer(palette = "Blues")
dev.off()
