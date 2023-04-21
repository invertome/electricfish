# Load packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Import Salmon output
samples <- metadata$SampleID
files <- file.path("salmon_output", samples, "quant.sf")
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, metadata, ~Tissue + Injection + Feeding)
dds <- DESeq(dds)

# Results
res_deseq2 <- results(dds, contrast = c("Injection", "leptin", "saline"))

# PCA plot
plotPCA(dds)

# Histogram
hist(res_deseq2$pvalue, breaks = 50, col = "skyblue", border = "slateblue", main = "Histogram of p-values")

# Volcano plot
EnhancedVolcano(res_deseq2, title = "Volcano plot")

# MA plot
plotMA(res_deseq2, main = "MA plot")

# Heatmap
top_genes_deseq2 <- rownames(res_deseq2)[order(res_deseq2$pvalue)[1:50]]
heatmap_data_deseq2 <- assay(vsd)[top_genes_deseq2,]
pheatmap(heatmap_data_deseq2, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, main = "Heatmap")
