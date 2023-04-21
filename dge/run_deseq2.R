if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "tximport", "ggplot2", "pheatmap", "EnhancedVolcano", "limma"))


# Load packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(limma)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Import Salmon output
samples <- metadata$SampleID
files <- file.path("salmon_output", samples, "quant.sf")
txi <- tximport(files, type = "salmon", txOut = TRUE)


# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, metadata, ~Tissue + Injection + Feeding + Tissue:Injection + Tissue:Feeding + Injection:Feeding + Tissue:Injection:Feeding)
dds <- DESeq(dds)

# Define contrasts
contrasts <- list(
  EO_leptin_fooddep_vs_saline_fooddep = c(0, 0, 1, -1, 0, 0, 1, -1),
  EO_leptin_adlib_vs_saline_adlib = c(0, 0, 1, 0, 0, 0, 0, 0),
  SM_leptin_fooddep_vs_saline_fooddep = c(0, 0, 1, -1, 0, 0, 1, 0),
  SM_leptin_adlib_vs_saline_adlib = c(0, 0, 1, 0, 0, 0, 0, 1)
)

# Apply contrasts and obtain results
res_EO_leptin_fooddep_vs_saline_fooddep <- results(dds, contrast = contrasts[["EO_leptin_fooddep_vs_saline_fooddep"]])
res_EO_leptin_adlib_vs_saline_adlib <- results(dds, contrast = contrasts[["EO_leptin_adlib_vs_saline_adlib"]])
res_SM_leptin_fooddep_vs_saline_fooddep <- results(dds, contrast = contrasts[["SM_leptin_fooddep_vs_saline_fooddep"]])
res_SM_leptin_adlib_vs_saline_adlib <- results(dds, contrast = contrasts[["SM_leptin_adlib_vs_saline_adlib"]])

# Save results
dir.create("deseq2_output")
write.csv(res_EO_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/EO_leptin_fooddep_vs_saline_fooddep.csv")
write.csv(res_EO_leptin_adlib_vs_saline_adlib, file = "deseq2_output/EO_leptin_adlib_vs_saline_adlib.csv")
write.csv(res_SM_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/SM_leptin_fooddep_vs_saline_fooddep.csv")
write.csv(res_SM_leptin_adlib_vs_saline_adlib, file = "deseq2_output/SM_leptin_adlib_vs_saline_adlib.csv")

# Save summaries
summary_EO_leptin_fooddep_vs_saline_fooddep <- summary(res_EO_leptin_fooddep_vs_saline_fooddep)
summary_EO_leptin_adlib_vs_saline_adlib <- summary(res_EO_leptin_adlib_vs_saline_adlib)
summary_SM_leptin_fooddep_vs_saline_fooddep <- summary(res_SM_leptin_fooddep_vs_saline_fooddep)
summary_SM_leptin_adlib_vs_saline_adlib <- summary(res_SM_leptin_adlib_vs_saline_adlib)

write.csv(summary_EO_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/summary_EO_leptin_fooddep_vs_saline_fooddep.csv")
write.csv(summary_EO_leptin_adlib_vs_saline_adlib, file = "deseq2_output/summary_EO_leptin_adlib_vs_saline_adlib.csv")
write.csv(summary_SM_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/summary_SM_leptin_fooddep_vs_saline_fooddep.csv")
write.csv(summary_SM_leptin_adlib_vs_saline_adlib, file = "deseq2_output/summary_SM_leptin_adlib_vs_saline_adlib.csv")

# PCA plot
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("Tissue", "Injection", "Feeding"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = Tissue, shape = Injection)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
print(pca_plot)
ggsave("deseq2_output/pca_plot.png", plot = pca_plot)

# Histograms
hist_EO_leptin_fooddep_vs_saline_fooddep <- ggplot(data.frame(x=res_EO_leptin_fooddep_vs_saline_fooddep$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), col = "slateblue", fill = "skyblue") +
  labs(title = "Histogram of p-values (EO Leptin Fooddep vs Saline Fooddep)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_EO_leptin_fooddep_vs_saline_fooddep.png", plot = hist_EO_leptin_fooddep_vs_saline_fooddep)

hist_EO_leptin_adlib_vs_saline_adlib <- ggplot(data.frame(x=res_EO_leptin_adlib_vs_saline_adlib$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), col = "slateblue", fill = "skyblue") +
  labs(title = "Histogram of p-values (EO Leptin Adlib vs Saline Adlib)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_EO_leptin_adlib_vs_saline_adlib.png", plot = hist_EO_leptin_adlib_vs_saline_adlib)

hist_SM_leptin_fooddep_vs_saline_fooddep <- ggplot(data.frame(x=res_SM_leptin_fooddep_vs_saline_fooddep$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), col = "slateblue", fill = "skyblue") +
  labs(title = "Histogram of p-values (SM Leptin Fooddep vs Saline Fooddep)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_SM_leptin_fooddep_vs_saline_fooddep.png", plot = hist_SM_leptin_fooddep_vs_saline_fooddep)

hist_SM_leptin_adlib_vs_saline_adlib <- ggplot(data.frame(x=res_SM_leptin_adlib_vs_saline_adlib$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), col = "slateblue", fill = "skyblue") +
  labs(title = "Histogram of p-values (SM Leptin Adlib vs Saline Adlib)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_SM_leptin_adlib_vs_saline_adlib.png", plot = hist_SM_leptin_adlib_vs_saline_adlib)



# Volcano plots
vol_EO_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_EO_leptin_fooddep_vs_saline_fooddep, title = "Volcano plot (EO Leptin Fooddep vs Saline Fooddep)", pCutoff = 0.05, FCcutoff = 1)
ggsave("deseq2_output/volcano_EO_leptin_fooddep_vs_saline_fooddep.png", plot = vol_EO_leptin_fooddep_vs_saline_fooddep)

vol_EO_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_EO_leptin_adlib_vs_saline_adlib, title = "Volcano plot (EO Leptin Adlib vs Saline Adlib)", pCutoff = 0.05, FCcutoff = 1)
ggsave("deseq2_output/volcano_EO_leptin_adlib_vs_saline_adlib.png", plot = vol_EO_leptin_adlib_vs_saline_adlib)

vol_SM_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_SM_leptin_fooddep_vs_saline_fooddep, title = "Volcano plot (SM Leptin Fooddep vs Saline Fooddep)", pCutoff = 0.05, FCcutoff = 1)
ggsave("deseq2_output/volcano_SM_leptin_fooddep_vs_saline_fooddep.png", plot = vol_SM_leptin_fooddep_vs_saline_fooddep)

vol_SM_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_SM_leptin_adlib_vs_saline_adlib, title = "Volcano plot (SM Leptin Adlib vs Saline Adlib)", pCutoff = 0.05, FCcutoff = 1)
ggsave("deseq2_output/volcano_SM_leptin_adlib_vs_saline_adlib.png", plot = vol_SM_leptin_adlib_vs_saline_adlib)

# MA plots
ma_EO_leptin_fooddep_vs_saline_fooddep <- plotMA(res_EO_leptin_fooddep_vs_saline_fooddep, main = "MA plot (EO Leptin Fooddep vs Saline Fooddep)")
ggsave("deseq2_output/ma_EO_leptin_fooddep_vs_saline_fooddep.png", plot = ma_EO_leptin_fooddep_vs_saline_fooddep)

ma_EO_leptin_adlib_vs_saline_adlib <- plotMA(res_EO_leptin_adlib_vs_saline_adlib, main = "MA plot (EO Leptin Adlib vs Saline Adlib)")
ggsave("deseq2_output/ma_EO_leptin_adlib_vs_saline_adlib.png", plot = ma_EO_leptin_adlib_vs_saline_adlib)

ma_SM_leptin_fooddep_vs_saline_fooddep <- plotMA(res_SM_leptin_fooddep_vs_saline_fooddep, main = "MA plot (SM Leptin Fooddep vs Saline Fooddep)")
ggsave("deseq2_output/ma_SM_leptin_fooddep_vs_saline_fooddep.png", plot = ma_SM_leptin_fooddep_vs_saline_fooddep)

ma_SM_leptin_adlib_vs_saline_adlib <- plotMA(res_SM_leptin_adlib_vs_saline_adlib, main = "MA plot (SM Leptin Adlib vs Saline Adlib)")
ggsave("deseq2_output/ma_SM_leptin_adlib_vs_saline_adlib.png", plot = ma_SM_leptin_adlib_vs_saline_adlib)

# Heatmaps
log2FoldChange = lfcShrink(dds, coef = "Injection_leptin_vs_saline")
top_genes <- head(order(log2FoldChange), 50)
heatmap_data <- assay(dds)[ top_genes, ]
pheatmap(heatmap_data, cluster_cols = TRUE, cluster_rows = TRUE)
ggsave("deseq2_output/heatmap_top_genes.png")

# Print summary
cat("Generated files for DESeq2 analysis:\n")
cat(paste("deseq2_output/EO_leptin_fooddep_vs_saline_fooddep.csv\n",
          "deseq2_output/EO_leptin_adlib_vs_saline_adlib.csv\n",
          "deseq2_output/SM_leptin_fooddep_vs_saline_fooddep.csv\n",
          "deseq2_output/SM_leptin_adlib_vs_saline_adlib.csv\n",
          "deseq2_output/summary_EO_leptin_fooddep_vs_saline_fooddep.csv\n",
          "deseq2_output/summary_EO_leptin_adlib_vs_saline_adlib.csv\n",
          "deseq2_output/summary_SM_leptin_fooddep_vs_saline_fooddep.csv\n",
          "deseq2_output/summary_SM_leptin_adlib_vs_saline_adlib.csv\n",
          "deseq2_output/pca_plot.png\n",
          "deseq2_output/histogram_EO_leptin_fooddep_vs_saline_fooddep.png\n",
          "deseq2_output/histogram_EO_leptin_adlib_vs_saline_adlib.png\n",
          "deseq2_output/histogram_SM_leptin_fooddep_vs_saline_fooddep.png\n",
          "deseq2_output/histogram_SM_leptin_adlib_vs_saline_adlib.png\n",
          "deseq2_output/volcano_EO_leptin_fooddep_vs_saline_fooddep.png\n",
          "deseq2_output/volcano_EO_leptin_adlib_vs_saline_adlib.png\n",
          "deseq2_output/volcano_SM_leptin_fooddep_vs_saline_fooddep.png\n",
          "deseq2_output/volcano_SM_leptin_adlib_vs_saline_adlib.png\n",
          "deseq2_output/ma_EO_leptin_fooddep_vs_saline_fooddep.png\n",
          "deseq2_output/ma_EO_leptin_adlib_vs_saline_adlib.png\n",
          "deseq2_output/ma_SM_leptin_fooddep_vs_saline_fooddep.png\n",
          "deseq2_output/ma_SM_leptin_adlib_vs_saline_adlib.png\n",
          "deseq2_output/heatmap_top_genes.png\n",
          sep = ""))


