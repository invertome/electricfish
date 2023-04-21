if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "tximport", "ggplot2", "pheatmap", "EnhancedVolcano"))


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
res_EO_leptin_fooddep_vs_saline_fooddep <- results(dds, contrast = list("Injection_leptin_vs_saline", "Feeding_fooddep"), listValues = c(1,-1), test = "LRT")
res_EO_leptin_adlib_vs_saline_adlib <- results(dds, contrast = list("Injection_leptin_vs_saline", "Feeding_adlib"), listValues = c(1,-1), test = "LRT")
res_SM_leptin_fooddep_vs_saline_fooddep <- results(dds, contrast = list("Injection_leptin_vs_saline", "Feeding_fooddep"), listValues = c(-1,1), test = "LRT")
res_SM_leptin_adlib_vs_saline_adlib <- results(dds, contrast = list("Injection_leptin_vs_saline", "Feeding_adlib"), listValues = c(-1,1), test = "LRT")

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
plotPCA(dds)
ggsave("deseq2_output/pca_plot.png")

# Histograms
hist(res_EO_leptin_fooddep_vs_saline_fooddep$pvalue, breaks = 50, col = "skyblue", border = "slateblue", main = "Histogram of p-values (EO Leptin Fooddep vs Saline Fooddep)")
ggsave("deseq2_output/histogram_EO_leptin_fooddep_vs_saline_fooddep.png")
hist(res_EO_leptin_adlib_vs_saline_adlib$pvalue, breaks = 50, col = "skyblue", border = "slateblue", main = "Histogram of p-values (EO Leptin Adlib vs Saline Adlib)")
ggsave("deseq2_output/histogram_EO_leptin_adlib_vs_saline_adlib.png")
hist(res_SM_leptin_fooddep_vs_saline_fooddep$pvalue, breaks = 50, col = "skyblue", border = "slateblue", main = "Histogram of p-values (SM Leptin Fooddep vs Saline Fooddep)")
ggsave("deseq2_output/histogram_SM_leptin_fooddep_vs_saline_fooddep.png")
hist(res_SM_leptin_adlib_vs_saline_adlib$pvalue, breaks = 50, col = "skyblue", border = "slateblue", main = "Histogram of p-values (SM Leptin Adlib vs Saline Adlib)")
ggsave("deseq2_output/histogram_SM_leptin_adlib_vs_saline_adlib.png")

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


