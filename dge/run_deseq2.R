# THIS SCRIPT WILL:
# Install necessary packages (DESeq2, tximport, ggplot2, pheatmap, EnhancedVolcano, and limma).
# Load packages.
# Read metadata from the "Sample_Metadata.csv" file.
# Set fold-change and p-value thresholds.
# Import Salmon output.
# Create DESeq2 object.
# Define contrasts for various comparisons.
# Apply contrasts and obtain results for each comparison.
# Save the results and summaries to CSV files in the "deseq2_output" directory.
# Create and save PCA plot.
# Create and save histograms for each comparison.
# Create and save volcano plots for each comparison.
# Create and save enhanced MA plots for each comparison.
# Create and save heatmaps.


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "tximport", "ggplot2", "pheatmap", "EnhancedVolcano", "limma"))
install.packages("scales")
install.packages("reshape2")
install.packages("viridis")

# Load packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(limma)
library(scales) # needed for oob parameter
library(viridis)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Set fold-change and p-value thresholds
foldchange_threshold <- 4
pvalue_threshold <- 0.000001
expression_threshold <- 2

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
  EO_leptin_adlib_vs_saline_adlib = c(0, 0, 1, 0, 0, 0, 1, 0),
  SM_leptin_fooddep_vs_saline_fooddep = c(0, 0, 1, -1, 0, 0, 1, 0),
  SM_leptin_adlib_vs_saline_adlib = c(0, 0, 1, 0, 0, 0, 0, 1),
  EO_leptin_fooddep_vs_EO_leptin_adlib = c(0, 0, 1, -1, 0, 0, 0, 0),
  SM_leptin_fooddep_vs_SM_leptin_adlib = c(0, 0, 1, 0, 0, 0, 0, -1),
  EO_saline_fooddep_vs_SM_saline_fooddep = c(0, 0, 0, 0, 1, -1, 0, 0),
  EO_saline_adlib_vs_SM_saline_adlib = c(0, 0, 0, 0, 1, 0, 0, -1),
  Interaction_EO_vs_SM_leptin_fooddep = c(0, 0, 1, -1, -1, 1, 0, 0),
  Interaction_EO_vs_SM_leptin_adlib = c(0, 0, 1, 0, -1, 0, 0, 1)
  
)


# Apply contrasts and obtain results
res_EO_leptin_fooddep_vs_saline_fooddep <- results(dds, contrast = contrasts[["EO_leptin_fooddep_vs_saline_fooddep"]])
res_EO_leptin_adlib_vs_saline_adlib <- results(dds, contrast = contrasts[["EO_leptin_adlib_vs_saline_adlib"]])
res_SM_leptin_fooddep_vs_saline_fooddep <- results(dds, contrast = contrasts[["SM_leptin_fooddep_vs_saline_fooddep"]])
res_SM_leptin_adlib_vs_saline_adlib <- results(dds, contrast = contrasts[["SM_leptin_adlib_vs_saline_adlib"]])
res_EO_leptin_fooddep_vs_EO_leptin_adlib <- results(dds, contrast = contrasts[["EO_leptin_fooddep_vs_EO_leptin_adlib"]])
res_SM_leptin_fooddep_vs_SM_leptin_adlib <- results(dds, contrast = contrasts[["SM_leptin_fooddep_vs_SM_leptin_adlib"]])
res_EO_saline_fooddep_vs_SM_saline_fooddep <- results(dds, contrast = contrasts[["EO_saline_fooddep_vs_SM_saline_fooddep"]])
res_EO_saline_adlib_vs_SM_saline_adlib <- results(dds, contrast = contrasts[["EO_saline_adlib_vs_SM_saline_adlib"]])
res_Interaction_EO_vs_SM_leptin_fooddep <- results(dds, contrast = contrasts[["Interaction_EO_vs_SM_leptin_fooddep"]])
res_Interaction_EO_vs_SM_leptin_adlib <- results(dds, contrast = contrasts[["Interaction_EO_vs_SM_leptin_adlib"]])

# Save results
dir.create("deseq2_output")
write.csv(res_EO_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/EO_leptin_fooddep_vs_saline_fooddep.csv")
write.csv(res_EO_leptin_adlib_vs_saline_adlib, file = "deseq2_output/EO_leptin_adlib_vs_saline_adlib.csv")
write.csv(res_SM_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/SM_leptin_fooddep_vs_saline_fooddep.csv")
write.csv(res_SM_leptin_adlib_vs_saline_adlib, file = "deseq2_output/SM_leptin_adlib_vs_saline_adlib.csv")

write.csv(res_EO_leptin_fooddep_vs_EO_leptin_adlib, file = "deseq2_output/EO_leptin_fooddep_vs_EO_leptin_adlib.csv")
write.csv(res_SM_leptin_fooddep_vs_SM_leptin_adlib, file = "deseq2_output/SM_leptin_fooddep_vs_SM_leptin_adlib.csv")
write.csv(res_EO_saline_fooddep_vs_SM_saline_fooddep, file = "deseq2_output/EO_saline_fooddep_vs_SM_saline_fooddep.csv")
write.csv(res_EO_saline_adlib_vs_SM_saline_adlib, file = "deseq2_output/EO_saline_adlib_vs_SM_saline_adlib.csv")
write.csv(res_Interaction_EO_vs_SM_leptin_fooddep, file = "deseq2_output/Interaction_EO_vs_SM_leptin_fooddep.csv")
write.csv(res_Interaction_EO_vs_SM_leptin_adlib, file = "deseq2_output/Interaction_EO_vs_SM_leptin_adlib.csv")


# Save summaries
summary_EO_leptin_fooddep_vs_saline_fooddep <- summary(res_EO_leptin_fooddep_vs_saline_fooddep)
summary_EO_leptin_adlib_vs_saline_adlib <- summary(res_EO_leptin_adlib_vs_saline_adlib)
summary_SM_leptin_fooddep_vs_saline_fooddep <- summary(res_SM_leptin_fooddep_vs_saline_fooddep)
summary_SM_leptin_adlib_vs_saline_adlib <- summary(res_SM_leptin_adlib_vs_saline_adlib)

summary_EO_leptin_fooddep_vs_EO_leptin_adlib <- summary(res_EO_leptin_fooddep_vs_EO_leptin_adlib)
summary_SM_leptin_fooddep_vs_SM_leptin_adlib <- summary(res_SM_leptin_fooddep_vs_SM_leptin_adlib)
summary_EO_saline_fooddep_vs_SM_saline_fooddep <- summary(res_EO_saline_fooddep_vs_SM_saline_fooddep)
summary_EO_saline_adlib_vs_SM_saline_adlib <- summary(res_EO_saline_adlib_vs_SM_saline_adlib)
summary_Interaction_EO_vs_SM_leptin_fooddep <- summary(res_Interaction_EO_vs_SM_leptin_fooddep)
summary_Interaction_EO_vs_SM_leptin_adlib <- summary(res_Interaction_EO_vs_SM_leptin_adlib)

write.table(summary_EO_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/summary_EO_leptin_fooddep_vs_saline_fooddep.txt")
write.table(summary_EO_leptin_adlib_vs_saline_adlib, file = "deseq2_output/summary_EO_leptin_adlib_vs_saline_adlib.txt")
write.table(summary_SM_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/summary_SM_leptin_fooddep_vs_saline_fooddep.txt")
write.table(summary_SM_leptin_adlib_vs_saline_adlib, file = "deseq2_output/summary_SM_leptin_adlib_vs_saline_adlib.txt")

write.table(summary_EO_leptin_fooddep_vs_EO_leptin_adlib, file = "deseq2_output/summary_EO_leptin_fooddep_vs_EO_leptin_adlib.txt")
write.table(summary_SM_leptin_fooddep_vs_SM_leptin_adlib, file = "deseq2_output/summary_SM_leptin_fooddep_vs_SM_leptin_adlib.txt")
write.table(summary_EO_saline_fooddep_vs_SM_saline_fooddep, file = "deseq2_output/summary_EO_saline_fooddep_vs_SM_saline_fooddep.txt")
write.table(summary_EO_saline_adlib_vs_SM_saline_adlib, file = "deseq2_output/summary_EO_saline_adlib_vs_SM_saline_adlib.txt")
write.table(summary_Interaction_EO_vs_SM_leptin_fooddep, file = "deseq2_output/summary_Interaction_EO_vs_SM_leptin_fooddep.txt")
write.table(summary_Interaction_EO_vs_SM_leptin_adlib, file = "deseq2_output/summary_Interaction_EO_vs_SM_leptin_adlib.txt")

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
ggsave("deseq2_output/pca_plot.pdf", plot = pca_plot)


# Histograms
hist_EO_leptin_fooddep_vs_saline_fooddep <- ggplot(data.frame(x=res_EO_leptin_fooddep_vs_saline_fooddep$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), col = "slateblue", fill = "skyblue") +
  geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
  labs(title = "Histogram of p-values (EO Leptin Fooddep vs Saline Fooddep)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_EO_leptin_fooddep_vs_saline_fooddep.png", plot = hist_EO_leptin_fooddep_vs_saline_fooddep)
ggsave("deseq2_output/histogram_EO_leptin_fooddep_vs_saline_fooddep.pdf", plot = hist_EO_leptin_fooddep_vs_saline_fooddep)


hist_EO_leptin_adlib_vs_saline_adlib <- ggplot(data.frame(x=res_EO_leptin_adlib_vs_saline_adlib$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), col = "slateblue", fill = "skyblue") +
  geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
  labs(title = "Histogram of p-values (EO Leptin Adlib vs Saline Adlib)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_EO_leptin_adlib_vs_saline_adlib.png", plot = hist_EO_leptin_adlib_vs_saline_adlib)
ggsave("deseq2_output/histogram_EO_leptin_adlib_vs_saline_adlib.pdf", plot = hist_EO_leptin_adlib_vs_saline_adlib)


hist_SM_leptin_fooddep_vs_saline_fooddep <- ggplot(data.frame(x=res_SM_leptin_fooddep_vs_saline_fooddep$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), col = "slateblue", fill = "skyblue") +
  geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
  labs(title = "Histogram of p-values (SM Leptin Fooddep vs Saline Fooddep)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_SM_leptin_fooddep_vs_saline_fooddep.png", plot = hist_SM_leptin_fooddep_vs_saline_fooddep)
ggsave("deseq2_output/histogram_SM_leptin_fooddep_vs_saline_fooddep.pdf", plot = hist_SM_leptin_fooddep_vs_saline_fooddep)


hist_SM_leptin_adlib_vs_saline_adlib <- ggplot(data.frame(x=res_SM_leptin_adlib_vs_saline_adlib$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), col = "slateblue", fill = "skyblue") +
  geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
  labs(title = "Histogram of p-values (SM Leptin Adlib vs Saline Adlib)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_SM_leptin_adlib_vs_saline_adlib.png", plot = hist_SM_leptin_adlib_vs_saline_adlib)
ggsave("deseq2_output/histogram_SM_leptin_adlib_vs_saline_adlib.pdf", plot = hist_SM_leptin_adlib_vs_saline_adlib)


# Volcano plots

# EO Leptin Fooddep vs Saline Fooddep
volcano_EO_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_EO_leptin_fooddep_vs_saline_fooddep,
  lab = rownames(res_EO_leptin_fooddep_vs_saline_fooddep),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'EO Leptin Fooddep vs Saline Fooddep',
  pCutoff = pvalue_threshold,
  FCcutoff = foldchange_threshold,
  pointSize = 1.8,
  labSize = 2.7)
ggsave("deseq2_output/volcano_EO_leptin_fooddep_vs_saline_fooddep.png", plot = volcano_EO_leptin_fooddep_vs_saline_fooddep)
ggsave("deseq2_output/volcano_EO_leptin_fooddep_vs_saline_fooddep.pdf", plot = volcano_EO_leptin_fooddep_vs_saline_fooddep)


# EO Leptin Adlib vs Saline Adlib
volcano_EO_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_EO_leptin_adlib_vs_saline_adlib,
  lab = rownames(res_EO_leptin_adlib_vs_saline_adlib),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'EO Leptin Adlib vs Saline Adlib',
  pCutoff = pvalue_threshold,
  FCcutoff = foldchange_threshold,
  pointSize = 1.8,
  labSize = 2.7)
ggsave("deseq2_output/volcano_EO_leptin_adlib_vs_saline_adlib.png", plot = volcano_EO_leptin_adlib_vs_saline_adlib)
ggsave("deseq2_output/volcano_EO_leptin_adlib_vs_saline_adlib.pdf", plot = volcano_EO_leptin_adlib_vs_saline_adlib)


# SM Leptin Fooddep vs Saline Fooddep
volcano_SM_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_SM_leptin_fooddep_vs_saline_fooddep,
  lab = rownames(res_SM_leptin_fooddep_vs_saline_fooddep),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'SM Leptin Fooddep vs Saline Fooddep',
  pCutoff = pvalue_threshold,
  FCcutoff = foldchange_threshold,
  pointSize = 1.8,
  labSize = 2.7)
ggsave("deseq2_output/volcano_SM_leptin_fooddep_vs_saline_fooddep.png", plot = volcano_SM_leptin_fooddep_vs_saline_fooddep)
ggsave("deseq2_output/volcano_SM_leptin_fooddep_vs_saline_fooddep.pdf", plot = volcano_SM_leptin_fooddep_vs_saline_fooddep)


# SM Leptin Adlib vs Saline Adlib
volcano_SM_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_SM_leptin_adlib_vs_saline_adlib,
  lab = rownames(res_SM_leptin_adlib_vs_saline_adlib),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'SM Leptin Adlib vs Saline Adlib',
  pCutoff = pvalue_threshold,
  FCcutoff = foldchange_threshold,
  pointSize = 1.8,
  labSize = 2.7)
ggsave("deseq2_output/volcano_SM_leptin_adlib_vs_saline_adlib.png", plot = volcano_SM_leptin_adlib_vs_saline_adlib)
ggsave("deseq2_output/volcano_SM_leptin_adlib_vs_saline_adlib.pdf", plot = volcano_SM_leptin_adlib_vs_saline_adlib)


# MA plots

# EO Leptin Fooddep vs Saline Fooddep
res_EO_leptin_fooddep_vs_saline_fooddep$baseMeanNew <- 1 / (10^log(res_EO_leptin_fooddep_vs_saline_fooddep$baseMean + 1))
enhanced_ma_EO_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_EO_leptin_fooddep_vs_saline_fooddep,
    lab = rownames(res_EO_leptin_fooddep_vs_saline_fooddep),
    title = 'MA plot: EO Leptin Fooddep vs Saline Fooddep',
    x = 'log2FoldChange',
    y = 'baseMeanNew',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'base mean + 1'),
    xlim = c(-40, 40),
    ylim = c(0, 12),
    pCutoff = pvalue_threshold,
    FCcutoff = foldchange_threshold,
    pointSize = 1.8,
    labSize = 2.7,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0) + coord_flip()

ggsave("deseq2_output/enhanced_ma_EO_leptin_fooddep_vs_saline_fooddep.png", plot = enhanced_ma_EO_leptin_fooddep_vs_saline_fooddep)
ggsave("deseq2_output/enhanced_ma_EO_leptin_fooddep_vs_saline_fooddep.pdf", plot = enhanced_ma_EO_leptin_fooddep_vs_saline_fooddep)


# EO Leptin Adlib vs Saline Adlib
res_EO_leptin_adlib_vs_saline_adlib$baseMeanNew <- 1 / (10^log(res_EO_leptin_adlib_vs_saline_adlib$baseMean + 1))
enhanced_ma_EO_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_EO_leptin_adlib_vs_saline_adlib,
    lab = rownames(res_EO_leptin_adlib_vs_saline_adlib),
    title = 'MA plot: EO Leptin Adlib vs Saline Adlib',
    x = 'log2FoldChange',
    y = 'baseMeanNew',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'base mean + 1'),
    xlim = c(-40, 40),
    ylim = c(0, 12),
    pCutoff = pvalue_threshold,
    FCcutoff = foldchange_threshold,
    pointSize = 1.8,
    labSize = 2.7,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0) + coord_flip()

ggsave("deseq2_output/enhanced_ma_EO_leptin_adlib_vs_saline_adlib.png", plot = enhanced_ma_EO_leptin_adlib_vs_saline_adlib)
ggsave("deseq2_output/enhanced_ma_EO_leptin_adlib_vs_saline_adlib.pdf", plot = enhanced_ma_EO_leptin_adlib_vs_saline_adlib)


# SM Leptin Fooddep vs Saline Fooddep
res_SM_leptin_fooddep_vs_saline_fooddep$baseMeanNew <- 1 / (10^log(res_SM_leptin_fooddep_vs_saline_fooddep$baseMean + 1))
enhanced_ma_SM_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_SM_leptin_fooddep_vs_saline_fooddep,
    lab = rownames(res_SM_leptin_fooddep_vs_saline_fooddep),
    title = 'MA plot: SM Leptin Fooddep vs Saline Fooddep',
    x = 'log2FoldChange',
    y = 'baseMeanNew',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'base mean + 1'),
    xlim = c(-40, 40),
    ylim = c(0, 12),
    pCutoff = pvalue_threshold,
    FCcutoff = foldchange_threshold,
    pointSize = 1.8,
    labSize = 2.7,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0) + coord_flip()

ggsave("deseq2_output/enhanced_ma_SM_leptin_fooddep_vs_saline_fooddep.png", plot = enhanced_ma_SM_leptin_fooddep_vs_saline_fooddep)
ggsave("deseq2_output/enhanced_ma_SM_leptin_fooddep_vs_saline_fooddep.pdf", plot = enhanced_ma_SM_leptin_fooddep_vs_saline_fooddep)


# SM Leptin Adlib vs Saline Adlib
res_SM_leptin_adlib_vs_saline_adlib$baseMeanNew <- 1 / (10^log(res_SM_leptin_adlib_vs_saline_adlib$baseMean + 1))
enhanced_ma_SM_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_SM_leptin_adlib_vs_saline_adlib,
    lab = rownames(res_SM_leptin_adlib_vs_saline_adlib),
    title = 'MA plot: SM Leptin Adlib vs Saline Adlib',
    x = 'log2FoldChange',
    y = 'baseMeanNew',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'base mean + 1'),
    xlim = c(-40, 40),
    ylim = c(0, 12),
    pCutoff = pvalue_threshold,
    FCcutoff = foldchange_threshold,
    pointSize = 1.8,
    labSize = 2.7,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0) + coord_flip()

ggsave("deseq2_output/enhanced_ma_SM_leptin_adlib_vs_saline_adlib.png", plot = enhanced_ma_SM_leptin_adlib_vs_saline_adlib)
ggsave("deseq2_output/enhanced_ma_SM_leptin_adlib_vs_saline_adlib.pdf", plot = enhanced_ma_SM_leptin_adlib_vs_saline_adlib)


# Heatmaps
# Load the reshape2 package
library(reshape2)

# Create a function to generate and save heatmaps using ggplot2
generate_ggplot_heatmap <- function(result_data, contrast_name, foldchange_threshold, pvalue_threshold, expression_threshold) {
  
  # Subset the results to meet the thresholds
  sig_genes <- subset(result_data, abs(log2FoldChange) >= log2(foldchange_threshold) & padj <= pvalue_threshold)
  
  # Extract the gene names
  sig_gene_names <- rownames(sig_genes)
  
  # Transform count data using the variance stabilizing transform
  ddsVST <- vst(dds)
  
  # Convert the DESeq transformed object to a data frame
  ddsVST_df <- as.data.frame(assay(ddsVST))
  ddsVST_df$Gene <- rownames(ddsVST_df)
  
  # Keep only the significantly differentiated genes
  ddsVST_filtered <- ddsVST_df[ddsVST_df$Gene %in% sig_gene_names,]
  
  # Convert the VST counts to long format for ggplot2
  ddsVST_long <- melt(ddsVST_filtered, id.vars=c("Gene"))
  
  # Create and save the heatmap using ggplot2
  heatmap_plot <- ggplot(ddsVST_long, aes(x=variable, y=Gene, fill=value)) +
    geom_raster() +
    scale_fill_viridis(trans = "sqrt") +
    theme(axis.text.x = element_text(angle = 65, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = paste("Heatmap of Differentially Expressed Genes:", contrast_name),
         x = "Samples",
         y = "Genes")
  
  # Save heatmap as PDF and PNG
  ggsave(file.path("deseq2_output", paste(contrast_name, "heatmap.pdf", sep = "_")), heatmap_plot, width = 8, height = 6, units = "in")
  ggsave(file.path("deseq2_output", paste(contrast_name, "heatmap.png", sep = "_")), heatmap_plot, width = 8, height = 6, units = "in")
}
# Transform count data using the variance stabilizing transform
ddsVST <- vst(dds)

# Convert the DESeq transformed object to a data frame
ddsVST_df <- as.data.frame(assay(ddsVST))
ddsVST_df$Gene <- rownames(ddsVST_df)

# Convert the VST counts to long format for ggplot2
ddsVST_long <- melt(ddsVST_df, id.vars=c("Gene"))

# Create and save the heatmap for all samples using ggplot2
all_samples_heatmap <- ggplot(ddsVST_long, aes(x=variable, y=Gene, fill=value)) +
  geom_raster() +
  scale_fill_viridis(trans = "sqrt") +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(title = "Heatmap of Differentially Expressed Genes: All Samples",
       x = "Samples",
       y = "Genes")

# Convert the significant genes back to a matrix for clustering
ddsVST_matrix <- dcast(ddsVST_long, Gene ~ variable)
rownames(ddsVST_matrix) <- ddsVST_matrix$Gene
ddsVST_matrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(ddsVST_matrix)
distanceSample <- dist(t(ddsVST_matrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construct a dendrogram for samples
install.packages("ggdendro")
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
ddsVST_long$variable <- factor(ddsVST_long$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(ddsVST_long, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combine the dendrogram and the heatmap
install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))

# Load in libraries necessary for modifying plots
library(gtable)
library(grid)

# Modify the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# Convert both grid-based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Draw the plot
grid.draw(finalGrob)

# Save the heatmap as PDF and PNG
pdf("deseq2_output/heatmap_with_dendrograms.pdf", height = 9, width = 6)
grid.draw(finalGrob)
dev.off()

png("deseq2_output/heatmap_with_dendrograms.png", height = 9, width = 6, units = "in", res = 300)
grid.draw(finalGrob)
dev.off()
