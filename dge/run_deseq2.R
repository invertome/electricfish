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


# Load packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(limma)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Set fold-change and p-value thresholds
foldchange_threshold <- 2
pvalue_threshold <- 0.000001

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
  EO_leptin_adlib_vs_saline_adlib = c(0, 0, 1, 0, 0, 0, 1, 0), # <-- Updated this line
  SM_leptin_fooddep_vs_saline_fooddep = c(0, 0, 1, -1, 0, 0, 1, 0),
  SM_leptin_adlib_vs_saline_adlib = c(0, 0, 1, 0, 0, 0, 0, 1)  # <-- Updated this line
  EO_leptin_fooddep_vs_EO_leptin_adlib = c(0, 0, 1, -1, 0, 0, 0, 0)
  SM_leptin_fooddep_vs_SM_leptin_adlib = c(0, 0, 1, 0, 0, 0, 0, -1)
  EO_saline_fooddep_vs_SM_saline_fooddep = c(0, 0, 0, 0, 1, -1, 0, 0)
  EO_saline_adlib_vs_SM_saline_adlib = c(0, 0, 0, 0, 1, 0, 0, -1)
  Interaction_EO_vs_SM_leptin_fooddep = c(0, 0, 1, -1, -1, 1, 0, 0)
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
write.csv(res_SM_leptin_fooddep_vs_SM_leptin_adlib, file = "deseq2_output/SM_leptin_fooddep_vs_SM_leptin_adlib")
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

write.csv(summary_EO_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/summary_EO_leptin_fooddep_vs_saline_fooddep.csv")
write.csv(summary_EO_leptin_adlib_vs_saline_adlib, file = "deseq2_output/summary_EO_leptin_adlib_vs_saline_adlib.csv")
write.csv(summary_SM_leptin_fooddep_vs_saline_fooddep, file = "deseq2_output/summary_SM_leptin_fooddep_vs_saline_fooddep.csv")
write.csv(summary_SM_leptin_adlib_vs_saline_adlib, file = "deseq2_output/summary_SM_leptin_adlib_vs_saline_adlib.csv")

write.csv(summary_EO_leptin_fooddep_vs_EO_leptin_adlib, file = "deseq2_output/summary_EO_leptin_fooddep_vs_EO_leptin_adlib.csv")
write.csv(summary_SM_leptin_fooddep_vs_SM_leptin_adlib, file = "deseq2_output/summary_SM_leptin_fooddep_vs_SM_leptin_adlib.csv")
write.csv(summary_EO_saline_fooddep_vs_SM_saline_fooddep, file = "deseq2_output/summary_EO_saline_fooddep_vs_SM_saline_fooddep.csv")
write.csv(summary_EO_saline_adlib_vs_SM_saline_adlib, file = "deseq2_output/summary_EO_saline_adlib_vs_SM_saline_adlib")
write.csv(summary_Interaction_EO_vs_SM_leptin_fooddep, file = "deseq2_output/summary_Interaction_EO_vs_SM_leptin_fooddep.csv")
write.csv(summary_Interaction_EO_vs_SM_leptin_adlib, file = "deseq2_output/summary_Interaction_EO_vs_SM_leptin_adlib.csv")

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
  geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
  labs(title = "Histogram of p-values (EO Leptin Fooddep vs Saline Fooddep)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_EO_leptin_fooddep_vs_saline_fooddep.png", plot = hist_EO_leptin_fooddep_vs_saline_fooddep)

hist_EO_leptin_adlib_vs_saline_adlib <- ggplot(data.frame(x=res_EO_leptin_adlib_vs_saline_adlib$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), col = "slateblue", fill = "skyblue") +
  geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
  labs(title = "Histogram of p-values (EO Leptin Adlib vs Saline Adlib)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_EO_leptin_adlib_vs_saline_adlib.png", plot = hist_EO_leptin_adlib_vs_saline_adlib)

hist_SM_leptin_fooddep_vs_saline_fooddep <- ggplot(data.frame(x=res_SM_leptin_fooddep_vs_saline_fooddep$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), col = "slateblue", fill = "skyblue") +
  geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
  labs(title = "Histogram of p-values (SM Leptin Fooddep vs Saline Fooddep)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_SM_leptin_fooddep_vs_saline_fooddep.png", plot = hist_SM_leptin_fooddep_vs_saline_fooddep)

hist_SM_leptin_adlib_vs_saline_adlib <- ggplot(data.frame(x=res_SM_leptin_adlib_vs_saline_adlib$pvalue), aes(x)) +
  geom_histogram(breaks = seq(0, 1, by = 0.02), col = "slateblue", fill = "skyblue") +
  geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
  labs(title = "Histogram of p-values (SM Leptin Adlib vs Saline Adlib)", x = "p-value", y = "Count")
ggsave("deseq2_output/histogram_SM_leptin_adlib_vs_saline_adlib.png", plot = hist_SM_leptin_adlib_vs_saline_adlib)

# Volcano plots

# EO Leptin Fooddep vs Saline Fooddep
volcano_EO_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_EO_leptin_fooddep_vs_saline_fooddep,
  lab = rownames(res_EO_leptin_fooddep_vs_saline_fooddep),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'EO Leptin Fooddep vs Saline Fooddep',
  pCutoff = 0.05,
  FCcutoff = 1)
ggsave("deseq2_output/volcano_EO_leptin_fooddep_vs_saline_fooddep.png", plot = volcano_EO_leptin_fooddep_vs_saline_fooddep)

# EO Leptin Adlib vs Saline Adlib
volcano_EO_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_EO_leptin_adlib_vs_saline_adlib,
  lab = rownames(res_EO_leptin_adlib_vs_saline_adlib),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'EO Leptin Adlib vs Saline Adlib',
  pCutoff = 0.05,
  FCcutoff = 1)
ggsave("deseq2_output/volcano_EO_leptin_adlib_vs_saline_adlib.png", plot = volcano_EO_leptin_adlib_vs_saline_adlib)

# SM Leptin Fooddep vs Saline Fooddep
volcano_SM_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_SM_leptin_fooddep_vs_saline_fooddep,
  lab = rownames(res_SM_leptin_fooddep_vs_saline_fooddep),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'SM Leptin Fooddep vs Saline Fooddep',
  pCutoff = 0.05,
  FCcutoff = 1)
ggsave("deseq2_output/volcano_SM_leptin_fooddep_vs_saline_fooddep.png", plot = volcano_SM_leptin_fooddep_vs_saline_fooddep)

# SM Leptin Adlib vs Saline Adlib
volcano_SM_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_SM_leptin_adlib_vs_saline_adlib,
  lab = rownames(res_SM_leptin_adlib_vs_saline_adlib),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'SM Leptin Adlib vs Saline Adlib',
  pCutoff = 0.05,
  FCcutoff = 1)
ggsave("deseq2_output/volcano_SM_leptin_adlib_vs_saline_adlib.png", plot = volcano_SM_leptin_adlib_vs_saline_adlib)

# MA plots

# EO Leptin Fooddep vs Saline Fooddep
enhanced_ma_EO_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_EO_leptin_fooddep_vs_saline_fooddep,
    lab = rownames(res_EO_leptin_fooddep_vs_saline_fooddep),
    title = 'MA plot: EO Leptin Fooddep vs Saline Fooddep',
    x = 'log2FoldChange',
    y = 'baseMean',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'base mean + 1'),
    ylim = c(0, 12),
    pCutoff = 1 / (10 ^ 7),
    FCcutoff = 2.0,
    pointSize = 3.5,
    labSize = 4.0,
    boxedLabels = TRUE,
    colAlpha = 1,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75) + coord_flip()

ggsave("deseq2_output/enhanced_ma_EO_leptin_fooddep_vs_saline_fooddep.png", plot = enhanced_ma_EO_leptin_fooddep_vs_saline_fooddep)

# EO Leptin Adlib vs Saline Adlib
enhanced_ma_EO_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_EO_leptin_adlib_vs_saline_adlib,
    lab = rownames(res_EO_leptin_adlib_vs_saline_adlib),
    title = 'MA plot: EO Leptin Adlib vs Saline Adlib',
    x = 'log2FoldChange',
    y = 'baseMean',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'base mean + 1'),
    ylim = c(0, 12),
    pCutoff = 1 / (10 ^ 7),
    FCcutoff = 2.0,
    pointSize = 3.5,
    labSize = 4.0,
    boxedLabels = TRUE,
    colAlpha = 1,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75) + coord_flip()

ggsave("deseq2_output/enhanced_ma_EO_leptin_adlib_vs_saline_adlib.png", plot = enhanced_ma_EO_leptin_adlib_vs_saline_adlib)

# SM Leptin Fooddep vs Saline Fooddep
enhanced_ma_SM_leptin_fooddep_vs_saline_fooddep <- EnhancedVolcano(res_SM_leptin_fooddep_vs_saline_fooddep,
    lab = rownames(res_SM_leptin_fooddep_vs_saline_fooddep),
    title = 'MA plot: SM Leptin Fooddep vs Saline Fooddep',
    x = 'log2FoldChange',
    y = 'baseMean',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'base mean + 1'),
    ylim = c(0, 12),
    pCutoff = 1 / (10 ^ 7),
    FCcutoff = 2.0,
    pointSize = 3.5,
    labSize = 4.0,
    boxedLabels = TRUE,
    colAlpha = 1,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75) + coord_flip()

ggsave("deseq2_output/enhanced_ma_SM_leptin_fooddep_vs_saline_fooddep.png", plot = enhanced_ma_SM_leptin_fooddep_vs_saline_fooddep)

# SM Leptin Adlib vs Saline Adlib
enhanced_ma_SM_leptin_adlib_vs_saline_adlib <- EnhancedVolcano(res_SM_leptin_adlib_vs_saline_adlib,
    lab = rownames(res_SM_leptin_adlib_vs_saline_adlib),
    title = 'MA plot: SM Leptin Adlib vs Saline Adlib',
    x = 'log2FoldChange',
    y = 'baseMean',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'base mean + 1'),
    ylim = c(0, 12),
    pCutoff = 1 / (10 ^ 7),
    FCcutoff = 2.0,
    pointSize = 3.5,
    labSize = 4.0,
    boxedLabels = TRUE,
    colAlpha = 1,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75) + coord_flip()

ggsave("deseq2_output/enhanced_ma_SM_leptin_adlib_vs_saline_adlib.png", plot = enhanced_ma_SM_leptin_adlib_vs_saline_adlib)


# Heatmaps
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
log_norm_counts <- log2(norm_counts + 1)

# Identify rows with missing or infinite values
problematic_rows <- apply(log_norm_counts, 1, function(x) any(is.na(x) | is.infinite(x)))

# Remove problematic rows from the data
filtered_log_norm_counts <- log_norm_counts[!problematic_rows, ]

# Remove rows with near-zero variance
filtered_log_norm_counts <- filtered_log_norm_counts[apply(filtered_log_norm_counts, 1, var) > 1e-10,]

# Remove columns with all-zero values
all_zero_columns <- apply(filtered_log_norm_counts, 2, function(x) all(x == 0))
filtered_log_norm_counts <- filtered_log_norm_counts[, !all_zero_columns]

# Update top variable genes matrix
top_var_genes <- head(order(rowVars(filtered_log_norm_counts), decreasing = TRUE), 100)
mat_top_var_genes <- filtered_log_norm_counts[top_var_genes, ]

# Heatmap for EO Leptin Fooddep vs Saline Fooddep
EO_leptin_fooddep_vs_saline_fooddep <- mat_top_var_genes[, colData(dds)$group %in% c("EO_leptin_fooddep", "EO_saline_fooddep")]
pheatmap(EO_leptin_fooddep_vs_saline_fooddep, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", annotation_col = colData(dds)[, "group", drop = FALSE])

# Heatmap for EO Leptin Adlib vs Saline Adlib
EO_leptin_adlib_vs_saline_adlib <- mat_top_var_genes[, colData(dds)$group %in% c("EO_leptin_adlib", "EO_saline_adlib")]
pheatmap(EO_leptin_adlib_vs_saline_adlib, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", annotation_col = colData(dds)[, "group", drop = FALSE])

# Heatmap for SM Leptin Fooddep vs Saline Fooddep
SM_leptin_fooddep_vs_saline_fooddep <- mat_top_var_genes[, colData(dds)$group %in% c("SM_leptin_fooddep", "SM_saline_fooddep")]
pheatmap(SM_leptin_fooddep_vs_saline_fooddep, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", annotation_col = colData(dds)[, "group", drop = FALSE])

# Heatmap for SM Leptin Adlib vs Saline Adlib
SM_leptin_adlib_vs_saline_adlib <- mat_top_var_genes[, colData(dds)$group %in% c("SM_leptin_adlib", "SM_saline_adlib")]
pheatmap(SM_leptin_adlib_vs_saline_adlib, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", annotation_col = colData(dds)[, "group", drop = FALSE])

# General heatmap with hierarchical clustering of all the samples together
pheatmap(mat_top_var_genes, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", annotation_col = colData(dds)[, "group", drop = FALSE])



