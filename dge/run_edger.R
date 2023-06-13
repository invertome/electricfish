# THIS SCRIPT WILL:
# Install necessary packages (edgeR, tximport, ggplot2, EnhancedVolcano, and limma).
# Load packages.
# Read metadata from the "Sample_Metadata.csv" file.
# Set fold-change and p-value thresholds.
# Import Salmon output.
# Create DGEList object.
# Define contrasts for various comparisons.
# Apply contrasts and obtain results for each comparison.
# Save the results and summaries to CSV files in the "edger_output" directory.
# Create and save PCA plot.
# Create and save histograms for each comparison.
# Create and save volcano plots for each comparison.
# Create and save enhanced MA plots for each comparison.
# Create and save a heatmap

# Load packages
library(edgeR)
library(tximport)
library(ggplot2)
library(EnhancedVolcano)
library(limma)
library(scales) # needed for oob parameter
library(viridis)
library(reshape2)
library(gplots)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Modify metadata to have a combined condition column
metadata$condition <- paste(metadata$Tissue, metadata$Injection, metadata$Feeding, sep = "_")

# Set fold-change and p-value thresholds
lfc_threshold <- log2(4)  # Convert fold-change threshold to log2 scale
pvalue_threshold <- 0.00001
cpm_threshold <- log2(2)  # Convert counts per million (CPM) threshold to log2 scale

# Import Salmon output
samples <- metadata$SampleID
files <- file.path("salmon_output", samples, "quant.sf")
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Create DGEList object
dds <- DGEList(counts = txi$counts, group = metadata$condition)
dds <- calcNormFactors(dds)

# Define contrasts
contrasts <- makeContrasts(
  "conditionEO_saline_fooddep - conditionEO_saline_adlib",
  "conditionEO_saline_adlib - conditionSM_saline_fooddep",
  "conditionEO_saline_fooddep - conditionSM_saline_fooddep",
  "conditionEO_leptin_adlib - conditionEO_saline_adlib",
  "conditionEO_leptin_adlib - conditionSM_leptin_adlib",
  "conditionEO_leptin_fooddep - conditionEO_saline_fooddep",
  "conditionEO_leptin_fooddep - conditionSM_leptin_fooddep",
  "conditionEO_leptin_fooddep - conditionEO_leptin_adlib",
  levels = design
)

# Apply contrasts and get results
dds <- estimateDisp(dds, design)
fit <- glmQLFit(dds, design)
res_list <- lapply(contrasts, function(cntrst) {
  res <- glmQLFTest(fit, contrast = cntrst)
  return(topTags(res, n = Inf))
})

names(res_list) <- c(
  "EO_saline_fooddep_vs_EO_saline_adlib",
  "EO_saline_adlib_vs_SM_saline_fooddep",
  "EO_saline_fooddep_vs_SM_saline_fooddep",
  "EO_leptin_adlib_vs_EO_saline_adlib",
  "EO_leptin_adlib_vs_SM_leptin_adlib",
  "EO_leptin_fooddep_vs_EO_saline_fooddep",
  "EO_leptin_fooddep_vs_SM_leptin_fooddep",
  "EO_leptin_fooddep_vs_EO_leptin_adlib"
)

# Save results and summaries to CSV files
if (!dir.exists("edger_output")) dir.create("edger_output")
for (i in 1:length(res_list)) {
  write.csv(res_list[[i]], file = paste0("edger_output/", names(res_list[i]), ".csv"))
}

# Plot PCA
logCPM <- cpm(dds, log = TRUE)
pca <- prcomp(t(logCPM))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- metadata$condition
p <- ggplot(pca_df, aes(PC1, PC2, color = condition)) + geom_point(size = 2) + theme_bw()
ggsave(filename = "edger_output/PCA_plot.png", plot = p)

# Plot histograms
for (i in 1:length(res_list)) {
  p <- ggplot(res_list[[i]], aes(x = logFC)) + geom_histogram(bins = 50) + theme_bw() + xlab("log2 Fold Change")
  ggsave(filename = paste0("edger_output/", names(res_list[i]), "_histogram.png"), plot = p)
}

# Plot volcano plots
for (i in 1:length(res_list)) {
  p <- EnhancedVolcano(res_list[[i]],
                       lab = rownames(res_list[[i]]),
                       x = 'logFC',
                       y = 'PValue',
                       title = names(res_list[i]),
                       pCutoff = pvalue_threshold,
                       FCcutoff = lfc_threshold)
  ggsave(filename = paste0("edger_output/", names(res_list[i]), "_volcano.png"), plot = p)
}

# Enhanced MA plots
for (i in 1:length(res_list)) {
  p <- plotSmear(res_list[[i]], de.tags = rownames(res_list[[i]][abs(res_list[[i]]$logFC) > lfc_threshold & res_list[[i]]$PValue < pvalue_threshold])) + geom_hline(yintercept = lfc_threshold, col = "red") + geom_hline(yintercept = -lfc_threshold, col = "red")
  ggsave(filename = paste0("edger_output/", names(res_list[i]), "_MA_plot.png"), plot = p)
}

# Create a heatmap
# Subset the dds data by removing genes that have an average log CPM less than the threshold.
keep <- rowMeans(cpm(dds) > cpm_threshold) >= 0.5
dds_subset <- dds[keep,]
logCPM <- cpm(dds_subset, log = TRUE)

# Do a hierarchical clustering of the genes and samples.
hc_gene <- hclust(dist(logCPM), method = "average")
hc_sample <- hclust(dist(t(logCPM)), method = "average")

# Draw the heatmap.
heatmap.2(
  as.matrix(logCPM),
  Rowv = as.dendrogram(hc_gene),
  Colv = as.dendrogram(hc_sample),
  scale = "none",
  col = viridis(256),
  trace = "none",
  dendrogram = "row",
  margin = c(12, 12)
)

# Save the heatmap
ggsave(filename = "edger_output/heatmap.png")

