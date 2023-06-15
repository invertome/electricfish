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


# Set fold-change and p-value thresholds
lfc_threshold <- 4  
pvalue_threshold <- 0.001
cpm_threshold <- (log(10)+1)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Import Salmon output
samples <- metadata$SampleID
files <- file.path("salmon_output", samples, "quant.sf")
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Modify metadata to have a combined condition column
metadata$condition <- paste(metadata$Tissue, metadata$Injection, metadata$Feeding, sep = "_")

# Create the design matrix
design <- model.matrix(~0 + condition, data = metadata)
colnames(design) <- levels(metadata$condition)

# Define the group and factor
group <- factor(metadata$condition)

# Generate the design matrix
design <- model.matrix(~0 + group)

# Make sure to assign column names to your design matrix that match your conditions
colnames(design) <- levels(group)

# Create DGEList object
dds <- DGEList(counts = txi$counts)
dds <- calcNormFactors(dds)

# Estimate dispersions
dds <- estimateDisp(dds, design)


# Define contrasts
contrasts <- makeContrasts(
  EO_saline_fooddep_vs_EO_saline_adlib = "EO_saline_fooddep - EO_saline_adlib",
  EO_saline_adlib_vs_SM_saline_fooddep = "EO_saline_adlib - SM_saline_fooddep",
  EO_saline_fooddep_vs_SM_saline_fooddep = "EO_saline_fooddep - SM_saline_fooddep",
  EO_leptin_adlib_vs_EO_saline_adlib = "EO_leptin_adlib - EO_saline_adlib",
  EO_leptin_adlib_vs_SM_leptin_adlib = "EO_leptin_adlib - SM_leptin_adlib",
  EO_leptin_fooddep_vs_EO_saline_fooddep = "EO_leptin_fooddep - EO_saline_fooddep",
  EO_leptin_fooddep_vs_SM_leptin_fooddep = "EO_leptin_fooddep - SM_leptin_fooddep",
  EO_leptin_fooddep_vs_EO_leptin_adlib = "EO_leptin_fooddep - EO_leptin_adlib",
  levels = colnames(design)
)


# Apply contrasts and get results
dds <- estimateDisp(dds, design)
fit <- glmQLFit(dds, design)
res_list <- lapply(colnames(contrasts), function(cntrst) {
  res <- glmQLFTest(fit, contrast = contrasts[,cntrst])
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
  res_df <- as.data.frame(res_list[[i]])  # Convert TopTags object to a data frame
  p <- ggplot(res_df, aes(x = logFC)) + geom_histogram(bins = 50) + theme_bw() + xlab("log2 Fold Change")
  ggsave(filename = paste0("edger_output/", names(res_list[i]), "_histogram.png"), plot = p)
}


# Plot volcano plots
for (i in 1:length(res_list)) {
  table_df <- res_list[[i]]@.Data[[1]]  # Access the data frame in the 'table' slot
  p <- EnhancedVolcano(table_df,
                       lab = rownames(table_df),
                       x = 'logFC',
                       y = 'PValue',
                       title = names(res_list[i]),
                       pCutoff = pvalue_threshold,
                       FCcutoff = lfc_threshold,
                       pointSize = 1.8,
                       labSize = 2.7)
  ggsave(filename = paste0("edger_output/", names(res_list[i]), "_volcano.png"), plot = p)
}



##################

# Enhanced MA plots
for (i in 1:length(res_list)) {
  contrast_name <- names(res_list)[i]
  res <- res_list[[i]]@.Data[[1]] # get the data frame from the 'table' slot
#  res$baseMeanNew <- 1 / (10^log(res$logCPM + 1)) # compute baseMeanNew
  enhanced_ma <- EnhancedVolcano(res,
    lab = rownames(res),
    title = paste0('MA plot: ', contrast_name),
    x = 'logFC',
    y = 'logCPM',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[e]~ 'mean expression'),
    xlim = c(-20, 20),
    ylim = c(0, 5),
    pCutoff = 0.01,
    FCcutoff = lfc_threshold,
    pointSize = 1.8,
    labSize = 2.7,
    legendLabels = c('NS', expression(Log[2]~FC),
      'Mean expression', expression(Mean-expression~and~log[2]~FC)),
    legendPosition = 'bottom',
    legendLabSize = 16,
    legendIconSize = 4.0) + coord_flip()
  
  ggsave(paste0("edger_output/enhanced_ma_", contrast_name, ".png"), plot = enhanced_ma)
  ggsave(paste0("edger_output/enhanced_ma_", contrast_name, ".pdf"), plot = enhanced_ma)
}



# Create a heatmap
# Subset the dds data by removing genes that have an average log CPM less than the threshold.
keep <- rowMeans(cpm(dds) > cpm_threshold) >= 0.5
dds_subset <- dds[keep,]
logCPM <- cpm(dds_subset, log = TRUE)

hc_gene <- hclust(dist(logCPM), method = "average")
hc_sample <- hclust(dist(t(logCPM)), method = "average")

# PNG Output
png(filename = "edger_output/heatmap.png")
heatmap.2(
  as.matrix(logCPM),
  Rowv = as.dendrogram(hc_gene),
  labRow = FALSE,
  Colv = as.dendrogram(hc_sample),
  scale = "none",
  col = viridis(256),
  trace = "none",
  dendrogram = "column",
  margin = c(12, 12)
)
dev.off()

# PDF Output
pdf(file = "edger_output/heatmap.pdf")
heatmap.2(
  as.matrix(logCPM),
  Rowv = as.dendrogram(hc_gene),
  labRow = FALSE,
  Colv = as.dendrogram(hc_sample),
  scale = "none",
  col = viridis(256),
  trace = "none",
  dendrogram = "column",
  margin = c(12, 12)
)
dev.off()
