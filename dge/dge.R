# THIS SCRIPT WILL:
# Install necessary packages (DESeq2, tximport, ggplot2, EnhancedVolcano, and limma).
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

# Load packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(EnhancedVolcano)
library(limma)
library(scales) # needed for oob parameter
library(viridis)
library(reshape2)
library(grid)
library(gtable)
library(gridExtra)
library(ggdendro)
library(cluster)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Modify metadata to have a combined condition column
metadata$condition <- paste(metadata$Tissue, metadata$Injection, metadata$Feeding, sep = "_")

# Set fold-change and p-value thresholds
foldchange_threshold <- 4
pvalue_threshold <- 0.001
expression_threshold <- 2

# Import Salmon output
samples <- metadata$SampleID
files <- file.path("salmon_output", samples, "quant.sf")
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, metadata, ~ condition)
dds <- DESeq(dds)

# Define contrasts
contrasts <- list(
  c("condition", "EO_saline_fooddep", "EO_saline_adlib"),
  c("condition", "EO_saline_adlib", "SM_saline_fooddep"),
  c("condition", "EO_saline_fooddep", "SM_saline_fooddep"),
  c("condition", "EO_leptin_adlib", "EO_saline_adlib"),
  c("condition", "EO_leptin_adlib", "SM_leptin_adlib"),
  c("condition", "EO_leptin_fooddep", "EO_saline_fooddep"),
  c("condition", "EO_leptin_fooddep", "SM_leptin_fooddep"),
  c("condition", "EO_leptin_fooddep", "EO_leptin_adlib")
)

# Apply contrasts and get results
res_list <- lapply(contrasts, function(cntrst) {
  res <- results(dds, contrast=cntrst)
  return(res)
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


# Save results
lapply(names(res_list), function(x) {
  write.csv(as.data.frame(res_list[[x]]), file = paste0("deseq2_output/", x, ".csv"))
})

# Function to filter genes and save to file
save_filtered_genes <- function(res, contrast_name, pvalue_threshold, foldchange_threshold) {
  filtered <- subset(res, padj < pvalue_threshold & abs(log2FoldChange) > foldchange_threshold)
  output_directory <- "deseq2_output"
  output_filename <- paste0(output_directory, "/", contrast_name, "_filtered_genes.csv")
  
  # Extract the required columns
  filtered_data <- data.frame(
    Gene = row.names(filtered),
    log2FoldChange = filtered$log2FoldChange,
    pvalue = filtered$pvalue,
    padj = filtered$padj,
    baseMean = filtered$baseMean
  )
  
  # Print information for debugging
  print(paste("Saving data for contrast:", contrast_name))
  print(paste("Number of filtered genes:", nrow(filtered_data)))
  print(paste("Output file:", output_filename))
  
  # Save the data frame as a .csv file
  write.table(filtered_data, file = output_filename, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Return filtered data for further analysis
  return(filtered_data)
}

# Apply contrasts, filter genes, and save to files
filtered_res_list <- lapply(names(res_list), function(x) {
  filtered_res <- save_filtered_genes(res_list[[x]], x, pvalue_threshold, foldchange_threshold)
  return(filtered_res)
})

names(filtered_res_list) <- names(res_list)  # Naming the filtered results similar to original results

# Finding common genes
filtered_EO_saline_fooddep_vs_leptin_fooddep <- filtered_res_list[["EO_leptin_fooddep_vs_EO_saline_fooddep"]]
filtered_SM_leptin_fooddep_vs_EO_leptin_fooddep <- filtered_res_list[["EO_leptin_fooddep_vs_SM_leptin_fooddep"]]

common_genes <- intersect(filtered_EO_saline_fooddep_vs_leptin_fooddep$Gene, filtered_SM_leptin_fooddep_vs_EO_leptin_fooddep$Gene)
write.table(common_genes, file = "deseq2_output/common_filtered_genes-EO_leptin_fooddep_vs_EO_saline_fooddep-EO_leptin_fooddep_vs_SM_leptin_fooddep.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# PCA plot (note that this is not specific to a contrast, so is only done once)
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

# Histograms, volcano plots, MA plots
for (contrast_name in names(res_list)) {
  res <- res_list[[contrast_name]]
  filtered_res <- filtered_res_list[[contrast_name]]

  # Histogram
  hist <- ggplot(data.frame(x=res$pvalue), aes(x)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), col = "slateblue", fill = "skyblue") +
    geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
    labs(title = paste0("Histogram of p-values (", contrast_name, ")"), x = "p-value", y = "Count")
  ggsave(paste0("deseq2_output/histogram_", contrast_name, ".png"), plot = hist)
  ggsave(paste0("deseq2_output/histogram_", contrast_name, ".pdf"), plot = hist)

  # Volcano plot
  volcano <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = contrast_name,
    pCutoff = pvalue_threshold,
    FCcutoff = foldchange_threshold,
    pointSize = 1.8,
    labSize = 2.7)
  ggsave(paste0("deseq2_output/volcano_", contrast_name, ".png"), plot = volcano)
  ggsave(paste0("deseq2_output/volcano_", contrast_name, ".pdf"), plot = volcano)

  # MA plot
  res$baseMeanNew <- 1 / (10^log(res$baseMean + 1))
  enhanced_ma <- EnhancedVolcano(res,
    lab = rownames(res),
    title = paste0('MA plot: ', contrast_name),
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

  ggsave(paste0("deseq2_output/enhanced_ma_", contrast_name, ".png"), plot = enhanced_ma)
  ggsave(paste0("deseq2_output/enhanced_ma_", contrast_name, ".pdf"), plot = enhanced_ma)
}


