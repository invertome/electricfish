# Load required libraries
library(DESeq2)
library(readr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tximport)

# Read the sample metadata file
sample_metadata <- read.csv("Sample_Metadata.csv", stringsAsFactors = FALSE)

# Check if the sample metadata is read correctly
print(sample_metadata)

# Define the path to the Salmon output files for each sample
salmon_files <- file.path("salmon_output", sample_metadata$SampleID, "quant.sf")

# Check if the salmon_files variable is populated correctly
print(salmon_files)

# Load the first quant.sf file
first_quant_file <- read_tsv(salmon_files[1], col_types = cols())

# Extract transcript and gene IDs
transcript_ids <- first_quant_file$Name
gene_ids <- gsub("(t\\d+)$", "", transcript_ids)

# Create the tx2gene object
tx2gene <- data.frame(transcript = transcript_ids, gene = gene_ids)

# Import the transcript-level data using tximport
txi <- tximport(salmon_files, type = "salmon", txOut = TRUE, tx2gene = tx2gene)

# Remove intermediate objects and perform garbage collection
rm(first_quant_file, transcript_ids, gene_ids)
gc()

# Convert the summarized data to a DESeqDataSet object
dds <- DESeqDataSetFromTximport(txi, colData = sample_metadata, design = ~ Tissue + Injection * Feeding)

# Normalize and perform DESeq2 analysis
dds <- DESeq(dds)

# PCA plot
vsd <- vst(dds)
pca_data <- plotPCA(vsd, intgroup = c("Tissue", "Injection", "Feeding"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = Tissue, shape = Injection, linetype = Feeding)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA plot") +
  scale_color_discrete(name = "Tissue") +
  scale_shape_discrete(name = "Injection") +
  scale_linetype_discrete(name = "Feeding") +
  guides(shape = guide_legend(title = "Injection"), linetype = guide_legend(title = "Feeding")) +
  geom_text_repel(aes(label = rownames(pca_data)), force = 5)
ggsave("PCA_plot.png", width = 10, height = 8)

# Remove intermediate objects and perform garbage collection
rm(vsd, pca_data, percentVar)
gc()

# Define the comparisons of interest
contrasts <- list(
  "EO_leptin_fooddep_vs_EO_saline_fooddep" = c("Injection", "leptin", "saline"),
  "EO_leptin_adlib_vs_EO_saline_adlib" = c("Injection", "leptin", "saline"),
  "EO_leptin_fooddep_vs_EO_leptin_adlib" = c("Feeding", "fooddep", "adlib"),
  "SM_leptin_fooddep_vs_SM_saline_fooddep" = c("Injection", "leptin", "saline"),
  "SM_leptin_adlib_vs_SM_saline_adlib" = c("Injection", "leptin", "saline"),
  "SM_leptin_fooddep_vs_SM_leptin_adlib" = c("Feeding", "fooddep", "adlib"),
  "EO_leptin_fooddep_vs_SM_leptin_fooddep" = c("Tissue", "EO", "SM"),
  "EO_saline_fooddep_vs_SM_saline_fooddep" = c("Tissue", "EO", "SM"),
  "EO_leptin_adlib_vs_SM_leptin_adlib" = c("Tissue", "EO", "SM"),
  "EO_saline_adlib_vs_SM_saline_adlib" = c("Tissue", "EO", "SM")
)

# Perform the comparisons and store the results
results_list <- lapply(contrasts, function(contrast) {
  res <- results(dds, contrast = contrast)
  return(res)
})

names(results_list) <- names(contrasts)

# Combine all the results into a single data frame
results_list_df <- lapply(results_list, function(res) {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  return(res_df)
})

combined_results <- dplyr::bind_rows(results_list_df, .id = "comparison")
write.table(combined_results, file = "combined_DESeq2_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Write the results to a CSV file
write.csv(combined_results, "DESeq2_results.csv")


# Generate the heatmap
rld <- rlog(dds)
select_conditions <- colnames(rld)[colData(rld)$condition %in% c("EO_leptin_fooddep", "EO_saline_fooddep", "EO_leptin_adlib", "EO_saline_adlib",
                                                                 "SM_leptin_fooddep", "SM_saline_fooddep", "SM_leptin_adlib", "SM_saline_adlib")]

annotation <- as.data.frame(colData(rld)[, c("Tissue", "Injection", "Feeding")])
rownames(annotation) <- rownames(colData(rld))
pheatmap(assay(rld[, select_conditions]), annotation_col = annotation, show_colnames = FALSE)

# Generate a volcano plot for each comparison
volcano_plot <- function(res, comparison_name) {
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = padj < 0.05)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    ggtitle(paste0("Volcano plot: ", comparison_name)) +
    xlab("Log2 fold change") +
    ylab("-Log10 p-value") +
    scale_color_manual(values = c("gray", "red"), name = "Significant") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = subset(res, padj < 0.05 & abs(log2FoldChange) > 1),
      aes(label = rownames(res)),
      force = 5
    )
  ggsave(paste0("Volcano_plot_", comparison_name, ".png"), plot = p, width = 10, height = 8)
}


lapply(names(results_list), function(comparison_name) {
  volcano_plot(results_list[[comparison_name]], comparison_name)
})

# Remove intermediate objects and perform garbage collection
rm(combined_results, rld, select_conditions, annotation, results_list)
gc()

# Generate an MA plot for each comparison
ma_plot <- function(res, comparison_name) {
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  p <- ggplot(res, aes(x = baseMean, y = log2FoldChange, color = padj < 0.05)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    ggtitle(paste0("MA plot: ", comparison_name)) +
    xlab("Mean of normalized counts") +
    ylab("Log2 fold change") +
    scale_color_manual(values = c("gray", "red"), name = "Significant") +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = subset(res, padj < 0.05 & abs(log2FoldChange) > 1),
      aes(label = rownames(res)),
      force = 5
    )
  ggsave(paste0("MA_plot_", comparison_name, ".png"), plot = p, width = 10, height = 8)
}


lapply(names(contrasts), function(comparison_name) {
  ma_plot(results_list[[comparison_name]], comparison_name)
})

# Remove intermediate objects and perform garbage collection
rm(contrasts, dds)
gc()
