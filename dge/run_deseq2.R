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
# Create and save a heatmap

# Load packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(EnhancedVolcano)
library(limma)
library(scales) # needed for oob parameter
library(viridis)
library(reshape2)
library(gplots)
library(ggrepel)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Modify metadata to have a combined condition column
metadata$condition <- paste(metadata$Tissue, metadata$Injection, metadata$Feeding, sep = "_")

# Set fold-change and p-value thresholds
foldchange_threshold <- 4
pvalue_threshold <- 0.001
expression_threshold <- 10

# Import Salmon output
samples <- metadata$SampleID
files <- file.path("salmon_output", samples, "quant.sf")
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, metadata, ~ condition)
dds <- DESeq(dds)

# Define contrasts
contrasts <- list(
  c("condition", "EO_saline_adlib", "SM_saline_adlib"),
  c("condition", "EO_saline_fooddep", "SM_saline_fooddep"),
  c("condition", "EO_leptin_adlib", "SM_leptin_adlib"),
  c("condition", "EO_leptin_fooddep", "SM_leptin_fooddep"),
  c("condition", "EO_saline_fooddep", "EO_saline_adlib"),
  c("condition", "EO_leptin_fooddep", "EO_leptin_adlib"),
  c("condition", "SM_leptin_adlib", "SM_saline_adlib"),
  c("condition", "SM_leptin_fooddep", "SM_saline_fooddep")
)

# Apply contrasts and get results
res_list <- lapply(contrasts, function(cntrst) {
  res <- results(dds, contrast=cntrst)
  return(res)
})

names(res_list) <- c(
  "EO_saline_adlib_vs_SM_saline_adlib",
  "EO_saline_fooddep_vs_SM_saline_fooddep",
  "EO_leptin_adlib_vs_SM_leptin_adlib",
  "EO_leptin_fooddep_vs_SM_leptin_fooddep",
  "EO_saline_fooddep_vs_EO_saline_adlib",
  "EO_leptin_fooddep_vs_EO_leptin_adlib",
  "SM_leptin_adlib_vs_SM_saline_adlib",
  "SM_leptin_fooddep_vs_SM_saline_fooddep"
)





# Save results
lapply(names(res_list), function(x) {
  write.csv(as.data.frame(res_list[[x]]), file = paste0("deseq2_output/", x, ".csv"))
})

# Function to filter genes and save to file
save_filtered_genes <- function(res, contrast_name, pvalue_threshold, foldchange_threshold, expression_threshold) {
  filtered <- subset(res, padj < pvalue_threshold & abs(log2FoldChange) > foldchange_threshold & baseMean > expression_threshold)
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
  filtered_res <- save_filtered_genes(res_list[[x]], x, pvalue_threshold, foldchange_threshold, expression_threshold)
  return(filtered_res)
})


names(filtered_res_list) <- names(res_list)  # Naming the filtered results similar to original results

# Finding common genes
filtered_EO_saline_fooddep_vs_leptin_fooddep <- filtered_res_list[["EO_leptin_fooddep_vs_EO_saline_fooddep"]]
filtered_SM_leptin_fooddep_vs_EO_leptin_fooddep <- filtered_res_list[["EO_leptin_fooddep_vs_SM_leptin_fooddep"]]

common_genes <- intersect(filtered_EO_saline_fooddep_vs_leptin_fooddep$Gene, filtered_SM_leptin_fooddep_vs_EO_leptin_fooddep$Gene)
write.table(common_genes, file = "deseq2_output/common_filtered_genes-EO_leptin_fooddep_vs_EO_saline_fooddep-EO_leptin_fooddep_vs_SM_leptin_fooddep.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# PCA plot
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("Tissue", "Injection", "Feeding"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = Tissue, shape = Injection, fill = Feeding)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("adlib" = "black", "fooddep" = NA)) +
  geom_text_repel(aes(label = rownames(pcaData)), size = 1.8, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = c("black", "black"), fill = c("black", "white"))))  # Manually define legend aesthetics

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
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +
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

  # Heatmap
  # Heatmap
  tryCatch({
    # Extract condition names from the contrast
    cond1 <- strsplit(contrast_name, "_vs_")[[1]][1]
    cond2 <- strsplit(contrast_name, "_vs_")[[1]][2]
  
    # Subset the samples for this contrast
    samples_for_contrast <- colnames(assay(vsd)) %in% c(metadata[metadata$condition == cond1,]$SampleID, metadata[metadata$condition == cond2,]$SampleID)
  
    # Select genes with finite values only
    select_genes <- rowSums(is.finite(assay(vsd)[select_genes, samples_for_contrast])) == length(samples_for_contrast)
  
    # Check that there are genes to plot
    if (sum(select_genes) > 0) {
      top_variance_genes <- head(order(rowVars(assay(vsd)[select_genes, samples_for_contrast]), decreasing = TRUE), 100)
    
      # Check that there are top variance genes to plot
      if(length(top_variance_genes) > 0) {
        heatmap_data <- assay(vsd)[top_variance_genes, samples_for_contrast]
        rownames(heatmap_data) <- make.names(rownames(heatmap_data), unique = TRUE)
      
        heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
      
        heatmap <- pheatmap(heatmap_data,
                            color = heatmap_colors,
                            show_rownames = FALSE,
                            scale = "row",
                            annotation_col = as.data.frame(colData(vsd)[samples_for_contrast, "condition"]),
                            main = paste0("Heatmap for ", contrast_name))
      
        # Check if directory exists
        if(!dir.exists("deseq2_output")) {
          dir.create("deseq2_output")
        }
      
        ggsave(paste0("deseq2_output/heatmap_", contrast_name, ".png"), plot = heatmap)
        ggsave(paste0("deseq2_output/heatmap_", contrast_name, ".pdf"), plot = heatmap)
      } else {
        print(paste("No top variance genes found for contrast:", contrast_name))
      }
    } else {
      print(paste("No genes selected for contrast:", contrast_name))
    }
  }, error=function(e){print(paste("Error in creating heatmap for contrast:", contrast_name, "Error:", e))})



} 
  


# Apply contrasts, filter genes, and save to files
filtered_res_list <- lapply(names(res_list), function(x) {
  filtered_res <- save_filtered_genes(res_list[[x]], x, pvalue_threshold, foldchange_threshold)
  return(filtered_res)
})

names(filtered_res_list) <- names(res_list)  # Naming the filtered results similar to original results

# Finding all unique genes across all contrasts
all_genes <- Reduce(union, lapply(filtered_res_list, function(x) x$Gene))
write.table(all_genes, file = "deseq2_output/all_filtered_genes_all_contrasts.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Extract normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# Filter normalized counts for the all genes across all contrasts
filtered_normalized_counts <- normalized_counts[all_genes,]

# Z-scores and log-transformed counts for the filtered genes
filtered_log_transformed_counts <- log2(filtered_normalized_counts + 1)

# Log10-transformed counts for the filtered genes
filtered_log10_transformed_counts <- log10(filtered_normalized_counts + 1)

# Variance Stabilizing Transformed counts for the filtered genes
vst_transformed_counts <- assay(vst(dds, blind=FALSE))
filtered_vst_transformed_counts <- vst_transformed_counts[all_genes,]

# Heatmap function to avoid repetition
plot_heatmap <- function(counts_matrix, filename, title){
  # Save the heatmap to a .png file
  png(paste0("deseq2_output/", filename, ".png"))
  heatmap.2(as.matrix(counts_matrix),
            trace = "none",
            margins = c(5,10),
            col = viridis(100),
            dendrogram = 'both',
            Rowv = TRUE,
            Colv = TRUE,
            labCol = colnames(counts_matrix),
            labRow = FALSE,
            scale = "none",
            key = TRUE,
            keysize = 1.2,
            cexRow = 0.6,
            cexCol = 0.6,
            main = title)
  dev.off()
  
  # Save the heatmap to a .pdf file
  pdf(paste0("deseq2_output/", filename, ".pdf"))
  heatmap.2(as.matrix(counts_matrix),
            trace = "none",
            margins = c(5,10),
            col = viridis(100),
            dendrogram = 'both',
            Rowv = TRUE,
            Colv = TRUE,
            labCol = colnames(counts_matrix),
            labRow = FALSE,
            scale = "none",
            key = TRUE,
            keysize = 1.2,
            cexRow = 0.6,
            cexCol = 0.6,
            main = title)
  dev.off()
}

# Plot heatmaps for log2-transformed, log10-transformed, and vst-transformed counts
plot_heatmap(filtered_log_transformed_counts, "filtered_heatmap_log2_transformed_all_contrasts", "Gene Expression Heatmap (log2)")
plot_heatmap(filtered_log10_transformed_counts, "filtered_heatmap_log10_transformed_all_contrasts", "Gene Expression Heatmap (log10)")
plot_heatmap(filtered_vst_transformed_counts, "filtered_heatmap_vst_transformed_all_contrasts", "Gene Expression Heatmap (VST)")
                                  
# We will use the top 100 genes with the highest variance across samples for the heatmap.
# You can adjust this number according to your preference.
top_variance_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100)

heatmap_data <- assay(vsd)[top_variance_genes, ]
rownames(heatmap_data) <- make.names(rownames(heatmap_data), unique = TRUE)

heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

heatmap <- pheatmap(heatmap_data,
                    color = heatmap_colors,
                    show_rownames = FALSE,
                    scale = "row",
                    annotation_col = as.data.frame(colData(vsd)[, "condition"]),
                    main = paste0("VSD normalized Heatmap", contrast_name))
ggsave(paste0("deseq2_output/heatmap_", "top100", ".png"), plot = heatmap)
ggsave(paste0("deseq2_output/heatmap_", "top100", ".pdf"), plot = heatmap)



