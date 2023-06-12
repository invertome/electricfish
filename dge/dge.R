# Load packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(EnhancedVolcano)
library(limma)
library(scales) # needed for oob parameter
library(viridis)

# Read metadata
metadata <- read.csv("Sample_Metadata.csv", header = TRUE)

# Modify metadata to have a combined condition column
metadata$condition <- paste(metadata$Tissue, metadata$Injection, metadata$Feeding, sep = "_")

# Set fold-change and p-value thresholds
foldchange_threshold <- 7
pvalue_threshold <- 0.000001
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
  c("condition", "EO_saline_adlib", "EO_saline_fooddep"),
  c("condition", "SM_saline_fooddep", "EO_saline_adlib"),
  c("condition", "SM_saline_fooddep", "EO_saline_fooddep"),
  c("condition", "EO_saline_adlib", "EO_leptin_adlib"),
  c("condition", "SM_leptin_adlib", "EO_leptin_adlib"),
  c("condition", "EO_saline_fooddep", "EO_leptin_fooddep"),
  c("condition", "SM_leptin_fooddep", "EO_leptin_fooddep"),
  c("condition", "EO_leptin_adlib", "EO_leptin_fooddep")
)

# Apply contrasts and get results
res_list <- lapply(contrasts, function(cntrst) {
  res <- results(dds, contrast=cntrst)
  return(res)
})

names(res_list) <- c(
  "EO_saline_adlib_vs_EO_saline_fooddep",
  "SM_saline_fooddep_vs_EO_saline_adlib",
  "SM_saline_fooddep_vs_SM_saline_fooddep",
  "EO_saline_adlib_vs_EO_leptin_adlib",
  "SM_leptin_adlib_vs_EO_leptin_adlib",
  "EO_saline_fooddep_vs_EO_leptin_fooddep",
  "SM_leptin_fooddep_vs_EO_leptin_fooddep",
  "EO_leptin_adlib_vs_EO_leptin_fooddep"
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
filtered_EO_leptin_fooddep_vs_saline_fooddep <- filtered_res_list[["EO_saline_fooddep_vs_EO_leptin_fooddep"]]
filtered_EO_leptin_fooddep_vs_SM_leptin_fooddep <- filtered_res_list[["SM_leptin_fooddep_vs_EO_leptin_fooddep"]]

common_genes <- intersect(filtered_EO_leptin_fooddep_vs_saline_fooddep$Gene, filtered_EO_leptin_fooddep_vs_SM_leptin_fooddep$Gene)
write.table(common_genes, file = "deseq2_output/common_filtered_genes_EO_saline_fooddep_vs_EO_leptin_fooddep-SM_leptin_fooddep_vs_EO_leptin_fooddep.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
