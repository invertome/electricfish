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
library(dplyr)
library(ggrepel)



# Set fold-change and p-value thresholds
lfc_threshold <- 4  
pvalue_threshold <- 0.001
cpm_threshold <- 10
log_cpm_threshold <- (log(cpm_threshold)+1)

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

# Function to filter genes and save to file
save_filtered_genes <- function(res, contrast_name, pvalue_threshold, foldchange_threshold) {
  
  res_data <- res@.Data[[1]] # get the data frame from the 'table' slot
  
  # Add debugging print statements
  print(paste("Type of logFC: ", class(res_data$logFC)))

  
  # Check for non-numeric entries
  non_numeric_entries <- which(!is.numeric(res_data$logFC))
  if(length(non_numeric_entries) > 0) {
    print(paste("Non-numeric entries found at indices: ", non_numeric_entries))
    print("Values:")
    print(res_data$logFC[non_numeric_entries])
  }
  
  # Check for NA values
  na_entries <- which(is.na(res_data$logFC))
  if(length(na_entries) > 0) {
    print(paste("NA entries found at indices: ", na_entries))
    print("Values:")
    print(res_data$logFC[na_entries])
    # Replace NA values with 0
    res_data$logFC[na_entries] <- 0
  }
  
  # Check for non-finite values (Inf and -Inf)
  non_finite_entries <- which(!is.finite(res_data$logFC))
  if(length(non_finite_entries) > 0) {
    print(paste("Non-finite entries found at indices: ", non_finite_entries))
    # Replace non-finite values with 0
    res_data$logFC[non_finite_entries] <- 0
  }
  
  # Use dplyr::filter for more explicit NA handling
  filtered <- dplyr::filter(res_data, PValue < pvalue_threshold & abs(logFC) > foldchange_threshold)
  
  output_directory <- "edger_output"
  output_filename <- paste0(output_directory, "/", contrast_name, "_filtered_genes.csv")
  
  # Extract the required columns
  filtered_data <- data.frame(
    Gene = row.names(filtered),
    logFC = filtered$logFC,
    PValue = filtered$PValue
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
  filtered_res <- save_filtered_genes(res_list[[x]], x, pvalue_threshold, lfc_threshold)
  return(filtered_res)
})

names(filtered_res_list) <- names(res_list)  # Naming the filtered results similar to original results

# Finding common genes
filtered_EO_saline_fooddep_vs_leptin_fooddep <- filtered_res_list[["EO_leptin_fooddep_vs_EO_saline_fooddep"]]
filtered_SM_leptin_fooddep_vs_EO_leptin_fooddep <- filtered_res_list[["EO_leptin_fooddep_vs_SM_leptin_fooddep"]]

common_genes <- intersect(filtered_EO_saline_fooddep_vs_leptin_fooddep$Gene, filtered_SM_leptin_fooddep_vs_EO_leptin_fooddep$Gene)
write.table(common_genes, file = "edger_output/common_filtered_genes-EO_leptin_fooddep_vs_EO_saline_fooddep-EO_leptin_fooddep_vs_SM_leptin_fooddep.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


# Finding all unique genes across all contrasts
all_genes <- Reduce(union, lapply(filtered_res_list, function(x) x$Gene))
write.table(all_genes, file = "edger_output/all_filtered_genes_all_contrasts.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Calculate the logCPM values
logCPM <- cpm(dds, log = TRUE)

# Perform PCA on the transposed logCPM values
pcaData <- prcomp(t(logCPM))

# Create a data frame from the PCA results
pca_df <- as.data.frame(pcaData$x)

# Include the experimental conditions in the data frame
pca_df$Tissue <- metadata$Tissue
pca_df$Injection <- metadata$Injection
pca_df$Feeding <- metadata$Feeding
pca_df$SampleName <- metadata$SampleID

# Calculate the percent variance explained by each principal component
percentVar <- round(100 * (pcaData$sdev^2 / sum(pcaData$sdev^2)))

# Plot the PCA with ggplot2
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Tissue, shape = Injection)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = SampleName), size = 1.8, max.overlaps = 10) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
print(pca_plot)

# Save the plot
if (!dir.exists("edger_output")) dir.create("edger_output")
ggsave("edger_output/pca_plot.png", plot = pca_plot)
ggsave("edger_output/pca_plot.pdf", plot = pca_plot)


# Plot histograms
for (i in 1:length(res_list)) {
  res_df <- as.data.frame(res_list[[i]])  # Convert TopTags object to a data frame
  p <- ggplot(res_df, aes(x = logFC)) + geom_histogram(bins = 50) + theme_bw() + xlab("log2 Fold Change")
  ggsave(filename = paste0("edger_output/", names(res_list[i]), "_histogram.png"), plot = p)
}

# Add this new block of code after the 'Create and save histograms for each comparison.' comment
# Create and save histograms of p-values for each comparison
for (i in 1:length(res_list)) {
  hist <- ggplot(data.frame(x=res_list[[i]]$table$PValue), aes(x)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), col = "slateblue", fill = "skyblue") +
    geom_vline(xintercept = c(0.05, 0.001), linetype = "dashed", color = "red") +
    labs(title = paste0("Histogram of p-values (", names(res_list[i]), ")"), x = "p-value", y = "Count")
  ggsave(paste0("edger_output/pvalue_histogram_", names(res_list[i]), ".png"), plot = hist)
  ggsave(paste0("edger_output/pvalue_histogram_", names(res_list[i]), ".pdf"), plot = hist)
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
