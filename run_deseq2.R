# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tximport)

# Read the count matrix and sample metadata files
count_matrix <- read.csv("Count_Matrix.csv", row.names = 1)
sample_metadata <- read.csv("Sample_Metadata.csv")

# Create the tx2gene data frame
transcript_ids <- rownames(count_matrix)
gene_ids <- gsub("t\\d+$", "", transcript_ids)
tx2gene <- data.frame(transcript = transcript_ids, gene = gene_ids)

# Import the non-integer count matrix as a list
counts_list <- list(counts = as.matrix(count_matrix))

# Custom function to create a Tximport object
custom_tximport <- function(counts_list, tx2gene) {
  list(abundance = counts_list$counts,
       counts = counts_list$counts,
       length = NULL,
       tx2gene = tx2gene)
}

# Create a Tximport object using the custom function
txi <- custom_tximport(counts_list, tx2gene)

# Summarize the transcript-level data
txi <- summarizeToGene(txi, tx2gene)

# Convert the summarized data to a DESeqDataSet object with the simplified design formula
dds <- DESeqDataSetFromTximport(txi, colData = sample_metadata, design = ~ Tissue + Injection + Feeding + Tissue:Injection + Tissue:Feeding + Injection:Feeding)

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


# Define the comparisons of interest
contrasts <- list(
  "EO_leptin_fooddep_vs_EO_saline_fooddep" = c("Tissue_Interaction_Feeding", "EO_leptin_fooddep", "EO_saline_fooddep"),
  "EO_leptin_adlib_vs_EO_saline_adlib" = c("Tissue_Interaction_Feeding", "EO_leptin_adlib", "EO_saline_adlib"),
  "EO_leptin_fooddep_vs_EO_leptin_adlib" = c("Tissue_Interaction_Feeding", "EO_leptin_fooddep", "EO_leptin_adlib"),
  "SM_leptin_fooddep_vs_SM_saline_fooddep" = c("Tissue_Interaction_Feeding", "SM_leptin_fooddep", "SM_saline_fooddep"),
  "SM_leptin_adlib_vs_SM_saline_adlib" = c("Tissue_Interaction_Feeding", "SM_leptin_adlib", "SM_saline_adlib"),
  "SM_leptin_fooddep_vs_SM_leptin_adlib" = c("Tissue_Interaction_Feeding", "SM_leptin_fooddep", "SM_leptin_adlib"),
  "EO_leptin_fooddep_vs_SM_leptin_fooddep" = c("Tissue_Interaction_Feeding", "EO_leptin_fooddep", "SM_leptin_fooddep"),
  "EO_saline_fooddep_vs_SM_saline_fooddep" = c("Tissue_Interaction_Feeding", "EO_saline_fooddep", "SM_saline_fooddep"),
  "EO_leptin_adlib_vs_SM_leptin_adlib" = c("Tissue_Interaction_Feeding", "EO_leptin_adlib", "SM_leptin_adlib"),
  "EO_saline_adlib_vs_SM_saline_adlib" = c("Tissue_Interaction_Feeding", "EO_saline_adlib", "SM_saline_adlib")
)

# Perform the comparisons and store the results
results_list <- lapply(contrasts, function(contrast) {
  res <- results(dds, contrast = contrast)
  return(res)
})

names(results_list) <- names(contrasts)

# Combine the results into a single data frame
all_res <- bind_rows(results_list, .id = "comparison")

# Write the results to a CSV file
write.csv(all_res, "DESeq2_results.csv")

# Generate the heatmap
rld <- rlog(dds)
select_conditions <- colnames(rld)[colData(rld)$condition %in% c("EO_leptin_fooddep", "EO_saline_fooddep", "EO_leptin_adlib", "EO_saline_adlib",
                                                                 "SM_leptin_fooddep", "SM_saline_fooddep", "SM_leptin_adlib", "SM_saline_adlib")]

annotation <- as.data.frame(colData(rld)[, c("Tissue", "Injection", "Feeding")])
rownames(annotation) <- rownames(colData(rld))
pheatmap(assay(rld[, select_conditions]), annotation_col = annotation, show_colnames = FALSE)

