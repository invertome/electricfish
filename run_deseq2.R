# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Read the count matrix and sample metadata files
count_matrix <- read.csv("Count_Matrix.csv", row.names = 1)
sample_metadata <- read.csv("Sample_Metadata.csv")

# Convert the count matrix to a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_metadata,
                              design = ~ condition)

# Normalize and perform DESeq2 analysis
dds <- DESeq(dds)

# Define the comparisons of interest
contrast1 <- c("condition", "leptin_fooddep", "saline_fooddep")
contrast2 <- c("condition", "leptin_adlib", "saline_adlib")
contrast3 <- c("condition", "leptin_fooddep", "leptin_adlib")
contrast4 <- c("condition", "saline_fooddep", "saline_adlib")

# Perform the comparisons and store the results
res1 <- results(dds, contrast = contrast1)
res2 <- results(dds, contrast = contrast2)
res3 <- results(dds, contrast = contrast3)
res4 <- results(dds, contrast = contrast4)

# Combine the results into a single data frame
all_res <- list("leptin_fooddep_vs_saline_fooddep" = res1,
                "leptin_adlib_vs_saline_adlib" = res2,
                "leptin_fooddep_vs_leptin_adlib" = res3,
                "saline_fooddep_vs_saline_adlib" = res4) %>%
  bind_rows(.id = "comparison")

# Write the results to a CSV file
write.csv(all_res, "DESeq2_results.csv")

# PCA plot
vsd <- vst(dds)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA plot") +
  scale_color_discrete(name = "Condition") +
  guides(shape = guide_legend(title = "Condition")) +
  geom_text_repel(aes(label = rownames(pca_data)), force = 5)

# Save PCA plot
ggsave("PCA_plot.png", width = 10, height = 8)

# Function to create volcano plots
create_volcano_plot <- function(df, title, output_file) {
  # Define threshold lines
  signif_thr <- -log10(0.05)
  log2fc_thr <- 1
  
  df <- df %>% 
    mutate(sig = padj < 0.05 & abs(log2FoldChange) > log2fc_thr)
  
  ggplot(df, aes(log2FoldChange, -log10(padj), color = sig)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("gray", "red")) +
    xlab("log2(Fold Change)") +
    ylab("-log10(Adjusted p-value)") +
    theme_bw() +
    ggtitle(title) +
    geom_hline(yintercept = signif_thr, linetype = "dashed") +
    geom_vline(xintercept = c(-log2fc_thr, log2fc_thr), linetype = "dashed") +
    ggsave(output_file, width = 10, height = 8)
}

# Create and save volcano plots
create_volcano_plot(res1, "Leptin Food Dep vs. Saline Food Dep", "volcano_plot1.png")
create_volcano_plot(res2, "Leptin Adlib vs. Saline Adlib", "volcano_plot2.png")
create_volcano_plot(res3, "Leptin Food Dep vs. Leptin Adlib", "volcano_plot3.png")
create_volcano_plot(res4, "Saline Food Dep vs. Saline Adlib", "volcano_plot4.png")
