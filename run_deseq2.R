# Load required libraries
library(DESeq2)
library(tidyverse)

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
