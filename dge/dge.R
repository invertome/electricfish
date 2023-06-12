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

