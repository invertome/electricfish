# Load required libraries
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(heatmap3)

# Read the count matrix and sample metadata files
count_matrix <- read.csv("Count_Matrix.csv", row.names = 1)
sample_metadata <- read.csv("Sample_Metadata.csv")

# Create DGEList object
y <- DGEList(counts = count_matrix, group = factor(paste(sample_metadata$Tissue, sample_metadata$Injection, sample_metadata$Feeding, sep="_")))

# Filter out lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)

# Estimate dispersions
y <- estimateDisp(y)

# Perform quasi-likelihood F-test
fit <- glmQLFit(y, design = model.matrix(~0 + group(y)))
qlf <- glmQLFTest(fit, coef = 2)

# Define the contrasts of interest
contr.matrix <- makeContrasts(
  EO_leptin_fooddep_vs_EO_saline_fooddep = "EO_leptin_fooddep - EO_saline_fooddep",
  EO_leptin_adlib_vs_EO_saline_adlib = "EO_leptin_adlib - EO_saline_adlib",
  EO_leptin_fooddep_vs_EO_leptin_adlib = "EO_leptin_fooddep - EO_leptin_adlib",
  SM_leptin_fooddep_vs_SM_saline_fooddep = "SM_leptin_fooddep - SM_saline_fooddep",
  SM_leptin_adlib_vs_SM_saline_adlib = "SM_leptin_adlib - SM_saline_adlib",
  SM_leptin_fooddep_vs_SM_leptin_adlib = "SM_leptin_fooddep - SM_leptin_adlib",
  EO_leptin_fooddep_vs_SM_leptin_fooddep = "EO_leptin_fooddep - SM_leptin_fooddep",
  EO_saline_fooddep_vs_SM_saline_fooddep = "EO_saline_fooddep - SM_saline_fooddep",
  EO_leptin_adlib_vs_SM_leptin_adlib = "EO_leptin_adlib - SM_leptin_adlib",
  EO_saline_adlib_vs_SM_saline_adlib = "EO_saline_adlib - SM_saline_adlib",
  levels = group(y)
)

# Perform the comparisons and store the results
results_list <- list()
for (cname in colnames(contr.matrix)) {
  contrast <- contr.matrix[, cname]
  qlf_contrast <- glmQLFTest(fit, contrast = contrast)
  res <- topTags(qlf_contrast, n = Inf, adjust.method = "BH")
  results_list[[cname]] <- res
}

# Combine the results into a single data frame
all_res <- bind_rows(results_list, .id = "comparison")

# Write the results to a CSV file
write.csv(all_res, "edgeR_results.csv")

# PCA plot
logcpm <- cpm(y, log = TRUE)
percentVar <- round(100 * attr(logcpm, "percentVar"))

ggplot(logcpm, aes(PC1, PC2, color = Tissue, shape = Injection, linetype = Feeding)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA plot") +
  scale_color_manual(name = "Tissue", values = c("EO" = "blue", "SM" = "red")) +
  scale_shape_manual(name = "Injection", values = c("leptin" = 19, "saline" = 17)) +
  scale_linetype_discrete(name = "Feeding") +
  guides(shape = guide_legend(title = "Injection"), linetype = guide_legend(title = "Feeding")) +
  geom_text_repel(aes(label = rownames(logcpm)), force = 5)
ggsave("PCA_plot.png", width = 10, height = 8)


# Heatmap visualization
# Find top 50 differentially expressed genes across all contrasts
top_genes <- unique(unlist(lapply(results_list, function(x) rownames(head(x$table, 50)))))
top_genes_counts <- count_matrix[top_genes,]

# Create a heatmap
heatmap_data <- log2(top_genes_counts + 1)
colnames(heatmap_data) <- paste(sample_metadata$Tissue, sample_metadata$Injection, sample_metadata$Feeding, sep = "_")
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

heatmap3(heatmap_data, ColSideColors = as.matrix(sample_metadata[, c("Tissue", "Injection", "Feeding")]), col = heatmap_colors, scale = "row", margins = c(5,10), cexRow = 0.8, cexCol = 1.0, labCol = colnames(heatmap_data), hclustfun = function(x) hclust(x, method = "ward.D2"))

# Save the heatmap as a PNG file
ggsave("heatmap.png", width = 1000, height = 1000, dpi = 150)
