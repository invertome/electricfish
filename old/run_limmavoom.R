# Load required libraries
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(heatmap3)
library(cowplot)

# Read the count matrix and sample metadata files
count_matrix <- read.csv("Count_Matrix.csv", row.names = 1)
sample_metadata <- read.csv("Sample_Metadata.csv")

# Create DGEList object
y <- DGEList(counts = count_matrix, group = factor(paste(sample_metadata$Tissue, sample_metadata$Injection, sample_metadata$Feeding, sep="_")))

# Filter out lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y, method = "TMM")

# Perform voom transformation
v <- voom(y, normalize.method = "none", plot = TRUE)

# Fit linear model
design <- model.matrix(~0 + group(y))
fit <- lmFit(v, design)
fit2 <- eBayes(fit)

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
  fit2_contrast <- contrasts.fit(fit2, contrast)
  fit2_contrast <- eBayes(fit2_contrast)
  res <- topTable(fit2_contrast, n = Inf, adjust.method = "BH", sort.by = "none")
  results_list[[cname]] <- res
}

# Combine the results into a single data frame
all_res <- bind_rows(results_list, .id = "comparison")

# Write the results to a CSV file
write.csv(all_res, "limma_voom_results.csv")

# Use decideTests to summarize differentially expressed genes across contrasts
decideTests_res <- decideTests(fit2, method = "separate", adjust.method = "BH", p.value = 0.05)
summary(decideTests_res)

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

# Volcano plots for each comparison
volcano_plots <- list()
for (cname in colnames(contr.matrix)) {
  res <- results_list[[cname]]
  volcano_data <- data.frame(logFC = res$logFC, pvalue = res$P.Value, adjpvalue = res$adj.P.Val, Gene = rownames(res))
  volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(pvalue), color = adjpvalue < 0.05)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey", "red"), guide = FALSE) +
    theme_bw() +
    xlab("log2 Fold Change") +
    ylab("-log10 P-value") +
    ggtitle(paste0("Volcano plot: ", cname)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_text_repel(aes(label = ifelse(adjpvalue < 0.05 & (abs(logFC) > 1), Gene, "")), size = 3, box.padding = 0.5)
  volcano_plots[[cname]] <- volcano_plot
  ggsave(paste0("volcano_plot_", cname, ".png"), volcano_plot, width = 10, height = 8)
}

# Bar plot showing the number of differentially expressed genes in each comparison
num_de_genes <- sapply(results_list, function(x) sum(x$adj.P.Val < 0.05))
bar_plot <- ggplot(data.frame(Comparison = names(num_de_genes), Num_DE_Genes = num_de_genes), aes(x = Comparison, y = Num_DE_Genes)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  xlab("Comparison") +
  ylab("Number of Differentially Expressed Genes") +
  ggtitle("Number of Differentially Expressed Genes per Comparison") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("bar_plot_num_de_genes.png", bar_plot, width = 10, height = 8)


# MD plots for each comparison
md_plots <- list()
for (cname in colnames(contr.matrix)) {
  res <- results_list[[cname]]
  md_data <- data.frame(AveExpr = res$AveExpr, logFC = res$logFC, adjpvalue = res$adj.P.Val, Gene = rownames(res))
  md_plot <- plotMD(fit2, contrast = cname, status = decideTests_res[, cname], main = paste0("MD plot: ", cname), xlim = c(min(md_data$AveExpr) - 0.5, max(md_data$AveExpr) + 0.5))
  md_plots[[cname]] <- md_plot
  ggsave(paste0("md_plot_", cname, ".png"), md_plot, width = 10, height = 8)
}

# Heatmap visualization
# Find top 50 differentially expressed genes across all contrasts
top_genes <- unique(unlist(lapply(results_list, function(x) rownames(head(x, 50)))))
top_genes_counts <- count_matrix[top_genes,]

# Create a heatmap
heatmap_data <- log2(top_genes_counts + 1)
colnames(heatmap_data) <- paste(sample_metadata$Tissue, sample_metadata$Injection, sample_metadata$Feeding, sep = "_")
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

heatmap3(heatmap_data, ColSideColors = as.matrix(sample_metadata[, c("Tissue", "Injection", "Feeding")]), col = heatmap_colors, scale = "row", margins = c(5,10), cexRow = 0.8, cexCol = 1.0, labCol = colnames(heatmap_data), hclustfun = function(x) hclust(x, method = "ward.D2"))

# Save the heatmap as a PNG file
ggsave("heatmap.png", width = 1000, height = 1000, dpi = 150)

