library(DESeq2)
library(ggplot2)
# Read the count data and metadata from CSV files
count_data <- read.csv("C:\\Users\\ASUS\\Desktop\\annu\\assignment\\counts3.csv", header=TRUE, row.names=1)
meta_data <- read.csv("C:\\Users\\ASUS\\Desktop\\annu\\assignment\\meta3.csv", header=TRUE, row.names=1)

# Check if the sample names match between count_data and meta_data
if (!all(rownames(meta_data) %in% colnames(count_data))) {
  stop("Sample names in metadata and count data do not match.")
}

# Ensure the column names of count_data match the row names of meta_data
count_data <- count_data[, rownames(meta_data)]

# Check for NA values in the count_data
if (any(is.na(count_data))) {
  warning("NA values detected in count data. Replacing NA values with 0.")
  count_data[is.na(count_data)] <- 0
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta_data, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Save DESeq2 results to CSV
write.csv(as.data.frame(res), file = "deseq2_results.csv")

# Extract upregulated and downregulated genes
upregulated_genes <- subset(res, log2FoldChange > 2 & padj < 0.05)
downregulated_genes <- subset(res, log2FoldChange < (-2) & padj < 0.05)

# Save upregulated and downregulated genes to CSV files
write.csv(as.data.frame(upregulated_genes), file = "upregulated_genes.csv")
write.csv(as.data.frame(downregulated_genes), file = "downregulated_genes.csv")

# Prepare data for volcano plot
res_df <- as.data.frame(res)
res_df$logP <- -log10(res_df$pvalue)
res_df$category <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 2, "Upregulated",
                          ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -2, "Downregulated", "Not Significant"))

# Volcano plot with upregulated genes in blue and downregulated in red
ggplot(res_df, aes(x = log2FoldChange, y = logP, color = category)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "blue", "Downregulated" = "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90"))

# Save the volcano plot
ggsave("volcano_plot.png", width = 8, height = 6)

# Perform variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Perform PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate percent variance explained
percentVar <- attr(pca_data, "percentVar")

# PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", round(percentVar[1], 2), "% variance"),
       y = paste0("PC2: ", round(percentVar[2], 2), "% variance"),
       title = "PCA Plot") +
  theme_minimal()

# Save the PCA plot
ggsave("pca_plot.png", width = 8, height = 6)

