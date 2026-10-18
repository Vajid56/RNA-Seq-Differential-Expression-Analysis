library(data.table)

count_dir <- "/home/vajid/InternLab/ncbi_dataset/data/GCF_000001405.40/counts"

count_files <- list.files(count_dir, pattern = "\\.txt$", full.names = TRUE)

count_list <- lapply(count_files, function(file) {
  sample <- gsub("_Counts.txt", "", basename(file))
  counts <- fread(file, header = FALSE, col.names = c("Gene", sample))
  counts <- counts[!duplicated(counts$Gene), ]
  counts
})

merged_counts <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), count_list)
merged_counts[is.na(merged_counts)] <- 0
merged_counts <- merged_counts[!grepl("^__", merged_counts$Gene), ]

fwrite(merged_counts, "merged_counts.csv", row.names = FALSE)

library(DESeq2)

# Correct metadata preparation (2 samples per condition)
colData <- data.frame(
  sampleName = c(
    "SRR16638184",
    "SRR16638186",
    "SRR11806265",
    "SRR11806266",
    "SRR16638202",
    "SRR16638201"
  ),
  condition = factor(
    c(
      "stroke", "stroke",
      "Healthy_control", "Healthy_control",
      "DIABETES_TYPE1", "DIABETES_TYPE1"
    )
  )
)

rownames(colData) <- colData$sampleName
colData$sampleName <- NULL

write.csv(colData, "metadata.csv", row.names = TRUE)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

colData <- colData[colnames(countData), , drop = FALSE]
all(rownames(colData) == colnames(countData))
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ condition
)

exists("DESeqDataSetFromMatrix")
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ condition
)

colData <- colData[colnames(countData), , drop = FALSE]
all(rownames(colData) == colnames(countData))
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ condition
)

countData <- read.csv("merged_counts.csv", row.names = 1, check.names = FALSE)
colnames(countData) <- gsub("_Counts$", "", colnames(countData))

countData <- as.matrix(countData)
storage.mode(countData) <- "integer"
countData[is.na(countData)] <- 0L

colnames(countData) <- tolower(gsub("[^a-zA-Z0-9]", "", colnames(countData)))
colData <- read.csv("metadata.csv", row.names = 1)
colData$condition <- factor(colData$condition)
rownames(colData) <- tolower(gsub("[^a-zA-Z0-9]", "", rownames(colData)))

# Ensure Sample Matching
common_samples <- intersect(colnames(countData), rownames(colData))
countData <- countData[, common_samples, drop = FALSE]
colData <- colData[common_samples, , drop = FALSE]
if (length(common_samples) != ncol(countData)) {
  warning("Some samples are missing in metadata! Check sample names.")
}
# Print sample names in countData and colData
print("Sample names in countData:")
print(colnames(countData))

print("Sample names in colData:")
print(rownames(colData))

# Check for mismatches
setdiff(colnames(countData), rownames(colData))  
setdiff(rownames(colData), colnames(countData))  


countData <- read.csv("merged_counts.csv", row.names = 1, check.names = FALSE)

colnames(countData) <- gsub("_Counts$", "", colnames(countData))
colnames(countData) <- tolower(gsub("[^a-zA-Z0-9]", "", colnames(countData)))

countData <- read.csv("merged_counts.csv", row.names = 1, check.names = FALSE)
colnames(countData) <- gsub("_Counts$", "", colnames(countData))

countData <- as.matrix(countData)
storage.mode(countData) <- "integer"
countData[is.na(countData)] <- 0L

colnames(countData) <- tolower(gsub("[^a-zA-Z0-9]", "", colnames(countData)))
colData <- read.csv("metadata.csv", row.names = 1)
colData$condition <- factor(colData$condition)
rownames(colData) <- tolower(gsub("[^a-zA-Z0-9]", "", rownames(colData)))

# Ensure Sample Matching
common_samples <- intersect(colnames(countData), rownames(colData))
countData <- countData[, common_samples, drop = FALSE]
colData <- colData[common_samples, , drop = FALSE]
if (length(common_samples) != ncol(countData)) {
  warning("Some samples are missing in metadata! Check sample names.")
}
# Print sample names in countData and colData
print("Sample names in countData:")
print(colnames(countData))

print("Sample names in colData:")
print(rownames(colData))

# Check for mismatches
setdiff(colnames(countData), rownames(colData))  
setdiff(rownames(colData), colnames(countData))  

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ condition
)

# Set Healthy_control as reference
dds$condition <- relevel(dds$condition, ref = "Healthy_control")

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes retained after filtering:", nrow(dds), "\n")

# Remove zero-count samples (safety check)
zero_samples <- colSums(counts(dds)) == 0
if(any(zero_samples)) {
  warning("Removed zero-count samples: ",
          paste(names(zero_samples)[zero_samples], collapse = ", "))
  dds <- dds[, !zero_samples]
}


DESeq2 Setup & Filtering 


dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, lfcThreshold = 1)
res <- res[order(res$padj), ]
res <- na.omit(res)
write.csv(as.data.frame(res), "DE_results.csv")
res_stroke <- results(
  dds,
  contrast = c("condition", "stroke", "Healthy_control"),
  alpha = 0.05,
  lfcThreshold = 1
)

res_stroke <- res_stroke[order(res_stroke$padj), ]
res_stroke <- na.omit(res_stroke)

write.csv(as.data.frame(res_stroke), "DE_stroke_vs_healthy.csv")


dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, lfcThreshold = 1)
res <- res[order(res$padj), ]
res <- na.omit(res)
write.csv(as.data.frame(res), "DE_results.csv")
res_diabetes <- results(
  dds,
  contrast = c("condition", "DIABETES_TYPE1", "Healthy_control"),
  alpha = 0.05,
  lfcThreshold = 1
)

res_diabetes <- res_diabetes[order(res_diabetes$padj), ]
res_diabetes <- na.omit(res_diabetes)

write.csv(as.data.frame(res_diabetes), "DE_diabetes_vs_healthy.csv")


##### Volcano Plot 
library(ggplot2)

# Convert DESeq2 results to dataframe                           
res_df <- as.data.frame(res)

# Define significance using padj 
res_df$significance <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,  
  "Significant",
  "Not Significant"
)

# Select the top 10 most significant genes 
top_genes <- head(res_df[order(res_df$padj), ], 10)
top_genes$Gene <- rownames(top_genes)  

# Create volcano plot
volcano_plot <- ggplot(res_df, aes(
  x = log2FoldChange, 
  y = -log10(padj),
  color = significance
)) +  
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(
    name = "Significance",
    values = c("Not Significant" = "gray", "Significant" = "red"),
    labels = c("Not Significant", "Significant (padj < 0.05 & |LFC| > 1)")
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_text(data = top_genes, aes(label = Gene), vjust = -1, size = 3)  

print(volcano_plot)

library(DESeq2)
plotMA(res, main = "MA Plot", ylim = c(-5, 5)) 

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA Plot of RNA-seq Samples",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  )


library(pheatmap)
library(RColorBrewer)

vsd <- vst(dds, blind = FALSE)
vst_counts <- assay(vsd)
res_ordered <- res[order(res$padj), ]
top30_genes <- rownames(res_ordered)[10:40]
colData_heatmap <- colData[colnames(vst_counts), , drop = FALSE]
pheatmap(
  vst_counts[top30_genes, ],
  annotation_col = colData_heatmap,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap of Top 30 Differentially Expressed Genes",
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
)

 
####

library(pheatmap)
library(RColorBrewer)

vsd <- vst(dds, blind = FALSE)
vst_counts <- assay(vsd)

# Order genes by significance
res_ordered <- res[order(res$padj), ]

# Select 30 genes (skip top few if you want)
top30_genes <- rownames(res_ordered)[10:40]

# Match metadata order
colData_heatmap <- colData[colnames(vst_counts), , drop = FALSE]

# Draw heatmap with clarity controls
pheatmap(
  vst_counts[top30_genes, ],
  annotation_col = colData_heatmap,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  # 👇 KEY SETTINGS FOR CLARITY
  fontsize_row = 7,      # smaller gene text
  fontsize_col = 8,     # readable sample text
  cellheight = 8,      # more vertical space
  cellwidth  = 40,      # more horizontal space
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  
  main = "Heatmap of Top 30 Differentially Expressed Genes",
  
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(50)
)

library(DESeq2)
library(apeglm)

# Get variance-stabilized counts
vst_counts <- assay(vsd)

# Ensure colData matches the samples in vst_counts
colData_heatmap <- colData[colnames(vst_counts), , drop = FALSE]

# Select top 30 DEGs
top_genes <- rownames(sig_genes_ppi)[1:min(30, nrow(sig_genes_ppi))]

# Check if fewer than 10 DEGs are available and warn
if (length(top_genes) < 10) {
  warning("Fewer than 10 DEGs available for heatmap!")
}

# Ensure selected genes are present in vst_counts
top_genes <- top_genes[top_genes %in% rownames(vst_counts)]

# Check if we still have genes to plot
if (length(top_genes) > 0) {
  
  # Generate heatmap
  pheatmap(vst_counts[top_genes, ],
           annotation_col = colData_heatmap,  
           scale = "row",   
           cluster_rows = TRUE,     
           cluster_cols = TRUE,     
           show_rownames = TRUE,   
           show_colnames = TRUE,   
           angle_col = 45,         
           main = "Top 30 DEG Expression Profile",
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))
  
} else {
  warning("No valid DEGs available for heatmap plotting.")
}


padj_threshold <- 0.05
lfc_threshold  <- 1
print(top_10_sig)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

up_genes   <- subset(sig_genes, log2FoldChange > 1)
down_genes <- subset(sig_genes, log2FoldChange < -1)
top_10_sig
top_10_up
top_10_down

dds <- DESeq(dds)
res_stroke_vs_control <- results(
  dds,
  contrast = c("condition", "stroke", "Healthy_control")
)

res_control_vs_diabetes <- results(
  dds,
  contrast = c("condition", "Healthy_control", "DIABETES_TYPE1")
)


padj_threshold <- 0.05
lfc_threshold  <- 1

sig_stroke <- subset(
  res_stroke_vs_control,
  padj < padj_threshold & abs(log2FoldChange) > lfc_threshold
)

top_10_stroke <- head(
  sig_stroke[order(sig_stroke$padj), ],
  10
)

cat("\nTop 10 Genes: Stroke vs Control\n")
print(top_10_stroke)

sig_diabetes <- subset(
  res_control_vs_diabetes,
  padj < padj_threshold & abs(log2FoldChange) > lfc_threshold
)

top_10_diabetes <- head(
  sig_diabetes[order(sig_diabetes$padj), ],
  10
)

cat("\nTop 10 Genes: Control vs Diabetes Type 1\n")
print(top_10_diabetes)


sig_diabetes <- subset(
  res_control_vs_diabetes,
  padj < padj_threshold & abs(log2FoldChange) > lfc_threshold
)

top_10_diabetes <- head(
  sig_diabetes[order(sig_diabetes$padj), ],
  10
)

cat("\nTop 10 Genes: Control vs Diabetes Type 1\n")
print(top_10_diabetes)

write.csv(
  as.data.frame(top_10_stroke),
  "Top10_Stroke_vs_Control.csv"
)

write.csv(
  as.data.frame(top_10_diabetes),
  "Top10_Control_vs_DiabetesType1.csv"
)


library(ggplot2)

res_stroke <- results(
  dds,
  contrast = c("condition", "stroke", "Healthy_control")
)

res_stroke_df <- as.data.frame(res_stroke)
res_stroke_df <- res_stroke_df[!is.na(res_stroke_df$padj), ]

res_stroke_df$significance <- ifelse(
  res_stroke_df$padj < 0.05 & abs(res_stroke_df$log2FoldChange) > 1,
  "Significant", "Not Significant"
)

ggplot(res_stroke_df,
       aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Stroke vs Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

## valcono plote for strock vs diatyes type 1

res_diabetes <- results(
  dds,
  contrast = c("condition", "Healthy_control", "DIABETES_TYPE1")
)

res_diabetes_df <- as.data.frame(res_diabetes)
res_diabetes_df <- res_diabetes_df[!is.na(res_diabetes_df$padj), ]

res_diabetes_df$significance <- ifelse(
  res_diabetes_df$padj < 0.05 & abs(res_diabetes_df$log2FoldChange) > 1,
  "Significant", "Not Significant"
)

ggplot(res_diabetes_df,
       aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("gray", "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Control vs Diabetes Type-1",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

library(pheatmap)
library(RColorBrewer)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

# Stroke significant genes
sig_stroke <- subset(
  res_stroke,
  padj < 0.05 & abs(log2FoldChange) > 1
)

top20_stroke <- head(order(sig_stroke$padj), 20)
genes_stroke <- rownames(sig_stroke)[top20_stroke]

pheatmap(
  vst_mat[genes_stroke, ],
  annotation_col = colData,
  scale = "row",
  fontsize_row = 7,
  main = "Top 20 Stroke-specific DEGs",
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
)



library(DESeq2)
# Results: Stroke vs Control
res_stroke <- results(
  dds,
  contrast = c("condition", "stroke", "Healthy_control")
)

# MA plot
plotMA(
  res_stroke,
  main = "MA Plot: Stroke vs Control",
  ylim = c(-5, 5)
)

# Results: Control vs Diabetes Type-1
res_diabetes <- results(
  dds,
  contrast = c("condition", "Healthy_control", "DIABETES_TYPE1")
)

# MA plot
plotMA(
  res_diabetes,
  main = "MA Plot: Control vs Diabetes Type-1",
  ylim = c(-5, 5)
)

library(apeglm)

# Stroke vs Control (shrunk)
res_stroke_shrink <- lfcShrink(
  dds,
  contrast = c("condition", "stroke", "Healthy_control"),
  type = "apeglm"
)

plotMA(
  res_stroke_shrink,
  main = "MA Plot (LFC Shrunk): Stroke vs Control",
  ylim = c(-5, 5)
)

library(DESeq2)
res_stroke <- results(
  dds,
  contrast = c("condition", "stroke", "Healthy_control")
)
exists("results")


library(ggplot2)

res_stroke <- results(
  dds,
  contrast = c("condition", "stroke", "Healthy_control")
)

res_stroke_df <- as.data.frame(res_stroke)
res_stroke_df <- res_stroke_df[!is.na(res_stroke_df$padj), ]

res_stroke_df$significance <- ifelse(
  res_stroke_df$padj < 0.05 & abs(res_stroke_df$log2FoldChange) > 1,
  "Significant", "Not Significant"
)

ggplot(res_stroke_df,
       aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Stroke vs Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

# Results: Stroke vs Control
res_stroke <- results(
  dds,
  contrast = c("condition", "stroke", "Healthy_control")
)

# MA plot
plotMA(
  res_stroke,
  main = "MA Plot: Stroke vs Control",
  ylim = c(-5, 5)
)



# Results: Control vs Diabetes Type-1
res_diabetes <- results(
  dds,
  contrast = c("condition", "Healthy_control", "DIABETES_TYPE1")
)

# MA plot
plotMA(
  res_diabetes,
  main = "MA Plot: Control vs Diabetes Type-1",
  ylim = c(-5, 5)
)



vsd <- vst(dds, blind = FALSE)
samples_stroke_control <- colData$condition %in% c("stroke", "Healthy_control")

vsd_sc <- vsd[, samples_stroke_control]

colData_sc <- colData[samples_stroke_control, , drop = FALSE]
library(DESeq2)
library(ggplot2)

## valcono plot 


  
  library(ggplot2)
  
  res_stroke <- results(
    dds,
    contrast = c("condition", "stroke", "Healthy_control")
  )
  
  res_stroke_df <- as.data.frame(res_stroke)
  res_stroke_df <- res_stroke_df[!is.na(res_stroke_df$padj), ]
  
  res_stroke_df$significance <- ifelse(
    res_stroke_df$padj < 0.05 & abs(res_stroke_df$log2FoldChange) > 1,
    "Significant", "Not Significant"
  )
  
  ggplot(res_stroke_df,
         aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("gray", "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = "Volcano Plot: Stroke vs Control",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    )
  
  res_diabetes <- results(
    dds,
    contrast = c("condition", "Healthy_control", "DIABETES_TYPE1")
  )
  
  res_diabetes_df <- as.data.frame(res_diabetes)
  res_diabetes_df <- res_diabetes_df[!is.na(res_diabetes_df$padj), ]
  
  res_diabetes_df$significance <- ifelse(
    res_diabetes_df$padj < 0.05 & abs(res_diabetes_df$log2FoldChange) > 1,
    "Significant", "Not Significant"
  )
  
  ggplot(res_diabetes_df,
         aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("gray", "blue")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = "Volcano Plot: Control vs Diabetes Type-1",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    )
  

  # PCA PLOT STROKE VS CONTROL 
  
  library(DESeq2)
  library(ggplot2)
  
  vsd <- vst(dds, blind = FALSE)
  
  samples_stroke_control <- colData$condition %in% c("stroke", "Healthy_control")
  
  vsd_sc <- vsd[, samples_stroke_control]
  colData_sc <- colData[samples_stroke_control, , drop = FALSE]
  
  pca_sc <- plotPCA(vsd_sc, intgroup = "condition", returnData = TRUE)
  percentVar_sc <- round(100 * attr(pca_sc, "percentVar"))
  
  ggplot(pca_sc, aes(PC1, PC2, color = condition)) +
    geom_point(size = 5) +
    theme_minimal() +
    labs(
      title = "PCA Plot: Stroke vs Control",
      x = paste0("PC1: ", percentVar_sc[6], "% variance"),
      y = paste0("PC2: ", percentVar_sc[6], "% variance")
    )
  samples_control_diabetes <- colData$condition %in% c("Healthy_control", "DIABETES_TYPE1")
  
  vsd_cd <- vsd[, samples_control_diabetes]
  colData_cd <- colData[samples_control_diabetes, , drop = FALSE]
  
  ## PCA PLOT CONTROL VS DIABETES TYPE 1
  
  samples_control_diabetes <- colData$condition %in% c("Healthy_control", "DIABETES_TYPE1")
  
  vsd_cd <- vsd[, samples_control_diabetes]
  colData_cd <- colData[samples_control_diabetes, , drop = FALSE]
  pca_cd <- plotPCA(vsd_cd, intgroup = "condition", returnData = TRUE)
  percentVar_cd <- round(100 * attr(pca_cd, "percentVar"))
  ggplot(pca_cd, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA Plot: Control vs Diabetes Type-1",
    x = paste0("PC1: ", percentVar_cd[1], "% variance"),
    y = paste0("PC2: ", percentVar_cd[2], "% variance")
  )

##### heat map top 30 DEGS 
  
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)
  res_stroke <- results(
    dds,
    contrast = c("condition", "stroke", "Healthy_control")
  )
  
  res_stroke <- res_stroke[order(res_stroke$padj), ]
  res_stroke <- na.omit(res_stroke)
  
  top30_stroke <- rownames(res_stroke)[1:30]
  samples_sc <- colData$condition %in% c("stroke", "Healthy_control")
  mat_sc <- vst_mat[top30_stroke, samples_sc]
  anno_sc <- colData[samples_sc, , drop = FALSE]
  pheatmap(
    mat_sc,
    annotation_col = anno_sc,
    scale = "row",
    fontsize_row = 6,
    fontsize_col = 8,
    main = "Heatmap: Top 30 DEGs (Stroke vs Control)",
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(70)
  )
  
  
  ### heatmap s=control vs diabtesi type 1 
  res_diabetes <- results(
    dds,
    contrast = c("condition", "Healthy_control", "DIABETES_TYPE1")
  )
  
  res_diabetes <- res_diabetes[order(res_diabetes$padj), ]
  res_diabetes <- na.omit(res_diabetes)
  
  top30_diabetes <- rownames(res_diabetes)[1:30]
  samples_cd <- colData$condition %in% c("Healthy_control", "DIABETES_TYPE1")
  mat_cd <- vst_mat[top30_diabetes, samples_cd]
  anno_cd <- colData[samples_cd, , drop = FALSE]
  pheatmap(
    mat_cd,
    annotation_col = anno_cd,
    scale = "row",
    fontsize_row = 7,
    fontsize_col = 10,
    main = "Heatmap: Top 30 DEGs (Control vs Diabetes Type-1)",
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
  )
  
  ## top 10 gene 
  dds <- DESeq(dds)
  
  res_stroke_vs_control <- results(
    dds,
    contrast = c("condition", "stroke", "Healthy_control")
  )
  
  res_control_vs_diabetes <- results(
    dds,
    contrast = c("condition", "Healthy_control", "DIABETES_TYPE1")
  )
  
  padj_threshold <- 0.05
  lfc_threshold  <- 1
  
  sig_stroke <- subset(
    res_stroke_vs_control,
    padj < padj_threshold & abs(log2FoldChange) > lfc_threshold
  )
  
  top_10_stroke <- head(
    sig_stroke[order(sig_stroke$padj), ],
    10
  )
  
  cat("\nTop 10 Genes: Stroke vs Control\n")
  print(top_10_stroke)
  sig_diabetes <- subset(
    res_control_vs_diabetes,
    padj < padj_threshold & abs(log2FoldChange) > lfc_threshold
  )
  
  top_10_diabetes <- head(
    sig_diabetes[order(sig_diabetes$padj), ],
    10
  )
  
  cat("\nTop 10 Genes: Control vs Diabetes Type 1\n")
  print(top_10_diabetes)
  
  top_10_sig
  top_10_up
  top_10_down
 
  
  #### ppt 
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install(c(
    "STRINGdb",
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot"
  ))
  
  install.packages(c("igraph", "ggraph", "ggplot2", "dplyr"))
  
  
  library(STRINGdb)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(dplyr)
  
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install(c(
    "STRINGdb",
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot"
  ))
  
  install.packages(c("igraph", "ggraph", "ggplot2", "dplyr"))
 
  string_db <- STRINGdb$new(
    version = "11.5",
    species = 9606,   # Human
    score_threshold = 400
  )
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("STRINGdb")
  
 ##
  deg_df <- data.frame(gene = deg_genes)
  
  mapped <- string_db$map(
    deg_df,
    "gene",
    removeUnmappedRows = TRUE
  )
  
  ppi <- string_db$get_interactions(mapped$STRING_id)
  
ppi_graph <- graph_from_data_frame(
  ppi[, c("from", "to")],
  directed = FALSE
)  

ggraph(ppi_graph, layout = "fr")
 geom_edge_link(alpha = 0.3)
 geom_node_point(ase(size = ))
 
 ppi_graph <- graph_from_data_frame(
   ppi[, c("from", "to")],
   directed = FALSE
 )
 
 ggraph(ppi_graph, layout = "fr") +
   geom_edge_link(alpha = 0.3) +
   geom_node_point(aes(size = degree(ppi_graph)), color = "tomato") +
   geom_node_text(aes(label = name), size = 2, repel = TRUE) +
   theme_void() +
   labs(title = "Protein–Protein Interaction Network")
 
 
 
 
 library(DESeq2)
 library(STRINGdb)
 library(igraph)
 
 res_df <- as.data.frame(res)
 
 deg_genes <- rownames(
   subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
 )
 
 length(deg_genes)
 string_db <- STRINGdb$new(
   version = "11.5",
   species = 9606,      # Human
   score_threshold = 400
 )
 mapped <- string_db$map(
   data.frame(gene = deg_genes),
   "gene",
   removeUnmappedRows = TRUE
 )
 
 ppi <- string_db$get_interactions(mapped$STRING_id)
 head(ppi)
 
 g <- graph_from_data_frame(ppi[, c("from", "to")], directed = FALSE)
 
 plot(
   g,
   vertex.size = 5,
   vertex.label = NA,
   main = "Protein–Protein Interaction Network"
 )
 
 
 ppi_clean <- na.omit(ppi[, c("from", "to")])
 head(ppi)
 ppi_clean <- na.omit(ppi[, c("from", "to")])
 any(is.na(ppi_clean))
 # should be FALSE
 library(igraph)
 
 g <- graph_from_data_frame(ppi_clean, directed = FALSE)
 plot(
   g,
   vertex.size = 6,
   vertex.label = NA,
   vertex.color = "skyblue",
   edge.color = "gray",
   main = "Protein–Protein Interaction Network (STRING)"
 )
 set.seed(123)
 plot(
   g,
   layout = layout_with_fr,
   vertex.size = 6,
   vertex.label = NA,
   main = "PPI Network (Force-directed layout)"
 )
 
 
 
 
 
  # Convert DESeq results to data frame
  res_df <- as.data.frame(res)
  
  # Filter significant DEGs
  deg_genes <- rownames(
    subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
  )
  
  # Convert DESeq results to data frame
  res_df <- as.data.frame(res)
  
  # Filter significant DEGs
  deg_genes <- rownames(
    subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
  )
  
  # Check
  length(deg_genes)
  head(deg_genes)
  
  library(STRINGdb)
  
  string_db <- STRINGdb$new(
    version = "11.5",
    species = 9606,
    score_threshold = 400
  )
  
  # Map gene symbols to STRING IDs
  mapped <- string_db$map(
    data.frame(gene = deg_genes),
    "gene",
    removeUnmappedRows = TRUE
  )
  
  # Create PPI network
  string_db$plot_network(mapped$STRING_id)
  
  deg_genes <- rownames(
    subset(res_df, padj < 0.1 & abs(log2FoldChange) > 0.5)
  )
  

  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  gene_symbols <- top_hub_genes$Gene
  gene_symbols <- deg_genes
  
  
  sig_genes <- res[res$padj < 0.05 & abs(res$log2FoldChange) >= 1, ]
  nrow(sig_genes)
  
  up_genes <- res[res$padj < 0.05 & res$log2FoldChange >= 1, ]
  nrow(up_genes)
  down_genes <- res[res$padj < 0.05 & res$log2FoldChange <= -1, ]
  nrow(down_genes)
  
  top_hub_genes <- head(sig_genes[order(sig_genes$padj), ], 10)
  gene_symbols <- rownames(top_hub_genes)
  ls()
  
  res_clean <- na.omit(res)
  
  sig_genes <- res_clean[
    res_clean$padj < 0.05 &
      abs(res_clean$log2FoldChange) >= 1,
  ]
  
  nrow(sig_genes)
  head(sig_genes)

  # Remove NA
  res_clean <- na.omit(res)
  
  # Significant genes
  sig_genes <- res_clean[
    res_clean$padj < 0.05 &
      abs(res_clean$log2FoldChange) >= 1,
  ]
  
  # Top 10 genes
  top_hub_genes <- head(sig_genes[order(sig_genes$padj), ], 10)
  
  # Gene names
  gene_symbols <- rownames(top_hub_genes)
  nrow(sig_genes)
  head(sig_genes)
  head(sig_genes, 20)
  top_hub_genes
  View(top_hub_genes)
  print(gene_symbols)
  up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
  nrow(up_genes)
  head(up_genes)
  down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]
  nrow(down_genes)
  head(down_genes)
  write.csv(sig_genes, "significant_genes.csv")
  write.csv(top_hub_genes, "top_10_genes.csv")
  getwd()
  
  
  
  # Total significant
  nrow(sig_sc)
  
  # Upregulated in Stroke
  up_sc <- sig_sc[sig_sc$log2FoldChange > 0, ]
  nrow(up_sc)
  
  # Downregulated in Stroke
  down_sc <- sig_sc[sig_sc$log2FoldChange < 0, ]
  nrow(down_sc)
  
  
  
  # NA remove
  res_sc_clean <- na.omit(res_sc)
  
  # Significant genes (Stroke vs Control)
  sig_sc <- res_sc_clean[
    res_sc_clean$padj < 0.05 &
      abs(res_sc_clean$log2FoldChange) >= 1,
  ]
  

  only_stroke <- setdiff(genes_sc, union(genes_cd, genes_sd))
  
  
  nrow(sig_sc)
  head(sig_sc)
  
  
  # Stroke vs Control DESeq2 result
  res_sc <- results(
    dds_sc,
    contrast = c("condition", "stroke", "Healthy_control")
  )
  library(DESeq2)installed.packages()[,"Package"]

  
  # 1. Load library
  library(DESeq2)
  
  # 2. Check
  sessionInfo()
  
  
  res_sc <- results(
    dds_sc,
    contrast = c("condition", "stroke", "Healthy_control")
  )
  
  sig_sc   # Stroke vs Control significant genes
  sig_cd   # Control vs Diabetes significant genes
  sig_sd   # Stroke vs Diabetes significant genes
 
  
  
  
  