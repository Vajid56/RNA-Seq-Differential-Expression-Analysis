library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(data.table)

##### 1.Load & Merge Counts
library(data.table)

count_dir <- "/home/vajid/InternLab/ncbi_dataset/data/GCF_000001405.40/counts"

count_files <- list.files(count_dir, pattern = "_counts.txt$", full.names = TRUE)

count_list <- lapply(count_files, function(file) {
  
  sample <- gsub("_counts.txt", "", basename(file))
  
  counts <- fread(file, header = FALSE)
  colnames(counts) <- c("Gene", sample)
  
  counts <- counts[!duplicated(Gene)]
  
  counts
})

merged_counts <- Reduce(function(x,y) merge(x,y,by="Gene", all=TRUE), count_list)

merged_counts[is.na(merged_counts)] <- 0

merged_counts <- merged_counts[!grepl("^__", merged_counts$Gene),]

write.csv(merged_counts, "merged_counts.csv", row.names = FALSE)


### 2. metadata prepration 
  colData <- data.frame(
    sampleName = c(
      "SRR6216693",
      "SRR6216695",
      "SRR6216696",
      "SRR16638184",
      "SRR16638188",
      "SRR16638189",
      "SRR16638201",
      "SRR16638202",
      "SRR16638205"
    ),
    condition = factor(
      c(
        "DIABETES_TYPE1","DIABETES_TYPE1","DIABETES_TYPE1",
        "stroke","stroke","stroke",
        "Healthy_control","Healthy_control","Healthy_control"
      )
    )
  )

rownames(colData) <- colData$sampleName
colData$sampleName <- NULL

write.csv(colData, "metadata.csv", row.names = TRUE)


  ##### 3. load data 
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
if(any(zero_samples)) 
  warning("Removed zero-count samples: ",
          paste(names(zero_samples)[zero_samples], collapse = ", "))
  dds <- dds[, !zero_samples]
  

    
    # Run DESeq
    dds <- DESeq(dds)
  
  # Overall results
  res <- results(dds, alpha = 0.05, lfcThreshold = 1)
  
  # Sort by adjusted p-value
  res <- res[order(res$padj), ]
  
  # Remove NA values
  res <- na.omit(res)
  
  # Save results
  write.csv(as.data.frame(res), "DE_results.csv")
  
  
  # Stroke vs Healthy comparison
  res_stroke <- results(
    dds,
    contrast = c("condition", "stroke", "Healthy_control"),
    alpha = 0.05,
    lfcThreshold = 1
  )
  
  # Sort genes
  res_stroke <- res_stroke[order(res_stroke$padj), ]
  
  # Remove NA
  res_stroke <- na.omit(res_stroke)
  
  # Save comparison results
  write.csv(as.data.frame(res_stroke), "DE_stroke_vs_healthy.csv")
  

    # Diates vs control
  
    res_diabetes <- results(
      dds,
      contrast = c("condition", "DIABETES_TYPE1", "Healthy_control")
    )
  
  res_diabetes <- res_diabetes[order(res_diabetes$padj), ]
  res_diabetes <- na.omit(res_diabetes)
  
  write.csv(as.data.frame(res_diabetes), "DE_diabetes_vs_healthy.csv")



## only stroke 
  # Run DESeq
  dds <- DESeq(dds)

# Extract only stroke samples
stroke_samples <- colData$condition == "stroke"

stroke_counts <- countData[, stroke_samples]

# Save stroke counts
write.csv(stroke_counts, "only_stroke_samples.csv")


  # only Diabtes sample 
  # Diabetes samples identify
  diabetes_samples <- colData$condition == "DIABETES_TYPE1"

# Extract diabetes counts
diabetes_counts <- countData[, diabetes_samples]

# Save file
write.csv(diabetes_counts, "only_diabetes_samples.csv")

#####################
  # Run DESeq 
  dds <- DESeq(dds)

# Stroke vs Diabetes comparison
res_stroke_vs_diabetes <- results(
  dds,
  contrast = c("condition", "stroke", "DIABETES_TYPE1")
)

# Sort by adjusted p-value
res_stroke_vs_diabetes <- res_stroke_vs_diabetes[order(res_stroke_vs_diabetes$padj), ]

# Remove NA values
res_stroke_vs_diabetes <- na.omit(res_stroke_vs_diabetes)

# Save result file
write.csv(as.data.frame(res_stroke_vs_diabetes), "DE_stroke_vs_diabetes.csv")

############################################################################################
  
  # DIFFERENTIAL EXPRESION ANYALIS 
  # VISUALIZTION PLOT 
  
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
library(ggplot2)
library(ggrepel)

# Stroke vs Control results
res_stroke <- results(dds, contrast = c("condition","stroke","Healthy_control"))

# Convert to dataframe
res_df <- as.data.frame(res_stroke)

# Remove NA values
res_df <- res_df[!is.na(res_df$padj), ]

# Add gene names
res_df$Gene <- rownames(res_df)

# Define significance
res_df$significance <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)

# Select top 30 genes
top_genes <- head(res_df[order(res_df$padj), ], 30)

# Volcano Plot
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
  geom_point(alpha=0.6, size=2) +
  geom_text_repel(data=top_genes, aes(label=Gene), size=3) +
  scale_color_manual(values=c("gray","red")) +
  theme_minimal() +
  labs(
    title="Volcano Plot: Stroke vs Healthy Control",
    x="Log2 Fold Change",
    y="-Log10 Adjusted P-value"
  ) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed")

#####################3 valcno plote stroke vs Diabtes 
library(ggplot2)
library(ggrepel)

# Stroke vs Diabetes results
res_sd <- results(
  dds,
  contrast = c("condition","stroke","DIABETES_TYPE1")
)

# Convert to dataframe
res_df <- as.data.frame(res_sd)

# Remove NA
res_df <- res_df[!is.na(res_df$padj), ]

# Add gene names
res_df$Gene <- rownames(res_df)

# Define significance
res_df$significance <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)

# Select top 30 genes
top_genes <- head(res_df[order(res_df$padj), ], 30)

# Volcano plot
volcano_plot <- ggplot(res_df, aes(
  x = log2FoldChange, 
  y = -log10(padj),
  color = significance
)) +  
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(
    name = "Significance",
    values = c("Not Significant" = "gray", "Significant" = "blue")
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: Stroke vs Diabetes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3)

print(volcano_plot)


library(DESeq2)
library(ggplot2)
library(ggrepel)

# Diabetes vs Control results
res_diabetes <- results(
  dds,
  contrast = c("condition","DIABETES_TYPE1","Healthy_control")
)

# Convert to dataframe
res_df <- as.data.frame(res_diabetes)

# Remove NA values
res_df <- res_df[!is.na(res_df$padj), ]

# Add gene names
res_df$Gene <- rownames(res_df)

# Define significance
res_df$significance <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)

# Select top 10 genes
top_genes <- head(res_df[order(res_df$padj), ], 10)

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3) +
  scale_color_manual(values = c("gray","red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Diabetes vs Healthy Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

##############################
library(ggplot2)
library(ggrepel)

# Stroke vs Healthy Control results
res_stroke <- results(dds, contrast = c("condition", "stroke", "Healthy_control"))

# Convert to dataframe
res_df <- as.data.frame(res_stroke)

# Remove NA
res_df <- res_df[!is.na(res_df$padj), ]

# Define significance
res_df$significance <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)

# Select top 10 genes
top_genes <- head(res_df[order(res_df$padj), ], 10)
top_genes$Gene <- rownames(top_genes)

# Create volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: Stroke vs Healthy Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")


##################################A
  
    plotMA(res,
         main = "MA Plot",
         ylim = c(-5,5),
         alpha = 0.05)
    
#######################################333
 ## PCA PLOT FOR STROKE VS CONTROL VS DIABETEA TYPE ONE    
if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install("DESeq2")
    
    vsd <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pca_data, "percentVar"))
    ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = 5) +
      theme_minimal() +
      labs(
        title = "PCA Plot of RNA-seq Samples",
        x = paste0("PC1: ", percentVar[1], "% variance"),
        y = paste0("PC2: ", percentVar[2], "% variance")
      ) 

    plot_pca_condition(vsd, "Healthy_control", "stroke")  

    plot_pca_condition <- function(vsd, condition1, condition2) {
      
      # Subset samples
      keep <- vsd$condition %in% c(condition1, condition2)
      
      vsd_sub <- vsd[, keep]
      
      # PCA data
      pca_data <- plotPCA(vsd_sub, intgroup = "condition", returnData = TRUE)
      
      percentVar <- round(100 * attr(pca_data, "percentVar"))
      
      # Plot
      ggplot(pca_data, aes(PC1, PC2, color = condition)) +
        geom_point(size = 4) +
        theme_minimal(base_size = 14) +
        labs(
          title = paste(condition1, "vs", condition2),
          x = paste0("PC1: ", percentVar[1], "%"),
          y = paste0("PC2: ", percentVar[2], "%")
        )
    }
    plot_pca_condition(vsd, "stroke", "DIABETES_TYPE1") 
    plot_pca_condition(vsd, "Healthy_control", "DIABETES_TYPE1")

##################################################
    
    # ================================
    # HEATMAP - Clean & Attractive
    # Load libraries
    library(DESeq2)
    library(pheatmap)
    library(RColorBrewer)
    
    # Variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    
    # Extract normalized counts
    vst_counts <- assay(vsd)
    
    # Remove NA padj values
    res_clean <- res[!is.na(res$padj), ]
    
    # Order genes by significance
    res_ordered <- res_clean[order(res_clean$padj), ]
    
    # Select top 30 genes safely
    top30_genes <- rownames(res_ordered)[1:min(30, nrow(res_ordered))]
    
    # Ensure genes exist in matrix
    top30_genes <- top30_genes[top30_genes %in% rownames(vst_counts)]
    
    # Match metadata order
    colData_heatmap <- colData[colnames(vst_counts), , drop = FALSE]
    
    # Draw heatmap
    pheatmap(
      vst_counts[top30_genes, ],
      annotation_col = colData_heatmap,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      
      fontsize_row = 6,
      fontsize_col = 7,
      cellheight = 8,
      cellwidth  = 40,
      
      show_rownames = TRUE,
      show_colnames = TRUE,
      
      main = "Heatmap of Top 30 Differentially Expressed Genes",
      
      color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(50)
    )
    

    
    cat("✅ Saved: Heatmap_top30.png\n")

    ########
    
    # Load libraries
    library(DESeq2)
    library(pheatmap)
    library(RColorBrewer)
    
    # Variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    vst_counts <- assay(vsd)
    
    # Remove NA padj genes
    res_clean <- res[!is.na(res$padj), ]
    res_ordered <- res_clean[order(res_clean$padj), ]
    
    # Select top 30 genes safely
    top30_genes <- rownames(res_ordered)[1:min(30, nrow(res_ordered))]
    top30_genes <- top30_genes[top30_genes %in% rownames(vst_counts)]
    
    # Match metadata order
    colData_heatmap <- colData[colnames(vst_counts), , drop = FALSE]
    
    # Define bright colors for conditions
    annotation_colors <- list(
      condition = c(
        stroke = "#FF5733",          # bright orange-red
        Healthy_control = "#33C1FF", # bright blue
        DIABETES_TYPE1 = "#33FF99"   # bright green
      )
    )
    
    # Bright heatmap color palette
    heat_colors <- colorRampPalette(c("#053061", "#2166AC", "#F7F7F7", "#B2182B", "#67001F"))(100)
    
    # Save high resolution image
    png("Attractive_Top30_DEG_heatmap.png", width = 2200, height = 2600, res = 300)
    
    pheatmap(
      vst_counts[top30_genes, ],
      annotation_col = colData_heatmap,
      annotation_colors = annotation_colors,
      
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      
      fontsize_row = 7,
      fontsize_col = 10,
      
      cellheight = 12,
      cellwidth  = 35,
      
      border_color = NA,
      show_rownames = TRUE,
      show_colnames = TRUE,
      
      main = "Top 30 Differentially Expressed Genes",
      
      color = heat_colors
    )
    
    dev.off()
    
    
  #### heat map top 10\
  
    # Load DESeq2
    library(DESeq2)
    
    # Get results for Stroke vs Control
    res_stroke <- results(
      dds,
      contrast = c("condition", "stroke", "Healthy_control")
    )
    
    # Remove NA values
    res_stroke <- res_stroke[!is.na(res_stroke$padj), ]
    
    # MA plot
    plotMA(
      res_stroke,
      main = "MA Plot: Stroke vs Control",
      ylim = c(-5, 5),
      alpha = 0.05
    )
    
    # Add horizontal line for reference
    abline(h = 0, col = "blue", lwd = 2)
    
   
    plotMA(res_stroke)
    
#### MA PLOT STROKE VS DIABTES 
    
    res_stroke_diabetes <- results(
      dds,
      contrast = c("condition", "stroke", "DIABETES_TYPE1")
    )
    
    # remove NA
    res_stroke_diabetes <- res_stroke_diabetes[!is.na(res_stroke_diabetes$padj), ]
    
    # MA plot
    plotMA(
      res_stroke_diabetes,
      main = "MA Plot: Stroke vs Diabetes_type1",
      ylim = c(-5,5)
      
  
    )
    
#### 3 PLOT generete together ############################  
    library(DESeq2)
    
    # Get results
    res_stroke_control <- results(dds, contrast = c("condition","stroke","Healthy_control"))
    res_diabetes_control <- results(dds, contrast = c("condition","DIABETES_TYPE1","Healthy_control"))
    res_stroke_diabetes <- results(dds, contrast = c("condition","stroke","DIABETES_TYPE1"))
    
    # Remove NA
    res_stroke_control <- res_stroke_control[!is.na(res_stroke_control$padj), ]
    res_diabetes_control <- res_diabetes_control[!is.na(res_diabetes_control$padj), ]
    res_stroke_diabetes <- res_stroke_diabetes[!is.na(res_stroke_diabetes$padj), ]
    
    # Create 3 MA plots together
    par(mfrow=c(1,3))   # 1 row 3 plots
    
    plotMA(res_stroke_control,
           main="Stroke vs Control",
           ylim=c(-5,5))
    
    plotMA(res_diabetes_control,
           main="Diabetes vs Control",
           ylim=c(-5,5))
    
    plotMA(res_stroke_diabetes,
           main="Stroke vs Diabetes",
           ylim=c(-5,5))
    
    # Reset plotting layout
    par(mfrow=c(1,1))
    
###### olny  diabtes ma plote 
    library(DESeq2)
    
    # Diabetes vs Healthy Control results
    res_diabetes <- results(
      dds,
      contrast = c("condition", "DIABETES_TYPE1", "Healthy_control")
    )
    
    # Remove NA values
    res_diabetes <- res_diabetes[!is.na(res_diabetes$padj), ]
    
    # MA Plot
    plotMA(
      res_diabetes,
      main = "MA Plot: Diabetes Type-1 vs Healthy Control",
      ylim = c(-5, 5)
    )
    
    
################3 top 10 gene heatmap stroke vs diabtes 
  
    library(DESeq2)
    library(pheatmap)
    
    # Stroke vs Diabetes results
    res_stroke_diabetes <- results(
      dds,
      contrast = c("condition", "stroke", "DIABETES_TYPE1")
    )
    
    # Remove NA values
    res_stroke_diabetes <- res_stroke_diabetes[!is.na(res_stroke_diabetes$padj), ]
    
    # Order genes by significance
    res_ordered <- res_stroke_diabetes[order(res_stroke_diabetes$padj), ]
    
    # Select top 10 genes
    top10_genes <- rownames(res_ordered)[1:10]
    
    # Variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    vst_counts <- assay(vsd)
    
    # Heatmap
    pheatmap(
      vst_counts[top10_genes, ],
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      main = "Top 10 DEGs: Stroke vs Diabetes"
    )
    
    
    
############## TOP 10 DIFFERENTIALY EXPRESSED GENE DIABETES VA STROKE ###################
    
    library(pheatmap)
    library(RColorBrewer)
    
    # Variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    vst_counts <- assay(vsd)
    
    # Stroke vs Diabetes results
    res_sd <- results(dds, contrast = c("condition","stroke","DIABETES_TYPE1"))
    
    # Remove NA
    res_sd <- res_sd[!is.na(res_sd$padj), ]
    
    # Order genes by significance
    res_ordered <- res_sd[order(res_sd$padj), ]
    
    # Select TOP 10 genes
    top10_genes <- rownames(res_ordered)[1:10]
    
    # Select only stroke and diabetes samples
    samples_sd <- rownames(colData)[colData$condition %in% c("stroke","DIABETES_TYPE1")]
    
    # Heatmap
    pheatmap(
      vst_counts[top10_genes, samples_sd],
      annotation_col = colData[samples_sd, , drop = FALSE],
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      
      fontsize_row = 8,
      fontsize_col = 8,
      cellheight = 20,
      cellwidth  = 40,
      
      show_rownames = TRUE,
      show_colnames = TRUE,
      
      main = "Top 10 Differentially Expressed Genes: Stroke vs Diabetes",
      
      color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(50)
    )
    
############## heatmap stroke vs control  
    library(pheatmap)
    library(RColorBrewer)
    
    # Variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    vst_counts <- assay(vsd)
    
    # Stroke vs Diabetes results
    res_sd <- results(dds, contrast = c("condition","stroke","Healthy_control"))
    
    # Remove NA
    res_sd <- res_sd[!is.na(res_sd$padj), ]
    
    # Order genes by significance
    res_ordered <- res_sd[order(res_sd$padj), ]
    
    # Select TOP 10 genes
    top10_genes <- rownames(res_ordered)[1:10]
    
    # Select only stroke and diabetes samples
    samples_sd <- rownames(colData)[colData$condition %in% c("stroke","Healthy_control")]
    
    # Heatmap
    pheatmap(
      vst_counts[top10_genes, samples_sd],
      annotation_col = colData[samples_sd, , drop = FALSE],
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      
      fontsize_row = 8,
      fontsize_col = 8,
      cellheight = 20,
      cellwidth  = 40,
      
      show_rownames = TRUE,
      show_colnames = TRUE,
      
      main = "Top 10 Differentially Expressed Genes: Stroke vs Healthy_control",
      
      color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(50)
    )
#######################  TOP 10 DIFFERENTIALY EXPRESSED GENE DIABETES VA STROKE  ############
    library(pheatmap)
    library(RColorBrewer)
    
    # Variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    vst_counts <- assay(vsd)
    
    # Diabetes vs Stroke results
    res_ds <- results(dds, contrast = c("condition","stroke","Healthy_control"))
    
    # Remove NA
    res_ds <- res_ds[!is.na(res_ds$padj), ]
    
    # Order genes by significance
    res_ordered <- res_ds[order(res_ds$padj), ]
    
    # Select top 10 genes
    top10_genes <- rownames(res_ordered)[1:10]
    
    # Select only diabetes and stroke samples
    samples_ds <- rownames(colData)[colData$condition %in% c("stroke", "Healthy_control")]
    
    # Heatmap
    pheatmap(
      vst_counts[top10_genes, samples_ds],
      annotation_col = colData[samples_ds, , drop = FALSE],
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      
      fontsize_row = 8,
      fontsize_col = 8,
      cellheight = 20,
      cellwidth  = 40,
      
      show_rownames = TRUE,
      show_colnames = TRUE,
      
      main = "Top 10 Differentially Expressed Genes: stroke vs Healthy_control",
      
      color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(50)
    )
    
    
    
########################## top 10 gene all contion ############
  
    library(pheatmap)
    library(RColorBrewer)
    
    # Variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    vst_counts <- assay(vsd)
    
    # Stroke vs Diabetes results
    res_sd <- results(dds, contrast = c("condition","stroke","DIABETES_TYPE1"))
    
    # Remove NA
    res_sd <- res_sd[!is.na(res_sd$padj), ]
    
    # Order genes by significance
    res_ordered <- res_sd[order(res_sd$padj), ]
    
    # Select top 10 genes
    top10_genes <- rownames(res_ordered)[1:10]
    
    # Match metadata order
    colData_heatmap <- colData[colnames(vst_counts), , drop = FALSE]
    
    # Heatmap (All conditions)
    pheatmap(
      vst_counts[top10_genes, ],
      annotation_col = colData_heatmap,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      
      fontsize_row = 8,
      fontsize_col = 8,
      cellheight = 20,
      cellwidth  = 40,
      
      show_rownames = TRUE,
      show_colnames = TRUE,
      
      main = "Top 10 DEGs: Stroke vs Diabetes (All Samples)",
      
      color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(50)
    )
    
    levels(dds$condition)
    resultsNames(dds)
     dds$condition <- relevel(dds$condition, ref = "healthy.control")
     
   
     
##########################
## Expression of gene 
      library(ggplot2)
    
    gene <- ""
    
    df <- data.frame(
      expression = vst_counts[gene, ],
      condition = colData$condition
    )
    
    ggplot(df, aes(x=condition, y=expression, fill=condition)) +
      geom_violin(alpha=0.6) +
      geom_boxplot(width=0.2) +
      theme_minimal() +
      labs(title=paste("Expression of",gene))
    
###Top 10 genes ka combined plot (violin + boxplot) 

    # Stroke vs Control example
    res_sc <- results(dds, contrast = c("condition","stroke","Healthy_control"))
    
    # Remove NA
    res_sc <- res_sc[!is.na(res_sc$padj), ]
    
    # Sort by significance
    res_sc <- res_sc[order(res_sc$padj), ]
    
    # Top 10 genes
    top10_genes <- rownames(res_sc)[1:10]
##Expression data extract karo
    vsd_mat <- assay(vsd)
    
    # Subset matrix
    mat <- vsd_mat[top10_genes, ]
    
    # Convert to dataframe
    df <- as.data.frame(t(mat))
    df$condition <- colData(dds)$condition

##Long format (IMPORTANT)
    library(tidyr)
    
    df_long <- pivot_longer(
      df,
      cols = all_of(top10_genes),
      names_to = "gene",
      values_to = "expression"
    )
    
 ##Combined plot ##################
    library(ggplot2)
    
    ggplot(df_long, aes(x=condition, y=expression, fill=condition)) +
      
      geom_violin(trim=FALSE, alpha=0.6) +
      
      geom_boxplot(width=0.2, outlier.shape=NA, color="black") +
      
      facet_wrap(~gene, scales="free_y") +
      
      theme_minimal(base_size = 14) +
      
      labs(
        title = "Top 10 Differentially Expressed Genes",
        x = "Condition",
        y = "Expression"
      ) +
      
      theme(
        plot.title = element_text(hjust=0.5, face="bold"),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none"
      )
   
    ## improve colur 
    scale_fill_manual(values=c(
      "Healthy_control"="green",
      "DIABETES_TYPE1"="red",
      "stroke"="blue"
    ))
#####################################
    # Use SAME genes as heatmap
    top10_genes <- common_genes   # IMPORTANT ✔
    library(SummarizedExperiment)
    
    vsd_mat <- assay(vsd)
    
    # subset only top genes
    mat <- vsd_mat[top10_genes, ]
    df <- as.data.frame(t(mat))
    
    df$condition <- colData(dds)$condition
    
    # reshape manually (safe method)
    df_long <- data.frame()
    
    for(g in top10_genes){
      temp <- data.frame(
        gene = g,
        expression = df[[g]],
        condition = df$condition
      )
      df_long <- rbind(df_long, temp)
    }
    library(ggplot2)
    
    ggplot(df_long, aes(x=condition, y=expression, fill=condition)) +
      
      geom_violin(trim=FALSE, alpha=0.6) +
      
      geom_boxplot(width=0.15, outlier.shape=NA, color="black") +
      
      facet_wrap(~gene, scales="free", ncol=3) +
      
      theme_minimal(base_size = 14) +
      
      labs(
        title = "Top 10 Genes Expression (Violin + Boxplot)",
        x = "Condition",
        y = "Expression"
      ) +
      
      theme(
        plot.title = element_text(hjust=0.5, face="bold"),
        legend.position = "none"
      )
    scale_fill_manual(values=c(
      "Healthy_control"="green",
      "DIABETES_TYPE1"="red",
      "stroke"="blue"
    ))
############# 2 part 
    top10_genes <- rownames(res_sc)[1:10]
    
    # Split into 2 groups
    genes_1 <- top10_genes[1:5]
    genes_2 <- top10_genes[6:10]
#Common genes
    common1 <- intersect(genes_1, rownames(vsd_mat))
    common2 <- intersect(genes_2, rownames(vsd_mat))
# conbine plote 
    plot_top5 <- function(gene_list, title){
      
      library(tidyr)
      library(ggplot2)
      
      mat <- vsd_mat[gene_list, ]
      df <- as.data.frame(t(mat))
      df$condition <- colData(dds)$condition
      
      df_long <- pivot_longer(
        df,
        cols = all_of(gene_list),
        names_to = "gene",
        values_to = "expression"
      )
      
      ggplot(df_long, aes(x=condition, y=expression, fill=condition)) +
        
        geom_violin(trim=FALSE, alpha=0.6) +
        geom_boxplot(width=0.2, outlier.shape=NA) +
        
        facet_wrap(~gene, scales="free_y", ncol=2) +
        
        theme_minimal(base_size = 14) +
        
        labs(title = title, x="Condition", y="Expression") +
        
        scale_fill_manual(values=c(
          "Healthy_control"="green",
          "DIABETES_TYPE1"="red",
          "stroke"="blue"
        )) +
        
        theme(
          plot.title = element_text(hjust=0.5, face="bold"),
          axis.text.x = element_text(angle=45, hjust=1),
          legend.position="none"
        )
    }
#clean plots
    plot_top5(common1, "Top 5 Genes (Set 1)")
    plot_top5(common2, "Top 5 Genes (Set 2)")
##########################################################    
## gene are clear and clean### 
    head(top10_genes)
    head(rownames(vsd_mat))
    common_genes <- intersect(top10_genes, rownames(vsd_mat))
    
    length(common_genes)   # check kitne match hue
    # used common gene 
    mat <- vsd_mat[common_genes, ]
    
    df <- as.data.frame(t(mat))
    df$condition <- colData(dds)$condition
    
    ##
    library(tidyr)
    
    df_long <- pivot_longer(
      df,
      cols = all_of(common_genes),
      names_to = "gene",
      values_to = "expression"
    )
    ## 
    common_genes = intersect(top10_genes, rownames(vsd_mat))
   
    
    
# new graph scatter + regression plot #######################################
    df <- data.frame(
      gene1 = vst_counts["HSD11B1",],
      gene2 = vst_counts["COL17A1",]
    )
    
    ggplot(df,aes(x=gene1,y=gene2)) +
      geom_point() +
      geom_smooth(method="lm") +
      theme_minimal()

##scatter + regression plot
    res_sc <- results(dds, contrast = c("condition","stroke","Healthy_control"))
    
    res_sc <- res_sc[!is.na(res_sc$padj), ]
    res_sc <- res_sc[order(res_sc$padj), ]
    
    top10_genes <- rownames(res_sc)[1:10]    
##created all 
    library(tidyr)
    library(dplyr)
    
    # Make gene pairs
    pairs <- combn(common_genes, 2, simplify = FALSE)
    
    # Take only first 10 pairs (clean figure)
    pairs <- pairs[1:10]
#bulid data fram plot 
    plot_list <- lapply(pairs, function(p){
      data.frame(
        x = df[[p[1]]],
        y = df[[p[2]]],
        pair = paste(p[1], "vs", p[2])
      )
    })
    
    plot_df <- do.call(rbind, plot_list)
 ## final plot 
    library(ggplot2)
    
    ggplot(plot_df, aes(x=x, y=y)) +
      
      geom_point(color="black", size=2) +
      
      geom_smooth(method="lm", color="blue", fill="grey", alpha=0.3) +
      
      facet_wrap(~pair, scales="free", ncol=3) +
      
      theme_minimal(base_size = 14) +
      
      labs(
        title = "scatter+regression plot of Top 10 Genes",
        x = "Gene Expression",
        y = "Gene Expression"
      ) +
      
      theme(
        plot.title = element_text(hjust=0.5, face="bold")
      )
  
########### ander condition c for scatter + regression plot 
    library(ggplot2)
    
    # Loop for first 5 pairs (avoid too many plots)
    for(i in 1:5){
      
      gene_x <- common_genes[i]
      gene_y <- common_genes[i+1]
      
      p <- ggplot(df, aes_string(x=gene_x, y=gene_y)) +
        geom_point(color="black", size=2) +
        geom_smooth(method="lm", color="blue", fill="grey", alpha=0.3) +
        theme_minimal() +
        labs(
          title = paste("Correlation:", gene_x, "vs", gene_y),
          x = gene_x,
          y = gene_y
        )
      
      print(p)
    }
######3scatter + regression plot    
    library(tidyr)
    
    # Example: only first 3 pairs
    pairs_df <- data.frame(
      x = df[[common_genes[1]]],
      y = df[[common_genes[2]]],
      pair = paste(common_genes[1], "vs", common_genes[2])
    )
    
    ggplot(pairs_df, aes(x=x, y=y)) +
      geom_point() +
      geom_smooth(method="lm", color="blue") +
      facet_wrap(~pair) +
      theme_minimal()
    
#####################################################################################################################
  
    res_sc <- results(dds, contrast = c("condition","stroke","Healthy_control"))
    
    res_sc <- res_sc[!is.na(res_sc$padj), ]
    res_sc <- res_sc[order(res_sc$padj), ]
    
    top10_genes <- rownames(res_sc)[1:10]
    
  #Expression matrix
    vsd_mat <- assay(vsd)
    
    # Match genes
    common_genes <- intersect(top10_genes, rownames(vsd_mat))
    
    mat <- vsd_mat[common_genes, ]
 # transpose
    df <- as.data.frame(t(mat))

 ## correction plot 
    install.packages("GGally")
    library(GGally)
    
    ggpairs(
      df,
      title = "Correlation Plot of Top 10 Genes"
    )
    
  ## bonus 
    ggpairs(
      df,
      upper = list(continuous = wrap("cor", size=4)),  # show correlation value
      lower = list(continuous = wrap("points", alpha=0.6)),
      diag = list(continuous = wrap("densityDiag")),
      title = "Top 10 Gene Correlation"
    )
## altrenative 
    library(pheatmap)
    
    cor_mat <- cor(df)
    
    pheatmap(
      cor_mat,
      display_numbers = TRUE,
      main = "Correlation Heatmap (Top 10 Genes)"
    )
    
    ##############################################################################################################
############### 4.protien protien interation ppi network ##################################### 
    
    install.packages(c("igraph","ggraph","tidygraph","ggplot2"))
    
    if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
    
    BiocManager::install("STRINGdb")

# load labrary
library(STRINGdb)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
  
## selected de gene 
deg_genes <- rownames(res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ])

# Initialize STRING database 
string_db <- STRINGdb$new(
  version="11.5",
  species=9606,
  score_threshold=400,
  input_directory=""
)

ppi_genes <- rownames(sig_genes_ppi)

# Network stable rakhne ke liye limit
ppi_genes <- head(ppi_genes, 50)
##Map genes to STRING IDs
gene_df <- data.frame(gene=deg_genes)
mapped <- string_db$map(
  gene_df,
  "gene",
  removeUnmappedRows = TRUE
)
head(mapped)

ppi <- string_db$get_interactions(mapped$STRING_id)
## Get protein interactions
library(STRINGdb)
library(igraph)

# Get interactions
# Get interactions from STRING
ppi <- string_db$get_interactions(mapped$STRING_id)

# Keep required columns
ppi_edges <- ppi[, c("from", "to", "combined_score")]

# Remove NA rows
ppi_edges <- ppi_edges[complete.cases(ppi_edges), ]

# Remove duplicate interactions
ppi_edges <- unique(ppi_edges)

# Keep high-confidence interactions
ppi_edges <- subset(ppi_edges, combined_score > 300)

head(mapped)
mapped_ids <- mapped$STRING_id

mapped_ids <- mapped_ids[!is.na(mapped_ids)]

mapped_ids <- unique(mapped_ids)
length(mapped_ids)
head(mapped_ids)

ppi <- tryCatch(
  string_db$get_interactions(mapped_ids),
  error = function(e) NULL
)
ppi <- string_db$get_subnetwork(mapped_ids)
ppi <- as_data_frame(ppi, what="edges")




interactions <- string_db$interactions
# Check first rows
head(ppi_edges)

# Create graph
g <- graph_from_data_frame(
  ppi_edges[, c("from","to")],
  directed = FALSE
)

# Plot network
plot(g, vertex.size=5, vertex.label=NA)

### Convert to network graph
g <- graph_from_data_frame(
  ppi[, c("from","to")], 
  directed = FALSE
)

#####

# Define thresholds
padj_threshold <- 0.05
lfc_threshold  <- 1

# Select significant genes
sig_genes_ppi <- subset(
  res,
  padj < padj_threshold & abs(log2FoldChange) > lfc_threshold
)

# Count genes
num_sig <- nrow(sig_genes_ppi)

cat("Number of significant genes (padj <", padj_threshold,
    "& |LFC| >", lfc_threshold, "):", num_sig, "\n")


# If too few genes, relax threshold
if (num_sig < 50) {
  
  warning("Too few genes detected. Using relaxed threshold.")
  
  sig_genes_ppi <- subset(
    res,
    padj < 0.1 & abs(log2FoldChange) > 0.5
  )
  
}

# Sort by significance
sig_genes_ppi <- sig_genes_ppi[order(sig_genes_ppi$padj), ]

# Save gene list
ppi_gene_list <- data.frame(
  GeneID = rownames(sig_genes_ppi),
  log2FoldChange = sig_genes_ppi$log2FoldChange,
  padj = sig_genes_ppi$padj
)

write.csv(ppi_gene_list, "PPI_DEGs.csv", row.names = FALSE)

ppi_genes <- rownames(sig_genes_ppi)

########## PPI network with hub genes ####
deg <- degree(g)

V(g)$size <- deg

ggraph(g, layout="fr") +
  geom_edge_link(alpha=0.3) +
  geom_node_point(aes(size=size), color="red") +
  geom_node_text(aes(label=name), repel=TRUE, size=3) +
  theme_void() +
  ggtitle("PPI Network with Hub Genes") 


#################### KEGG Pathway Enrichment + Plot ######
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

### list gene prepar ##
# Use your significant genes
genes <- rownames(sig_genes_ppi)

# Convert gene symbols → Entrez IDs
gene_ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

### KEGG enrichment analysis

kegg_result <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

dotplot(
  kegg_result,
  showCategory = 10,
  title = "KEGG Pathway Enrichment"
)
head(kegg_result)

kegg_df <- as.data.frame(kegg_result)

ggplot(kegg_df[1:10, ], aes(
  x = GeneRatio,
  y = reorder(Description, GeneRatio),
  size = Count,
  color = p.adjust
)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(
    title = "KEGG Pathway Enrichment",
    x = "Gene Ratio",
    y = "Pathway"
  )
as.data.frame(kegg_result)

##
sig_genes_ppi <- subset(
  res,
  padj < 0.2 & abs(log2FoldChange) > 0.3
)

###
genes <- rownames(sig_genes_ppi)

gene_ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_ids <- na.omit(gene_ids)

# Check
head(gene_ids)
nrow(gene_ids)

kegg_result <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.5   # IMPORTANT (relax)
)

## Check result 
kegg_df <- as.data.frame(kegg_result)

nrow(kegg_df)

ggplot(kegg_df[1:min(10, nrow(kegg_df)), ], aes(
  x = GeneRatio,
  y = reorder(Description, GeneRatio),
  size = Count,
  color = p.adjust
)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(
    title = "KEGG Pathway Enrichment",
    x = "Gene Ratio",
    y = "Pathway"
  )


######### KEGG Bar Plot ## 
kegg_df <- as.data.frame(kegg_result)

# Check
nrow(kegg_df)

## Top 10 pathways select 
top_kegg <- kegg_df[1:min(10, nrow(kegg_df)), ]

## Bar plot (FINAL)

library(ggplot2)

ggplot(top_kegg, aes(
  x = Count,
  y = reorder(Description, Count),
  fill = p.adjust
)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(
    low = "red",
    high = "blue",
    name = "p.adjust"
  ) +
  theme_minimal() +
  labs(
    title = "Top 10 KEGG Pathways",
    x = "Count",
    y = ""
  )

head(kegg_res)
##### Reactome Pathway Analysis (Full Pipeline) ### 

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install(c("ReactomePA","clusterProfiler","org.Hs.eg.db","enrichplot"))

library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

 ####### Gene list prepare
# Use DE genes
genes <- rownames(sig_genes_ppi)

# Convert to ENTREZ
gene_ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
 
gene_ids <- na.omit(gene_ids)

## Reactome enrichment 
reactome_result <- enrichPathway(
  gene = gene_ids$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.1
)
## Bar Plot (Top pathways)
react_df <- as.data.frame(reactome_result)

top_react <- react_df[1:min(10, nrow(react_df)), ]

ggplot(top_react, aes(
  x = Count,
  y = reorder(Description, Count)
)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  labs(
    title = "Top Enriched Reactome Pathways",
    x = "Gene Count",
    y = ""
  )

## pathway for volcano plote ## 
ggplot(react_df, aes(
  x = Count,
  y = -log10(p.adjust)
)) +
  geom_point(color="red", size=2) +
  theme_minimal() +
  labs(
    title="Reactome Volcano Plot",
    x="Gene Count",
    y="-log10(p.adjust)"
  )
reactome_result

react_df <- as.data.frame(reactome_result)
nrow(react_df)
head(react_df)


################33
genes <- rownames(sig_genes_ppi)

gene_ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_ids <- na.omit(gene_ids)

reactome_result <- enrichPathway(
  gene = gene_ids$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.5   # keep relaxed
)

react_df <- as.data.frame(reactome_result)

# Check
nrow(react_df)
head(react_df)

gsea_res <- gsePathway(
  geneList = gene_list,
  organism = "human",
  pvalueCutoff = 0.5
)

gseaplot2(
  gsea_res,
  geneSetID = 1,   # top pathway
  title = gsea_res@result$Description[1]
)

gseaplot2(
  gsea_res,
  geneSetID = "REACTOME_INTERFERON_SIGNALING"
)
gene_list <- res$log2FoldChange

names(gene_list) <- rownames(res)

# Remove NA
gene_list <- gene_list[!is.na(gene_list)]

# Sort
gene_list <- sort(gene_list, decreasing = TRUE)

gene_df <- bitr(
  names(gene_list),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Match genes
gene_list <- gene_list[gene_df$SYMBOL]

names(gene_list) <- gene_df$ENTREZID

gsea_res <- gsePathway(
  geneList = gene_list,
  organism = "human",
  pvalueCutoff = 0.5
)

gseaplot2(
  gsea_res,
  geneSetID = 1,
  title = gsea_res@result$Description[1]
)
######################################################################

gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)

gene_list <- gene_list[!is.na(gene_list)]

gene_list <- sort(gene_list, decreasing = TRUE)

# Check
length(gene_list)

library(org.Hs.eg.db)
library(clusterProfiler)

gene_df <- bitr(
  names(gene_list),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Keep only mapped genes
gene_list <- gene_list[gene_df$SYMBOL]
names(gene_list) <- gene_df$ENTREZID

library(ReactomePA)

gsea_res <- gsePathway(
  geneList = gene_list,
  organism = "human",
  pvalueCutoff = 1   # IMPORTANT (max relax)
)

nrow(gsea_res@result)
head(gsea_res@result)

library(enrichplot)

gseaplot2(
  gsea_res,
  geneSetID = 1
)

length(gene_list)
head(gene_list)
head(gsea_res@result)

##########################################################################
 # 3 comparisons between pathway analysis (GSEA Reactome) 

levels(dds$condition)

library(DESeq2)
# Stroke vs Control
res_stroke <- results(dds, contrast=c("condition","stroke","Healthy_control"))

# Diabetes vs Control
res_diabetes <- results(dds, contrast=c("condition","DIABETES_TYPE1","Healthy_control"))

# Stroke vs Diabetes
res_sd <- results(dds, contrast=c("condition","stroke","DIABETES_TYPE1"))

head(res_stroke)
head(res_diabetes)
head(res_stroke_vs_diabetes)

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
run_gsea <- function(res_obj) {
  
  # Gene list
  gene_list <- res_obj$log2FoldChange
  names(gene_list) <- rownames(res_obj)
  
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Convert SYMBOL → ENTREZ
  gene_df <- bitr(
    names(gene_list),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # REMOVE DUPLICATES (IMPORTANT FIX)
  gene_df <- gene_df[!duplicated(gene_df$ENTREZID), ]
  
  # Match gene list
  gene_list <- gene_list[gene_df$SYMBOL]
  names(gene_list) <- gene_df$ENTREZID
  
  # Sort again
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run GSEA
  gsea_res <- gsePathway(
    geneList = gene_list,
    organism = "human",
    pvalueCutoff = 1
  )
  
  return(gsea_res)
}

gsea_stroke <- run_gsea(res_stroke)
gsea_diabetes <- run_gsea(res_diabetes)
gsea_sd <- run_gsea(res_sd)

gseaplot2(
  gsea_stroke,
  geneSetID = 1,
  title = "Stroke vs Control"
)

 

gseaplot2(
  gsea_diabetes,
  geneSetID = 1,
  title = "Stroke vs Diabetes",
  subplots = 1:3,
  pvalue_table = TRUE
)
# Top pathway genes
core_genes <- gsea_diabetes@result$core_enrichment[1]

core_genes
library(enrichplot)

gseaplot2(
  gsea_diabetes,
  geneSetID = 1,
  title = "Diabetes vs Control",
  ES_geom = "line"
)


gseaplot2(
  gsea_diabetes,
  geneSetID = 1,
  title = "Diabetes vs Control",
  color = "green"
)
gseaplot2(
  gsea_sd,
  geneSetID = 1,
  title = "Stroke vs Diabetes"
)
#####################################REACTOME GSEA + BARPLOT 
#reactome bar plot 
head(react_dc@result)
nrow(react_dc@result)

df <- as.data.frame(react_dc@result)

sum(is.na(df$NES))

dotplot(react_dc, showCategory=10)
plot_reactome <- function(res, title){
  
  library(ggplot2)
  
  df <- as.data.frame(res@result)
  
  # Remove NA
  df <- df[!is.na(df$NES) & is.finite(df$NES), ]
  
  # If empty → stop
  if(nrow(df) == 0){
    print("No valid pathways found")
    return(NULL)
  }
  
  # Top pathways
  df <- df[order(abs(df$NES), decreasing=TRUE), ]
  df <- df[1:min(10, nrow(df)), ]
  
  df$Description <- factor(df$Description, levels = df$Description[order(df$NES)])
  
  ggplot(df, aes(x=Description, y=NES)) +
    geom_col(fill="#4E79A7") +
    geom_vline(xintercept=0) +
    coord_flip() +
    theme_minimal() +
    ggtitle(title)
}
plot_reactome(react_dc, "Diabetes vs Control")
plot_reactome(react_sd, "Stroke vs Diabetes")
plot_reactome(react_sc, "stroke vs control")

head(react_dc@result)
summary(react_dc@result$NES)


## #################### Reactome Volcano Plot (3 Comparisons) 
library(ggplot2)
library(ggrepel)

reactome_volcano <- function(gsea_obj, title_name) {
  
  df <- as.data.frame(gsea_obj@result)
  
  if (nrow(df) == 0) {
    print(paste("No pathways for", title_name))
    return(NULL)
  }
  
  # Define significance
  df$significance <- ifelse(df$p.adjust < 0.05, "Significant", "Not Significant")
  
  # Select top pathways for labeling
  top_pathways <- head(df[order(df$p.adjust), ], 10)
  
  ggplot(df, aes(
    x = NES,
    y = -log10(p.adjust),
    color = significance
  )) +
    geom_point(size = 2) +
    
    #  Add pathway names
    geom_text_repel(
      data = top_pathways,
      aes(label = Description),
      size = 3,
      max.overlaps = 20
    ) +
    
    scale_color_manual(values = c("blue","red")) +
    theme_minimal() +
    labs(
      title = title_name,
      x = "NES",
      y = "-log10(p.adjust)"
    )
}
p1 <- reactome_volcano(gsea_stroke, "Stroke vs Control")
p2 <- reactome_volcano(gsea_diabetes, "Diabetes vs Control")
p3 <- reactome_volcano(gsea_sd, "Stroke vs Diabetes")

install.packages("patchwork")
library(patchwork)
p2
p1 + p2 + p3


################################################################################
 ## GO Pathway Analysis (ALL 3 COMPARISONS)
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler","org.Hs.eg.db"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)

##Function (GO enrichment)

run_go <- function(res_obj, ont_type="BP") {
  
  # Gene selection (relaxed for better result)
  genes <- rownames(
    subset(res_obj, padj < 0.1 & abs(log2FoldChange) > 0.5)
  )
  
  # Convert SYMBOL → ENTREZ
  gene_ids <- bitr(
    genes,
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )
  
  gene_ids <- na.omit(gene_ids)
  
  # GO enrichment
  go_res <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = ont_type,   # BP / MF / CC
    pvalueCutoff = 0.1,
    readable = TRUE
  )
  
  return(go_res)
}

####
# BP (Biological Process)

go_stroke <- run_go(res_stroke, "BP")
go_diabetes <- run_go(res_diabetes, "BP")
go_sd <- run_go(res_sd, "BP")

plot_go <- function(go_obj, title_name) {
  
  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) {
    return(ggplot() + ggtitle(paste("No GO terms:", title_name)))
  }
  
  dotplot(
    go_obj,
    showCategory = 10,
    title = title_name
  )
}

p1 <- plot_go(go_stroke, "Stroke vs Control")
p2 <- plot_go(go_diabetes, "Diabetes vs Control")
p3 <- plot_go(go_sd, "Stroke vs Diabetes")

p1


run_go <- function(res_obj, ont_type="BP") {
  
  # RELAXED FILTER (IMPORTANT)
  genes <- rownames(
    subset(res_obj, padj < 0.2 & abs(log2FoldChange) > 0.3)
  )
  
  # If still low, take top genes
  if (length(genes) < 20) {
    genes <- rownames(res_obj)[order(res_obj$padj)][1:100]
  }
  
  # Convert SYMBOL → ENTREZ
  gene_ids <- bitr(
    genes,
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )
  
  gene_ids <- na.omit(gene_ids)
  
  # GO enrichment
  go_res <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = ont_type,
    pvalueCutoff = 0.5, 
    readable = TRUE
  )
  
  return(go_res)
}

go_stroke <- run_go(res_stroke, "BP")
go_diabetes <- run_go(res_diabetes, "BP")
go_sd <- run_go(res_sd, "BP") 
nrow(as.data.frame(go_stroke))

p1 <- dotplot(go_stroke, showCategory=10, title="Stroke vs Control")
p2 <- dotplot(go_diabetes, showCategory=10, title="Diabetes vs Control")
p3 <- dotplot(go_sd, showCategory=10, title="Stroke vs Diabetes")

library(patchwork)
p1 + p2 + p3

###############################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(patchwork)
library(ggrepel)
make_GO_panel <- function(res_obj, title_name) {
  
  # 🔹 Gene list (GSEA)
  gene_list <- res_obj$log2FoldChange
  names(gene_list) <- rownames(res_obj)
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Convert IDs
  gene_df <- bitr(
    names(gene_list),
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )
  
  gene_list <- gene_list[gene_df$SYMBOL]
  names(gene_list) <- gene_df$ENTREZID
  
  # 🔹 GSEA GO
  gsea_go <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 1
  )
  
  go_df <- as.data.frame(gsea_go@result)
  
  # empty 
  if (nrow(go_df) == 0) {
    return(ggplot() + ggtitle(paste("No GO:", title_name)))
  }
  
  # 🔹 Volcano plot
  volcano <- ggplot(go_df, aes(
    x = NES,
    y = -log10(p.adjust)
  )) +
    geom_point(color="red", size=2) +
    theme_minimal() +
    labs(title=paste("Volcano:", title_name),
         x="NES",
         y="-log10(p.adjust)")
  
  # 🔹 Bar plot
  top_go <- go_df[1:10, ]
  
  bar <- ggplot(top_go, aes(
    x = NES,
    y = reorder(Description, NES)
  )) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_minimal() +
    labs(title="Top GO Pathways",
         x="NES",
         y="")
  
  # 🔹 GSEA curve
  gsea_plot <- gseaplot2(
    gsea_go,
    geneSetID = 1,
    title = go_df$Description[1]
  )
  
  # Combine 3 plots
  volcano | bar | gsea_plot
}

## 
p_stroke <- make_GO_panel(res_stroke, "Stroke vs Control")
p_diabetes <- make_GO_panel(res_diabetes, "Diabetes vs Control")
p_sd <- make_GO_panel(res_sd, "Stroke vs Diabetes")

p_stroke
p_diabetes
p_sd

######## 
plot_go_volcano <- function(gsea_obj, title_name) {
  
  df <- as.data.frame(gsea_obj@result)
  
  ggplot(df, aes(NES, -log10(p.adjust))) +
    geom_point(color="red") +
    theme_minimal() +
    labs(title=title_name,
         x="NES",
         y="-log10(p.adjust)")
}

# Run
plot_go_volcano(gsea_stroke, "Stroke vs Control")
plot_go_volcano(gsea_diabetes, "Diabetes vs Control")
plot_go_volcano(gsea_sd, "Stroke vs Diabetes")

###########3
library(ggplot2)
library(ggrepel)

plot_go_volcano <- function(gsea_obj, title_name) {
  
  df <- as.data.frame(gsea_obj@result)
  
  # 🔍 Check data
  print(head(df))
  print(nrow(df))
  
  if (nrow(df) == 0) {
    stop(_)
  }
  
  # Remove NA
  df <- df[!is.na(df$p.adjust) & !is.na(df$NES), ]
  
  # Significance
  df$significance <- ifelse(df$p.adjust < 0.05, "Significant", "Not Significant")
  
  # Force top pathways (guaranteed)
  top_pathways <- df[order(df$p.adjust), ][1:min(10, nrow(df)), ]
  
  ggplot(df, aes(
    x = NES,
    y = -log10(p.adjust)
  )) +
    geom_point(color="red", size=2) +
    
    # LABELS FORCE SHOW
    geom_text_repel(
      data = top_pathways,
      aes(label = Description),
      size = 4,
      color = "black",
      box.padding = 0.5,
      segment.color = "black",
      max.overlaps = Inf   # IMPORTANT
    ) +
    
    theme_minimal() +
    labs(
      title = title_name,
      x = "NES",
      y = "-log10(p.adjust)"
    )
}

top_pathways <- df[df$p.adjust < 0.05, ][1:5, ]
if ("Description" %in% colnames(df)) {
  df$Description <- substr(as.character(df$Description), 1, 40)
}
df$ShortName <- substr(as.character(df$Description), 1, 40)
aes(label = ShortName)

df$ShortName <- substr(as.character(df$Description), 1, 40)

geom_text_repel(
  data = top_pathways,
  aes(label = ShortName),
  size = 3.5
)

plot_go_volcano(gsea_stroke, "Stroke vs Control")
plot_go_volcano(gsea_diabetes, "Diabetes vs Control")
plot_go_volcano(gsea_sd, "Stroke vs Diabetes")

## 
library(enrichplot)

plot_go_gsea <- function(gsea_obj, title_name) {
  
  # Check data
  if (nrow(gsea_obj@result) == 0) {
    print(paste("No pathway:", title_name))
    return(NULL)
  }
  
  gseaplot2(
    gsea_obj,
    geneSetID = 1,   # top pathway
    title = title_name,
    color = "green",
    pvalue_table = FALSE
  )
}

plot_go_gsea(gsea_stroke, "Stroke vs Control")
plot_go_gsea(gsea_sd, "Stroke vs Diabetes")
plot_go_gsea(gsea_diabetes, "Diabetes vs Control")

#Stroke vs Control
top1 <- gsea_stroke@result$Description[1]

gseaplot2(
  gsea_stroke,
  geneSetID = 1,
  title = paste("GO Pathway Enrichment\n", top1),
  color = "green"
)

#Diabetes vs Control
top2 <- gsea_diabetes@result$Description[1]

gseaplot2(
  gsea_diabetes,
  geneSetID = 1,
  title = paste("GO Pathway Enrichment\n", top2),
  color = "green"
)

#Stroke vs Diabetes
top3 <- gsea_sd@result$Description[1]

gseaplot2(
  gsea_sd,
  geneSetID = 1,
  title = paste("GO Pathway Enrichment\n", top3),
  color = "green"
)
############################# og barplot 
res_sc <- results(dds, contrast = c("condition","stroke","Healthy_control"))
res_dc <- results(dds, contrast = c("condition","DIABETES_TYPE1","Healthy_control"))
res_sd <- results(dds, contrast = c("condition","stroke","DIABETES_TYPE1"))

library(clusterProfiler)
library(org.Hs.eg.db)

prepare_gene_list <- function(res){
  
  res <- res[!is.na(res$padj), ]
  
  gene_list <- res$log2FoldChange
  names(gene_list) <- rownames(res)
  
  gene_df <- bitr(names(gene_list),
                  fromType="SYMBOL",
                  toType="ENTREZID",
                  OrgDb=org.Hs.eg.db)
  
  gene_list <- gene_list[gene_df$SYMBOL]
  names(gene_list) <- gene_df$ENTREZID
  
  sort(gene_list, decreasing=TRUE)
}

gene_sc <- prepare_gene_list(res_sc)
gene_dc <- prepare_gene_list(res_ds)
gene_sd <- prepare_gene_list(res_sd)

kegg_sc <- gseKEGG(geneList = gene_sc, organism = "hsa")
kegg_dc <- gseKEGG(geneList = gene_dc, organism = "hsa")
kegg_sd <- gseKEGG(geneList = gene_sd, organism = "hsa")


library(ggplot2)

plot_kegg <- function(kegg_res, title){

  df <- as.data.frame(kegg_res@result)
  df <- df[order(df$NES, decreasing = TRUE), ][1:10, ]
  
  ggplot(df, aes(x=reorder(Description, NES), y=NES)) +
    geom_bar(stat="identity", fill="steelblue") +
    coord_flip() +
    theme_minimal() +
    ggtitle(title) +
    xlab("") +
    ylab("NES")
}

plot_kegg(kegg_sc, "Stroke vs Control")
plot_kegg(kegg_dc, "Diabetes vs Control")
plot_kegg(kegg_sd, "Stroke vs Diabetes")




#######################################################################################################
#### clean og pathway for only this condioton clean and mange data and row remove for stroke vs control 
plot_kegg <- function(kegg_res, title){
  
  df <- as.data.frame(kegg_res@result)
  
  #  CLEAN DATA
  df <- df[complete.cases(df), ]
  df <- df[is.finite(df$NES), ]
  df <- df[!is.na(df$Description) & df$Description != "", ]
  
  # Sort
  df <- df[order(df$NES, decreasing = TRUE), ]
  
  # Take top safely
  top_n <- min(10, nrow(df))
  df <- df[1:top_n, ]
  
  # Factor fix
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  library(ggplot2)
  
  ggplot(df, aes(x=Description, y=NES)) +
    geom_col(fill="#4E79A7", na.rm = TRUE) +   #  THIS FIXES WARNING
    coord_flip() +
    theme_minimal() +
    ggtitle(title) +
    xlab("") +
    ylab("NES")
}

plot_kegg(kegg_sc, "Stroke vs Control")
head(kegg_sc@result)

##
length(res_sc$log2FoldChange)
gene_list <- res_sc$log2FoldChange
names(gene_list) <- rownames(res_sc)

# Remove NA
gene_list <- gene_list[!is.na(gene_list)]
library(clusterProfiler)
library(org.Hs.eg.db)

gene_df <- bitr(names(gene_list),
                fromType="SYMBOL",
                toType="ENTREZID",
                OrgDb=org.Hs.eg.db)

# Check mapping
head(gene_df)
nrow(gene_df)
gene_list <- gene_list[gene_df$SYMBOL]
names(gene_list) <- gene_df$ENTREZID

gene_list <- sort(gene_list, decreasing=TRUE)
kegg_sc <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1   # IMPORTANT (no filtering)
)
nrow(kegg_sc@result)
head(kegg_sc@result)
dotplot(kegg_sc, showCategory=10)
kegg_sc <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 5,
  pvalueCutoff = 1
)

kegg_df <- as.data.frame(kegg_sc@result)


###############################################################################


##HALLMARK PATHWAY ############################################

# libreray 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr)
install.packages("msigdbr")
library(ggrepel)

## hallmark gene set 
hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

## GSEA FUNCTION GENE 
run_hallmark <- function(res_obj) {
  
  gene_list <- res_obj$log2FoldChange
  names(gene_list) <- rownames(res_obj)
  
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_res <- GSEA(
    geneList = gene_list,
    TERM2GENE = hallmark[, c("gs_name","gene_symbol")],
    pvalueCutoff = 1
  )
  
  return(gsea_res)
}

hall_stroke <- run_hallmark(res_stroke)
hall_diabetes <- run_hallmark(res_diabetes)
hall_sd <- run_hallmark(res_sd)

## hallmark volcano plot 
plot_hall_volcano <- function(gsea_obj, title_name) {
  
  df <- as.data.frame(gsea_obj@result)
  
  top <- df[order(df$p.adjust), ][1:6, ]
  
  ggplot(df, aes(NES, -log10(p.adjust))) +
    geom_point(color="red", size=2) +
    
    geom_text_repel(
      data = top,
      aes(label = Description),
      size = 3
    ) +
    
    theme_minimal() +
    labs(title = paste("Volcano:", title_name),
         x="NES", y="-log10(p.adjust)")
}

plot_hall_volcano(hall_stroke, "Stroke vs Control")
plot_hall_volcano(hall_diabetes, "Diabetes vs Control")
plot_hall_volcano(hall_sd, "Stroke vs Diabetes")

## 2. HALLMARK BAR PLOT
plot_hall_bar <- function(gsea_obj, title_name) {
  
  df <- as.data.frame(gsea_obj@result)
  top <- df[1:10, ]
  
  ggplot(top, aes(NES, reorder(Description, NES))) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_minimal() +
    labs(title="Top Hallmark Pathways",
         x="NES", y="")
}
plot_hall_bar(hall_stroke, "Stroke vs Control")
plot_hall_bar(hall_diabetes, "Diabetes vs Control")
plot_hall_bar(hall_sd, "Stroke vs Diabetes")

### 3. HALLMARK GSEA PLOT
plot_hall_gsea <- function(gsea_obj, title_name) {
  
  top_path <- gsea_obj@result$Description[1]
  
  gseaplot2(
    gsea_obj,
    geneSetID = 1,
    title = paste("Hallmark Pathway Enrichment\n", top_path),
    color = "green"
  )
}

plot_hall_gsea(hall_stroke, "Stroke vs Control")
plot_hall_gsea(hall_diabetes, "Diabetes vs Control")
plot_hall_gsea(hall_sd, "Stroke vs Diabetes")

####################################################################################


4.#### KEGG PATHWAY

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggrepel)

##2️⃣ KEGG GSEA function 
run_kegg <- function(res_obj) {
  
  gene_list <- res_obj$log2FoldChange
  names(gene_list) <- rownames(res_obj)
  
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Convert to ENTREZ
  gene_df <- bitr(
    names(gene_list),
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )
  
  gene_list <- gene_list[gene_df$SYMBOL]
  names(gene_list) <- gene_df$ENTREZID
  
  # KEGG GSEA
  gsea_kegg <- gseKEGG(
    geneList = gene_list,
    organism = "hsa",
    pvalueCutoff = 1
  )
  
  return(gsea_kegg)
}

## run all 
kegg_stroke <- run_kegg(res_stroke)
kegg_diabetes <- run_kegg(res_diabetes)
kegg_sd <- run_kegg(res_sd)

### 4️⃣ 1. KEGG VOLCANO

plot_kegg_volcano <- function(gsea_obj, title_name) {
  
  df <- as.data.frame(gsea_obj@result)
  
  df$significance <- ifelse(df$p.adjust < 0.05, "TRUE", "FALSE")
  
  top <- df[order(df$p.adjust), ][1:6, ]
  
  ggplot(df, aes(NES, -log10(p.adjust))) +
    geom_point(aes(color=significance), size=2) +
    
    geom_text_repel(
      data = top,
      aes(label = Description),
      size = 3
    ) +
    
    scale_color_manual(values=c("gray","red")) +
    
    theme_minimal() +
    labs(title="Volcano Plot of KEGG Pathway Enrichment",
         x="NES",
         y="-log10 Adjusted P-value")
}

## plot print
plot_kegg_volcano(kegg_stroke, "Stroke vs Control")
plot_kegg_volcano(kegg_diabetes, "Diabetes vs Control")
plot_kegg_volcano(kegg_sd, "Stroke vs Diabetes")

## 2. KEGG BAR PLOT
plot_kegg_bar <- function(gsea_obj) {
  
  df <- as.data.frame(gsea_obj@result)
  top <- df[1:10, ]
  
  ggplot(top, aes(NES, reorder(Description, NES))) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_minimal() +
    labs(title="Top Enriched KEGG Pathways",
         x="NES",
         y="")
}
plot_kegg_bar(kegg_stroke)
plot_kegg_bar(kegg_diabetes)
plot_kegg_bar(kegg_sd)

### 3. KEGG GSEA PLOT
plot_kegg_gsea <- function(gsea_obj) {
  
  top_path <- gsea_obj@result$Description[1]
  
  gseaplot2(
    gsea_obj,
    geneSetID = 1,
    title = paste("KEGG Pathway Enrichment\n", top_path),
    color = "green"
  )
}

plot_kegg_gsea(kegg_stroke)
plot_kegg_gsea(kegg_diabetes)
plot_kegg_gsea(kegg_sd)


##################################################################################

## PPi 
library(STRINGdb)
library(igraph)
library(ggraph)
library(tidygraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)
### 1.DEG genes select
deg_genes <- rownames(
  subset(res_stroke, padj < 0.05 & abs(log2FoldChange) > 1)
)

### 2. PPI NETWORK 
string_db <- STRINGdb$new(
  version="11.5",
  species=9606,
  score_threshold=400
)

gene_df <- data.frame(gene=deg_genes)

mapped <- string_db$map(gene_df, "gene", removeUnmappedRows=TRUE)

# interactions
ppi <- string_db$get_interactions(mapped$STRING_id)

ppi <- ppi[complete.cases(ppi), ]
ppi <- subset(ppi, combined_score > 400)

g <- graph_from_data_frame(ppi[,c("from","to")], directed=FALSE)

# plot
p0 <- ggraph(g, layout="fr") +
  geom_edge_link(alpha=0.3) +
  geom_node_point(aes(color=degree(g)), size=3) +
  geom_node_text(aes(label=name), size=2, repel=TRUE) +
  theme_void() +
  ggtitle("Protein-Protein Interaction Network")
p3

## KEGG BAR PLOT 

p2
p1


go <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.2
)

head(as.data.frame(go))


########################## summary ####
## top gene upregulated and downregulated 
# Define thresholds
 padj_threshold <- 0.05
 lfc_threshold <- 1  
 
 # Identify significant genes
 sig_genes <- subset(res, padj < padj_threshold & abs(log2FoldChange) > lfc_threshold)
  
    # Identify upregulated and downregulated genes
    up_genes <- subset(sig_genes, log2FoldChange > lfc_threshold)
   down_genes <- subset(sig_genes, log2FoldChange < -lfc_threshold)
  
    # Get top 10 genes from each category
    top_10_sig <- head(sig_genes[order(sig_genes$padj), ], 10)
   top_10_up <- head(up_genes[order(up_genes$padj), ], 10)
   top_10_down <- head(down_genes[order(down_genes$padj), ], 10)
   
   # Print total DE and significant gene counts
   cat("Total DE genes (padj < 0.05):", sum(res$padj < padj_threshold, na.rm = TRUE), "\n")
   cat("Total significant genes (padj < 0.05 & |LFC| > 1):", nrow(sig_genes), "\n")
   cat("Total number of significantly upregulated genes:", nrow(up_genes), "\n")
   cat("Total number of significantly downregulated genes:", nrow(down_genes), "\n")
   # Print the results
    cat("\nTop 10 Significant Genes (Sorted by padj):\n")
    print(top_10_sig)
    
    cat("\nTop 10 Significantly Downregulated Genes:\n")
    ##top 10 Significantly Downregulated Genes:
      print(top_10_down)   
      
      # Check number of significant genes
       num_sig <- nrow(sig_genes)
      cat("Number of significant genes (padj <", padj_threshold, 
          "& |LFC| >", lfc_threshold, "):", num_sig, "\n")
      
      #####################################################################
      
      #Stroke vs Diabetes (TOP 10 GENES)
      res_sd <- results(dds, contrast = c("condition","stroke","DIABETES_TYPE1"))
      
      # Remove NA
      res_sd <- res_sd[!is.na(res_sd$padj), ]
      
      # Order by significance
      res_sd <- res_sd[order(res_sd$padj), ]
      
      # Top 10 genes
      top10_sd <- res_sd[1:10, ]
      
      # Print
      print(top10_sd)
      
      sig_genes <- res_sd[res_sd$padj < 0.05 & abs(res_sd$log2FoldChange) > 1, ]
      up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
      
      up_genes <- up_genes[order(up_genes$padj), ]
      
      top10_up <- up_genes[1:10, ]
      
      print(top10_up)
      down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]
      
      down_genes <- down_genes[order(down_genes$padj), ]
      
      top10_down <- down_genes[1:10, ]
      
      print(top10_down)
#####################################################################
      #Diabetes vs Control (TOP 10 GENES)
      res_dc <- results(dds, contrast = c("condition","DIABETES_TYPE1","Healthy_control"))
      
      # Remove NA
      res_dc <- res_dc[!is.na(res_dc$padj), ]
      
      # Order
      res_dc <- res_dc[order(res_dc$padj), ]
      
      # Top 10
      top10_dc <- res_dc[1:10, ]
      
      # Print
      print(top10_dc)
    #
      res_dc <- results(dds, contrast = c("condition","DIABETES_TYPE1","Healthy_control"))
      
      res_dc <- res_dc[!is.na(res_dc$padj), ]
      sig_genes_dc <- res_dc[res_dc$padj < 0.05 & abs(res_dc$log2FoldChange) > 1, ]
      ##Upregulated (Diabetes ↑)
      up_dc <- sig_genes_dc[sig_genes_dc$log2FoldChange > 0, ]
      
      up_dc <- up_dc[order(up_dc$padj), ]
      
      top10_up_dc <- up_dc[1:10, ]
      
      print(top10_up_dc)
      #Downregulated (Diabetes ↓)
      down_dc <- sig_genes_dc[sig_genes_dc$log2FoldChange < 0, ]
      
      down_dc <- down_dc[order(down_dc$padj), ]
      
      top10_down_dc <- down_dc[1:10, ]
      
      print(top10_down_dc)
      
      ####################3
      as.data.frame(top10_sd)
      as.data.frame(top10_dc)
      
      ## csv file createdd 
      write.csv(as.data.frame(top10_sd), "top10_stroke_vs_diabetes.csv")
      write.csv(as.data.frame(top10_dc), "top10_diabetes_vs_control.csv")

###################################################@#########################
##
     s
########################################################################################
      ## stastical testand anova  
      ## t-test 
      library(DESeq2)
      library(SummarizedExperiment)
      
      vsd_mat <- assay(vsd)
      vsd_mat <- assay(vsd)
      
      # Example gene
      gene <- "NHERF2"
      
      group <- colData(dds)$condition
      
      # Extract values
      stroke_vals <- vsd_mat[gene, group == "stroke"]
      control_vals <- vsd_mat[gene, group == "healthy.control"]
      
      # t-test
      t.test(stroke_vals, control_vals)    
      
    ####### ANOVA 3 conditions   
      vsd_mat <- assay(vsd)
      group <- colData(dds)$condition
      
      # Example gene
      gene <- "NHERF2"
      
      expr <- vsd_mat[gene, ]
      
      # Data frame
      df <- data.frame(
        expression = expr,
        condition = group
      )
      
      # ANOVA
      anova_res <- aov(expression ~ condition, data = df)
      summary(anova_res)
      
      TukeyHSD(anova_res)

      ####Loop for multiple genes (Top10)
      top_genes <- rownames(res_ordered)[1:10]
      
      for (gene in top_genes) {
        
        expr <- vsd_mat[gene, ]
        df <- data.frame(expression = expr, condition = group)
        
        res <- summary(aov(expression ~ condition, data = df))
        
        cat("\nGene:", gene, "\n")
        print(res)
      }
      
######## Boxplot + stats (best for report)
      
      
#################################################################

      dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
      res_lrt <- results(dds_lrt)
      
      res_lrt <- res_lrt[!is.na(res_lrt$padj), ]
      res_lrt <- res_lrt[order(res_lrt$padj), ]
      
      top_genes <- rownames(res_lrt)[1:50]
      
      
      sum(is.na(mat))
      sum(is.nan(mat))
      sum(is.infinite(mat))
      mat <- mat[complete.cases(mat), ]
      mat[is.na(mat)] <- 0
      mat[is.nan(mat)] <- 0
      mat[is.infinite(mat)] <- 0
      
      dds <- dds[rowSums(counts(dds)) > 10, ]
      summary(res_lrt)
######## LTR VALCONO PLOT 
      library(ggplot2)
      
      df <- as.data.frame(res_lrt)
      
      ggplot(df, aes(log2FoldChange, -log10(pvalue))) +
        geom_point(alpha=0.5) +
        ggtitle("LRT Volcano Plot")   
      
      
########################################################################################################
### MACHINE LEARNING MODEL 
  # crated meta data for python export   
      # VST matrix
      vsd_mat <- assay(vsd)
      
      # Save as CSV
      write.csv(vsd_mat, "vst_matrix.csv")
      
      # Metadata
      write.csv(as.data.frame(colData(dds)), "metadata.csv")
################################################################################### 
      # Expression matrix
      vsd_mat <- assay(vsd)
      
      # Transpose (ML ke liye samples rows me hone chahiye)
      ml_data <- t(vsd_mat)
      
      # Labels
      labels <- colData(dds)$condition
      
      df <- data.frame(ml_data)
      df$condition <- labels
      

################################################################################################
#### correction heatmap and top correction genes 
   
      library(ggplot2)
      
      ggplot(plot_df, aes(x=x, y=y)) +
        
        geom_point(color="black", size=2) +
        
        geom_smooth(method="lm", color="blue", fill="grey", alpha=0.3) +
        
        facet_wrap(~pair, scales="free", ncol=3) +
        
        theme_minimal(base_size = 14) +
        
        labs(
          title = "Correlation of Top 10 Genes",
          x = "Gene Expression",
          y = "Gene Expression"
        ) +
        
        theme(
          plot.title = element_text(hjust=0.5, face="bold")
        )
  
      pairs <- combn(colnames(df), 2)
      df <- as.data.frame(t(vsd_mat[common_genes, ]))
      colnames(df)
      pairs <- combn(common_genes, 2, simplify = FALSE)
      pairs <- pairs[1:6]
      
    ##################Top genes (Stroke vs Diabetes)
     
      res_sd <- results(dds, contrast = c("condition","stroke","DIABETES_TYPE1"))
      
      res_sd <- res_sd[!is.na(res_sd$padj), ]
      res_sd <- res_sd[order(res_sd$padj), ]
      
      top10_sd <- rownames(res_sd)[1:10]
      
      vsd_mat <- assay(vsd)
      
      # match genes
      common_genes <- intersect(top10_sd, rownames(vsd_mat))
      
      # use ONLY these genes
      df <- as.data.frame(t(vsd_mat[common_genes, ]))
      
      ##############Create gene pairs###
      pairs <- combn(common_genes, 2, simplify = FALSE)
      
      # take only first 6 (clean plot)
      pairs <- pairs[1:6]
     
       ##Build plotting data
      plot_list <- lapply(pairs, function(p){
        data.frame(
          x = df[[p[1]]],
          y = df[[p[2]]],
          pair = paste(p[1], "vs", p[2])
        )
      })
      
      plot_df <- do.call(rbind, plot_list)
      ## scatter plot top 10 gene 
      library(ggplot2)
      
      ggplot(plot_df, aes(x=x, y=y)) +
        
        geom_point(color="black", size=2) +
        
        geom_smooth(method="lm", color="blue", fill="grey", alpha=0.3) +
        
        facet_wrap(~pair, scales="free", ncol=3) +
        
        theme_minimal(base_size = 14) +
        
        labs(
          title = "Scatter Plot of Top 10 Genes (Stroke vs Diabetes)",
          x = "Gene Expression",
          y = "Gene Expression"
        ) +
        
        theme(
          plot.title = element_text(hjust=0.5, face="bold")
        )
 #############################################################################################################
      # Use only genes that exist in matrix
      vsd_mat <- assay(vsd)
      
      available_genes <- rownames(vsd_mat)
      
      # Take top 10 from available genes (safe)
      top10_safe <- head(available_genes, 10)
      
      print(top10_safe)
      df <- as.data.frame(t(vsd_mat[top10_safe, ]))
      pairs <- combn(colnames(df), 2, simplify = FALSE)
      pairs <- pairs[1:6]   # keep it small
      
      plot_list <- lapply(pairs, function(p){
        data.frame(
          x = df[[p[1]]],
          y = df[[p[2]]],
          pair = paste(p[1], "vs", p[2])
        )
      })
      
      plot_df <- do.call(rbind, plot_list)
      
      head(plot_df)   # ✔ should show data
          
      library(ggplot2)
      
      ggplot(plot_df, aes(x=x, y=y)) +
        geom_point(color="black", size=2) +
        geom_smooth(method="lm", color="blue", fill="grey", alpha=0.3) +
        facet_wrap(~pair, scales="free", ncol=3) +
        theme_minimal(base_size = 14) +
        labs(
          title = "Correlation of Top Genes",
          x = "Gene Expression",
          y = "Gene Expression"
        ) +
        theme(
          plot.title = element_text(hjust=0.5, face="bold")
        )
##(Stroke vs Diabetes)
      res_sd <- results(dds, contrast = c("condition","stroke","DIABETES_TYPE1"))
      
      res_sd <- res_sd[!is.na(res_sd$padj), ]
      res_sd <- res_sd[order(res_sd$padj), ]
      
      top10_sd <- rownames(res_sd)[1:10]
      vsd_mat <- assay(vsd)
      
      common_genes <- intersect(top10_sd, rownames(vsd_mat))
      
      df <- as.data.frame(t(vsd_mat[common_genes, ]))
      vsd_mat <- assay(vsd)
      
      common_genes <- intersect(top10_sd, rownames(vsd_mat))
      
      df <- as.data.frame(t(vsd_mat[common_genes, ]))
      
      common_genes <- intersect(top10_sd, colnames(df))
      
      # or safer:
      common_genes <- intersect(top10_sd, names(df))
      
##top coreection heatmap 
      
      
      install.packages("pheatmap")  # if not installed
      library(pheatmap)
      
      pheatmap(
        cor_mat,
        
        # Show values
        display_numbers = TRUE,
        
        # Color scale (blue → white → red)
        color = colorRampPalette(c("blue","white","red"))(100),
        
        # Clustering
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        
        # Labels
        main = "Correlation Heatmap (Top 10 Genes)",
        
        fontsize = 10,
        border_color = NA
      )
##########################################################################################3
    
      common_genes <- intersect(top10_genes, rownames(vsd_mat))
      res_sd <- results(dds, contrast = c("condition","stroke","DIABETES_TYPE1"))
      
      res_sd <- res_sd[!is.na(res_sd$padj), ]
      res_sd <- res_sd[order(res_sd$padj), ]
      
      top10_sd <- rownames(res_sd)[1:10]
      
      print(top10_sd)  
      common_genes <- intersect(top10_sd, rownames(vsd_mat))
      
      print(common_genes)
      mat <- vsd_mat[common_genes, ]
      df <- as.data.frame(t(mat))
      
      cor_mat <- cor(df)
      
      pheatmap(cor_mat, display_numbers=TRUE)
      
  pheatmap(
  cor_mat,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue","white","red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation Heatmap (Stroke vs Diabetes - Top 10 Genes)",
  fontsize = 10,
  fontsize_number = 8,
  border_color = NA,
  angle_col = 45
)
#######################################################################################
  library(pheatmap)
      
      pheatmap(
        cor_mat,
        display_numbers = TRUE,
        color = colorRampPalette(c("blue","white","red"))(100),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        main = "Correlation Heatmap (Stroke vs Diabetes - Top 10 Genes)",
        fontsize = 10
      )
      
###################################################################################################
      install.packages("ggpubr")   # only first time
      library(ggpubr) 
      
      ggplot(plot_df, aes(x=x, y=y)) +
        
        geom_point(color="black", size=2) +
        
        geom_smooth(method="lm", color="blue", fill="grey", alpha=0.3) +
        
        stat_cor(
          method = "pearson",
          label.x.npc = "left",
          label.y.npc = "top",
          size = 4
        ) +
        
        facet_wrap(~pair, scales="free", ncol=3) +
        
        theme_minimal(base_size = 14) +
        
        labs(
          title = "Correlation of Top 10 Genes (Stroke vs Diabetes)",
          x = "Gene Expression",
          y = "Gene Expression"
        ) +
        
        theme(
          plot.title = element_text(hjust=0.5, face="bold")
        )
      
      ggplot(plot_df, aes(x=x, y=y)) +
        geom_point(color="black", size=2) +
        geom_smooth(method="lm", color="blue", fill="grey", alpha=0.3) +
        stat_cor(method="pearson") +
        facet_wrap(~pair)
      
      install.packages("ggpubr", dependencies = TRUE)
      library(ggpubr)
      