
count_dir <- "/home/vajid/InternLab/ncbi_dataset/data/GCF_000001405.40/counts/"

count_files <- list.files(count_dir, pattern="*_Counts.txt", full.names=TRUE)

count_list <- lapply(count_files, function(file) {
  sample <- gsub("_Counts.txt","",basename(file))
  fread(file, col.names=c("Gene", sample))
})

countData <- Reduce(function(x,y) merge(x,y,by="Gene",all=TRUE), count_list)
countData[is.na(countData)] <- 0
countData <- as.data.frame(countData)
rownames(countData) <- countData$Gene
countData$Gene <- NULL
countData <- as.matrix(countData)

# Metadata preparation (2 samples per condition)
colData <- data.frame(
  condition = factor(c(
    "stroke", "stroke",
    "DIABETES_TYPE1", "DIABETES_TYPE1",
    "Healthy_control", "Healthy_control"
  ))
)

rownames(colData) <- c(
  "SRR16638184",
  "SRR16638186",
  "SRR11806265",
  "SRR11806266",
  "SRR16638202",
  "SRR16638201"
)


colData <- colData[colnames(countData), , drop = FALSE]
all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ condition
)



