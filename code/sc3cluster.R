# ====================== Package Setup ======================
# Uncomment the lines below if you haven't installed the required packages
# BiocManager::install("scater")
# BiocManager::install("SC3")

library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library(scater)
library(Rtsne)

# ====================== File Path Setup ======================
# 使用当前工作目录，不再硬编码路径
cat("当前工作目录:", getwd(), "\n")
rm(list = ls())

# ====================== 检查输入文件 ======================
# 检查生成的数据文件是否存在
input_file <- "output/AD01103_generated.csv"
if (!file.exists(input_file)) {
    cat("❌ 生成的数据文件不存在:", input_file, "\n")
    cat("请先运行 scDDPM.py 生成数据\n")
    stop("文件不存在")
}

# ====================== Load Expression Data ======================
# Load generated gene expression matrix (rows: genes, columns: cells)
counts <- read.csv(input_file, row.names = 1)
cat("✅ Raw matrix dimensions: ", dim(counts)[1], "genes ×", dim(counts)[2], "cells\n")

# Select first 100 cells for visualization/clustering
if (ncol(counts) > 100) {
    counts <- counts[, 1:100]
    cat("✅ 选择前100个细胞进行分析\n")
}

# 如果没有标签文件，创建一个简单的标签
if (!file.exists("gen_label2.csv")) {
    cat("⚠️  标签文件不存在，创建默认标签\n")
    annotation <- data.frame(
        cell_type = rep("unknown", ncol(counts)),
        row.names = colnames(counts)
    )
    write.csv(annotation, "gen_label2.csv")
} else {
    annotation <- read.csv("gen_label2.csv", row.names = 1)
}

# Transpose expression matrix (SC3 requires genes as rows, cells as columns)
#counts <- t(counts)
cat("✅ Transposed matrix dimensions:", dim(counts)[1], "cells ×", dim(counts)[2], "genes\n")

# ====================== Create SingleCellExperiment Object ======================
# Input matrix must be normalized and log-transformed
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(counts),
    logcounts = as.matrix(counts)
  ),
  colData = annotation
)

# Remove duplicated gene symbols
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
cat("✅ SC3 input matrix after removing duplicates:", dim(sce), "\n")

# ====================== Run SC3 Clustering ======================
# Perform SC3 clustering with k = 6 clusters
sce <- sc3(sce, ks = 6, biology = TRUE, gene_filter = FALSE)
cat("✅ SC3 clustering complete.\n")

# ====================== Save Results ======================
output_file <- "output/sc3_clustering_labels.csv"
write.csv(colData(sce), output_file)
cat("✅ SC3 clustering labels saved to '", output_file, "'\n")
