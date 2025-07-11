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

# 加载目录结构创建函数
source("code/create_output_structure.R")

# 获取数据集路径信息
dataset_paths <- get_dataset_paths()
dataset_name <- dataset_paths$dataset_name

# 创建输出目录结构
output_paths <- create_output_structure()

# ====================== 检查输入文件 ======================
# 检查生成的数据文件是否存在
input_file <- dataset_paths$generated_data_file
if (!file.exists(input_file)) {
    cat("❌ 生成的数据文件不存在:", input_file, "\n")
    cat("请先运行 scDDPM.py 生成数据\n")
    stop("文件不存在")
}

# ====================== Load Expression Data ======================
# Load generated gene expression matrix (rows: genes, columns: cells)
counts <- read.csv(input_file, row.names = 1)
cat("✅ Raw matrix dimensions: ", dim(counts)[1], "genes ×", dim(counts)[2], "cells\n")

# 检查数据格式：原始数据是细胞为行，基因为列，需要转置为基因为行，细胞为列
cat("📊 原始数据维度:", dim(counts)[1], "行 ×", dim(counts)[2], "列\n")

# 检查是否有label列需要移除
if ("label" %in% colnames(counts)) {
    cat("📋 移除label列\n")
    labels <- counts$label
    counts <- counts[, !colnames(counts) %in% "label"]
} else {
    labels <- rep("unknown", nrow(counts))
}

# Select first 100 cells for visualization/clustering
if (nrow(counts) > 100) {
    counts <- counts[1:100, ]
    labels <- labels[1:100]
    cat("✅ 选择前100个细胞进行分析\n")
}

# 转置矩阵：SC3需要基因为行，细胞为列
counts <- t(counts)
cat("✅ 转置后矩阵维度:", dim(counts)[1], "基因 ×", dim(counts)[2], "细胞\n")

# 创建细胞注释（现在细胞是列）
annotation <- data.frame(
    cell_type = labels,
    row.names = colnames(counts)
)

# 保存注释文件
annotation_file <- get_output_path("clustering", "cell_annotation.csv")
write.csv(annotation, annotation_file)
cat("✅ 细胞注释已保存到:", annotation_file, "\n")

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
output_file <- get_output_path("clustering", "sc3_clustering_labels.csv")
write.csv(colData(sce), output_file)
cat("✅ SC3 clustering labels saved to '", output_file, "'\n")
