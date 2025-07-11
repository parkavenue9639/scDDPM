# ================== Environment Setup ==================
rm(list = ls())

# Required packages (uncomment to install if needed)
# install.packages("Seurat")
# devtools::install_version("Seurat", version = "4.1.4")
# BiocManager::install("Seurat")
# install.packages("Matrix")

# ================== Configuration ==================
# 设置要处理的数据集名称
dataset_name <- Sys.getenv("DATASET_NAME", "AD01103")  # 从环境变量读取数据集名称，默认AD01103

# 设置输入和输出路径
input_dir <- paste0("data/", dataset_name, "/")
output_dir <- "FD1000/"

# 创建输出目录
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ================== Load Data ==================
# Load gene expression matrix (rows: genes, columns: cells)
expr_file <- paste0(input_dir, dataset_name, "_expr.txt")
if (!file.exists(expr_file)) {
    stop(paste("❌ 文件不存在:", expr_file, "\n请先运行解压脚本或检查文件路径"))
}

data <- read.delim(expr_file, header = TRUE, row.names = 1)

# Calculate non-zero rate
completed_nz <- sum(data > 0) / (nrow(data) * ncol(data))
cat("✅ Non-zero rate:", round(completed_nz * 100, 2), "%\n")

library(Seurat)
library(SeuratObject)
library(magrittr)
library(dplyr)

# Basic statistics
cat("🔢 Dimensions: ", dim(data)[1], "genes ×", dim(data)[2], "cells\n")
cat("🧪 Total zero entries:", sum(data == 0), "\n")

# ================== Create Seurat Object ==================
# Rows: genes, Columns: cells. Apply minimum cell filter (5%)
scell <- CreateSeuratObject(counts = data, project = "Alzheimer0205",
                            min.cells = ceiling(0.05 * ncol(data)))
data <- 0  # Free memory

# Save gene and cell names
gene_names <- rownames(scell)
cell_names <- colnames(scell)
cat("🧬 Filtered dimensions:", dim(scell@assays$RNA@layers$counts), "\n")

# ================== Normalize Data ==================
# Log-normalization (CPM → log1p)
data <- NormalizeData(scell, normalization.method = "LogNormalize", scale.factor = 10000)
norm <- data@assays$RNA@layers$data
normalized_df <- as.data.frame(norm)
rownames(normalized_df) <- gene_names
colnames(normalized_df) <- cell_names

# ================== Select Top 1000 Highly Variable Genes ==================
# Select top 1000 variable genes (vst method)
data1000 <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 1000)
top1000 <- head(VariableFeatures(data1000), 1000)
pb1000 <- normalized_df[top1000, ]
pb1000 <- t(pb1000)
cat("✅ Selected top 1000 variable genes:", dim(pb1000), "\n")

# Save processed matrix (expression only)
output_file_1000 <- paste0(output_dir, dataset_name, "PrePro1000.csv")
write.csv(pb1000, output_file_1000)
cat("✅ 预处理数据已保存到:", output_file_1000, "\n")

# ================== Merge with Cell Labels ==================
label_file <- paste0(input_dir, dataset_name, "_cell_label.txt")
if (!file.exists(label_file)) {
    stop(paste("❌ 标签文件不存在:", label_file))
}

label <- read.delim(label_file, header = TRUE, row.names = 1)
merged_data_1000 <- merge(pb1000, label, by = "row.names", all.x = TRUE)

cat("🔖 Label distribution:\n")
print(table(label))

# Save final labeled matrix
output_file_labeled_1000 <- paste0(output_dir, dataset_name, "PreProLabel1000.csv")
write.csv(merged_data_1000, output_file_labeled_1000, row.names = FALSE)
cat("✅ 标记数据已保存到:", output_file_labeled_1000, "\n")

cat("✅ 处理完成。数据集:", dataset_name, "\n")
