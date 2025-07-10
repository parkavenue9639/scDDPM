# This script performs visualization to assess the overlap between real and generated scRNA-seq data.
# It includes PCA, t-SNE, and UMAP dimensionality reduction on real, generated, and combined datasets.

# Load required packages
# install.packages("BiocManager")
# BiocManager::install("scater")
# BiocManager::install("SC3")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("umap")

library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library(scater)
library(Rtsne)
library(umap)

# Set working directory and clean environment
rm(list = ls())
cat("当前工作目录:", getwd(), "\n")

# ====================== 检查文件是否存在 ======================
# 检查真实数据文件
real_data_file <- "FD1000/AD01103PreProLabel1000.csv"
if (!file.exists(real_data_file)) {
    cat("❌ 真实数据文件不存在:", real_data_file, "\n")
    cat("请先运行 Preprocess.R 生成预处理数据\n")
    stop("文件不存在")
}

# 检查生成数据文件
gen_data_file <- "output/AD01103_generated.csv"
if (!file.exists(gen_data_file)) {
    cat("❌ 生成数据文件不存在:", gen_data_file, "\n")
    cat("请先运行 scDDPM.py 生成数据\n")
    stop("文件不存在")
}

# Load real and generated data
org_counts <- read.csv(real_data_file, row.names = 1)
gen_counts <- read.csv(gen_data_file, row.names = 1, header = TRUE)

# 创建标签数据
org_label <- data.frame(
    cell_type = org_counts$label,
    row.names = rownames(org_counts)
)

gen_label <- data.frame(
    cell_type = gen_counts$label,
    row.names = rownames(gen_counts)
)

# Convert labels to integers
org_label$cell_type <- as.integer(org_label$cell_type) - 1
gen_label$cell_type <- as.integer(gen_label$cell_type) + 5

# Standardize column names and merge data
# 移除标签列，只保留表达数据
org_expr <- org_counts[, !colnames(org_counts) %in% c("label", "Cell")]
gen_expr <- gen_counts[, !colnames(gen_counts) %in% c("label", "Cell")]

# 确保列名一致
common_genes <- intersect(colnames(org_expr), colnames(gen_expr))
org_expr <- org_expr[, common_genes]
gen_expr <- gen_expr[, common_genes]

all_counts <- rbind(org_expr, gen_expr)
annotation <- rbind(org_label, gen_label)

cat("✅ 数据加载完成\n")
cat("  真实数据维度:", dim(org_expr), "\n")
cat("  生成数据维度:", dim(gen_expr), "\n")
cat("  合并数据维度:", dim(all_counts), "\n")

# -------- PCA visualization on real data --------
cat("📊 生成真实数据PCA图...\n")
counts <- t(org_expr[, 1:min(1000, ncol(org_expr))])
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts = as.matrix(counts)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
PCAsce <- runPCA(sce)

# 保存PCA图
pdf("output/real_data_pca.pdf")
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[, 1:2],
     col = org_label$cell_type[1:ncol(counts)] + 2, asp = 0.4, pch = 19, cex = 1.8, 
     main = "PCA - Real Data", xlab = "PC1", ylab = "PC2")
dev.off()

# -------- PCA visualization on generated data --------
cat("📊 生成合成数据PCA图...\n")
counts <- t(gen_expr[, 1:min(1000, ncol(gen_expr))])
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts = as.matrix(counts)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
PCAsce <- runPCA(sce)

pdf("output/generated_data_pca.pdf")
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[, 1:2],
     col = gen_label$cell_type[1:ncol(counts)] + 1, asp = 0.4, pch = 19, cex = 1.0,
     main = "PCA - Generated Data", xlab = "PC1", ylab = "PC2")
dev.off()

# -------- PCA on combined real and generated data --------
cat("📊 生成合并数据PCA图...\n")
counts <- t(all_counts[, 1:min(1000, ncol(all_counts))])
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts = as.matrix(counts)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
PCAsce <- runPCA(sce)

pdf("output/combined_data_pca.pdf")
plot(PCAsce@int_colData@listData$reducedDims@listData$PCA[, 1:2],
     col = annotation$cell_type[1:ncol(counts)], asp = 0.4, pch = 19, cex = 1.8,
     main = "PCA - Combined Data", xlab = "PC1", ylab = "PC2")
dev.off()

cat("✅ 可视化完成！结果保存在 output/ 目录\n")
