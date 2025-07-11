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

# 加载目录结构创建函数
source("code/create_output_structure.R")

# 获取数据集路径信息
dataset_paths <- get_dataset_paths()
dataset_name <- dataset_paths$dataset_name

# 创建输出目录结构
output_paths <- create_output_structure()

# ====================== 检查文件是否存在 ======================
# 检查真实数据文件
real_data_file <- dataset_paths$real_data_file
if (!file.exists(real_data_file)) {
    cat("❌ 真实数据文件不存在:", real_data_file, "\n")
    cat("请先运行 Preprocess.R 生成预处理数据\n")
    stop("文件不存在")
}

# 检查生成数据文件
gen_data_file <- dataset_paths$generated_data_file
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

# 检查原始标签
cat("🔍 原始标签检查:\n")
cat("  真实数据标签类型:", class(org_label$cell_type), "\n")
cat("  真实数据唯一标签:", unique(org_label$cell_type), "\n")
cat("  生成数据标签类型:", class(gen_label$cell_type), "\n")
cat("  生成数据唯一标签:", unique(gen_label$cell_type), "\n")

# 正确处理标签 - 转换为因子然后转换为数字
org_label$cell_type_factor <- as.factor(org_label$cell_type)
org_label$cell_type_numeric <- as.numeric(org_label$cell_type_factor)

gen_label$cell_type_factor <- as.factor(gen_label$cell_type)
gen_label$cell_type_numeric <- as.numeric(gen_label$cell_type_factor)

# 为了区分真实数据和生成数据，给生成数据标签加偏移
gen_label$cell_type_numeric <- gen_label$cell_type_numeric + max(org_label$cell_type_numeric)

cat("  处理后标签范围 - 真实数据:", range(org_label$cell_type_numeric), "\n")
cat("  处理后标签范围 - 生成数据:", range(gen_label$cell_type_numeric), "\n")

# Standardize column names and merge data
# 移除标签列，只保留表达数据
org_expr <- org_counts[, !colnames(org_counts) %in% c("label", "Cell")]
gen_expr <- gen_counts[, !colnames(gen_counts) %in% c("label", "Cell")]

# 确保列名一致
common_genes <- intersect(colnames(org_expr), colnames(gen_expr))
org_expr <- org_expr[, common_genes]
gen_expr <- gen_expr[, common_genes]

all_counts <- rbind(org_expr, gen_expr)
# 创建合并标签数据框
annotation <- data.frame(
    cell_type = c(org_label$cell_type, gen_label$cell_type),
    cell_type_numeric = c(org_label$cell_type_numeric, gen_label$cell_type_numeric),
    row.names = c(rownames(org_label), rownames(gen_label))
)

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

# 检查PCA结果
cat("  - PCA计算完成，细胞数:", ncol(sce), "基因数:", nrow(sce), "\n")
pca_data <- reducedDim(PCAsce, "PCA")
cat("  - PCA数据维度:", dim(pca_data), "\n")

# 详细诊断PCA数据
cat("  🔍 PCA数据诊断:\n")
cat("    - PC1 范围:", range(pca_data[, 1], na.rm = TRUE), "\n")
cat("    - PC2 范围:", range(pca_data[, 2], na.rm = TRUE), "\n")
cat("    - PC1 是否有NA:", any(is.na(pca_data[, 1])), "\n")
cat("    - PC2 是否有NA:", any(is.na(pca_data[, 2])), "\n")
cat("    - PC1 前5个值:", head(pca_data[, 1], 5), "\n")
cat("    - PC2 前5个值:", head(pca_data[, 2], 5), "\n")

# 准备颜色标签
cell_colors <- org_label$cell_type_numeric[1:ncol(counts)]
cat("  - 颜色标签长度:", length(cell_colors), "\n")
cat("  - 颜色标签范围:", range(cell_colors, na.rm = TRUE), "\n")
cat("  - 颜色标签前10个:", head(cell_colors, 10), "\n")

# 保存PCA图 - 方法1: 基础绘图
pdf(get_output_path("visualization", "real_data_pca.pdf"))
plot(pca_data[, 1], pca_data[, 2], 
     col = cell_colors, pch = 19, cex = 1.8, 
     main = "PCA - Real Data", xlab = "PC1", ylab = "PC2",
     xlim = range(pca_data[, 1], na.rm = TRUE),
     ylim = range(pca_data[, 2], na.rm = TRUE))
legend("topright", legend = unique(org_label$cell_type), 
       col = unique(cell_colors), pch = 19, title = "Cell Type")
dev.off()

# 保存PCA图 - 方法2: 使用ggplot2（更可靠）
pca_df <- data.frame(
  PC1 = pca_data[, 1],
  PC2 = pca_data[, 2],
  cell_type = as.factor(org_label$cell_type[1:ncol(counts)])
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA - Real Data", 
       x = "PC1", y = "PC2", 
       color = "Cell Type") +
  theme(legend.position = "right")

ggsave(get_output_path("visualization", "real_data_pca_ggplot.pdf"), plot = p, width = 10, height = 8)
cat("  ✅ ggplot2版本已保存: real_data_pca_ggplot.pdf\n")

# -------- PCA visualization on generated data --------
cat("📊 生成合成数据PCA图...\n")
counts <- t(gen_expr[, 1:min(1000, ncol(gen_expr))])
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts = as.matrix(counts)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
PCAsce <- runPCA(sce)

# 检查PCA结果
cat("  - PCA计算完成，细胞数:", ncol(sce), "基因数:", nrow(sce), "\n")
pca_data <- reducedDim(PCAsce, "PCA")
cat("  - PCA数据维度:", dim(pca_data), "\n")

# 详细诊断PCA数据
cat("  🔍 PCA数据诊断:\n")
cat("    - PC1 范围:", range(pca_data[, 1], na.rm = TRUE), "\n")
cat("    - PC2 范围:", range(pca_data[, 2], na.rm = TRUE), "\n")

# 准备颜色标签
cell_colors <- gen_label$cell_type_numeric[1:ncol(counts)]
cat("  - 颜色标签长度:", length(cell_colors), "\n")

pdf(get_output_path("visualization", "generated_data_pca.pdf"))
plot(pca_data[, 1], pca_data[, 2],
     col = cell_colors, pch = 19, cex = 1.0,
     main = "PCA - Generated Data", xlab = "PC1", ylab = "PC2",
     xlim = range(pca_data[, 1], na.rm = TRUE),
     ylim = range(pca_data[, 2], na.rm = TRUE))
legend("topright", legend = unique(gen_label$cell_type), 
       col = unique(cell_colors), pch = 19, title = "Cell Type")
dev.off()

# ggplot2版本
pca_df <- data.frame(
  PC1 = pca_data[, 1],
  PC2 = pca_data[, 2],
  cell_type = as.factor(gen_label$cell_type[1:ncol(counts)])
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA - Generated Data", 
       x = "PC1", y = "PC2", 
       color = "Cell Type") +
  theme(legend.position = "right")

ggsave(get_output_path("visualization", "generated_data_pca_ggplot.pdf"), plot = p, width = 10, height = 8)
cat("  ✅ ggplot2版本已保存: generated_data_pca_ggplot.pdf\n")

# -------- PCA on combined real and generated data --------
cat("📊 生成合并数据PCA图...\n")
counts <- t(all_counts[, 1:min(1000, ncol(all_counts))])
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts), logcounts = as.matrix(counts)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
PCAsce <- runPCA(sce)

# 检查PCA结果
cat("  - PCA计算完成，细胞数:", ncol(sce), "基因数:", nrow(sce), "\n")
pca_data <- reducedDim(PCAsce, "PCA")
cat("  - PCA数据维度:", dim(pca_data), "\n")

# 详细诊断PCA数据
cat("  🔍 PCA数据诊断:\n")
cat("    - PC1 范围:", range(pca_data[, 1], na.rm = TRUE), "\n")
cat("    - PC2 范围:", range(pca_data[, 2], na.rm = TRUE), "\n")

# 准备颜色标签
cell_colors <- annotation$cell_type_numeric[1:ncol(counts)]
cat("  - 颜色标签长度:", length(cell_colors), "\n")

pdf(get_output_path("visualization", "combined_data_pca.pdf"))
plot(pca_data[, 1], pca_data[, 2],
     col = cell_colors, pch = 19, cex = 1.8,
     main = "PCA - Combined Data", xlab = "PC1", ylab = "PC2",
     xlim = range(pca_data[, 1], na.rm = TRUE),
     ylim = range(pca_data[, 2], na.rm = TRUE))
legend("topright", legend = unique(annotation$cell_type), 
       col = unique(cell_colors), pch = 19, title = "Cell Type")
dev.off()

# ggplot2版本 - 带数据类型标识
pca_df <- data.frame(
  PC1 = pca_data[, 1],
  PC2 = pca_data[, 2],
  cell_type = as.factor(annotation$cell_type[1:ncol(counts)]),
  data_type = c(rep("Real", nrow(org_expr)), rep("Generated", nrow(gen_expr)))[1:ncol(counts)]
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type, shape = data_type)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA - Combined Data (Real vs Generated)", 
       x = "PC1", y = "PC2", 
       color = "Cell Type", shape = "Data Type") +
  theme(legend.position = "right") +
  scale_shape_manual(values = c("Real" = 16, "Generated" = 17))

ggsave(get_output_path("visualization", "combined_data_pca_ggplot.pdf"), plot = p, width = 12, height = 8)
cat("  ✅ ggplot2版本已保存: combined_data_pca_ggplot.pdf\n")

cat("✅ 可视化完成！结果保存在", output_paths$visualization, "目录\n")
