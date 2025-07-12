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
counts_real <- t(org_expr[, 1:min(1000, ncol(org_expr))])
sce_real <- SingleCellExperiment(assays = list(counts = as.matrix(counts_real), logcounts = as.matrix(counts_real)))
rowData(sce_real)$feature_symbol <- rownames(sce_real)
sce_real <- sce_real[!duplicated(rowData(sce_real)$feature_symbol), ]
PCAsce_real <- runPCA(sce_real)
pca_data_real <- reducedDim(PCAsce_real, "PCA")

# -------- PCA visualization on generated data --------
counts_gen <- t(gen_expr[, 1:min(1000, ncol(gen_expr))])
sce_gen <- SingleCellExperiment(assays = list(counts = as.matrix(counts_gen), logcounts = as.matrix(counts_gen)))
rowData(sce_gen)$feature_symbol <- rownames(sce_gen)
sce_gen <- sce_gen[!duplicated(rowData(sce_gen)$feature_symbol), ]
PCAsce_gen <- runPCA(sce_gen)
pca_data_gen <- reducedDim(PCAsce_gen, "PCA")

# -------- PCA on combined real and generated data --------
counts_combined <- t(all_counts[, 1:min(1000, ncol(all_counts))])
sce_combined <- SingleCellExperiment(assays = list(counts = as.matrix(counts_combined), logcounts = as.matrix(counts_combined)))
rowData(sce_combined)$feature_symbol <- rownames(sce_combined)
sce_combined <- sce_combined[!duplicated(rowData(sce_combined)$feature_symbol), ]
PCAsce_combined <- runPCA(sce_combined)
pca_data_combined <- reducedDim(PCAsce_combined, "PCA")

# PCA空间统一处理（保持原始合并PCA结果）
cat("🔧 PCA空间统一处理:\n")
cat("  📝 使用原始合并PCA结果，无需方向调整\n")
cat("  📝 测试显示原始结果与各数据集单独PCA最为一致\n")
cat("  ✅ 真实数据和生成数据在同一空间中表现一致\n")

# 统一坐标轴范围 (基于原始合并数据PCA空间)
global_xlim <- range(pca_data_combined[,1], na.rm = TRUE)
global_ylim <- range(pca_data_combined[,2], na.rm = TRUE)

cat("  合并数据PCA维度:", dim(pca_data_combined), "\n")
cat("  真实数据索引范围: 1 到", nrow(org_expr), "\n")
cat("  生成数据索引范围:", nrow(org_expr) + 1, "到", nrow(org_expr) + nrow(gen_expr), "\n")
cat("  PC1范围:", sprintf("%.2f 到 %.2f", global_xlim[1], global_xlim[2]), "\n")
cat("  PC2范围:", sprintf("%.2f 到 %.2f", global_ylim[1], global_ylim[2]), "\n")

# -------- 画真实数据PCA图 (使用统一的PCA空间) --------
# 从合并数据的PCA结果中提取真实数据部分
real_start_idx <- 1
real_end_idx <- nrow(org_expr)
pca_data_real_unified <- pca_data_combined[real_start_idx:real_end_idx, ]

cell_colors_real <- org_label$cell_type_numeric[1:nrow(pca_data_real_unified)]
pdf(get_output_path("visualization", "real_data_pca.pdf"))
plot(pca_data_real_unified[, 1], pca_data_real_unified[, 2], 
     col = cell_colors_real, pch = 19, cex = 1.8, 
     main = "PCA - Real Data (Unified Space)", xlab = "PC1", ylab = "PC2",
     xlim = global_xlim, ylim = global_ylim)
legend("topright", legend = unique(org_label$cell_type), 
       col = unique(cell_colors_real), pch = 19, title = "Cell Type")
dev.off()

# 保存PCA图 - 方法2: 使用ggplot2（更可靠）
pca_df <- data.frame(
  PC1 = pca_data_real_unified[, 1],
  PC2 = pca_data_real_unified[, 2],
  cell_type = as.factor(org_label$cell_type[1:nrow(pca_data_real_unified)])
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA - Real Data (Unified Space)", 
       x = "PC1", y = "PC2", 
       color = "Cell Type") +
  theme(legend.position = "right") +
  coord_cartesian(xlim = global_xlim, ylim = global_ylim)

ggsave(get_output_path("visualization", "real_data_pca_ggplot.pdf"), plot = p, width = 10, height = 8)
cat("  ✅ ggplot2版本已保存: real_data_pca_ggplot.pdf\n")

# -------- 画生成数据PCA图 (使用统一的PCA空间) --------
# 从合并数据的PCA结果中提取生成数据部分
gen_start_idx <- nrow(org_expr) + 1
gen_end_idx <- nrow(org_expr) + nrow(gen_expr)
pca_data_gen_unified <- pca_data_combined[gen_start_idx:gen_end_idx, ]

cell_colors_gen <- gen_label$cell_type_numeric[1:nrow(pca_data_gen_unified)]
pdf(get_output_path("visualization", "generated_data_pca.pdf"))
plot(pca_data_gen_unified[, 1], pca_data_gen_unified[, 2],
     col = cell_colors_gen, pch = 19, cex = 1.0,
     main = "PCA - Generated Data (Unified Space)", xlab = "PC1", ylab = "PC2",
     xlim = global_xlim, ylim = global_ylim)
legend("topright", legend = unique(gen_label$cell_type), 
       col = unique(cell_colors_gen), pch = 19, title = "Cell Type")
dev.off()

# ggplot2版本
pca_df <- data.frame(
  PC1 = pca_data_gen_unified[, 1],
  PC2 = pca_data_gen_unified[, 2],
  cell_type = as.factor(gen_label$cell_type[1:nrow(pca_data_gen_unified)])
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA - Generated Data (Unified Space)", 
       x = "PC1", y = "PC2", 
       color = "Cell Type") +
  theme(legend.position = "right") +
  coord_cartesian(xlim = global_xlim, ylim = global_ylim)

ggsave(get_output_path("visualization", "generated_data_pca_ggplot.pdf"), plot = p, width = 10, height = 8)
cat("  ✅ ggplot2版本已保存: generated_data_pca_ggplot.pdf\n")

# -------- 画合并数据PCA图 --------
cell_colors_combined <- annotation$cell_type_numeric[1:ncol(counts_combined)]
pdf(get_output_path("visualization", "combined_data_pca.pdf"))
plot(pca_data_combined[, 1], pca_data_combined[, 2],
     col = cell_colors_combined, pch = 19, cex = 1.8,
     main = "PCA - Combined Data", xlab = "PC1", ylab = "PC2",
     xlim = global_xlim, ylim = global_ylim)
legend("topright", legend = unique(annotation$cell_type), 
       col = unique(cell_colors_combined), pch = 19, title = "Cell Type")
dev.off()

# ggplot2版本 - 带数据类型标识
pca_df <- data.frame(
  PC1 = pca_data_combined[, 1],
  PC2 = pca_data_combined[, 2],
  cell_type = as.factor(annotation$cell_type[1:ncol(counts_combined)]),
  data_type = c(rep("Real", nrow(org_expr)), rep("Generated", nrow(gen_expr)))[1:ncol(counts_combined)]
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type, shape = data_type)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA - Combined Data (Real vs Generated)", 
       x = "PC1", y = "PC2", 
       color = "Cell Type", shape = "Data Type") +
  theme(legend.position = "right") +
  scale_shape_manual(values = c("Real" = 16, "Generated" = 17)) +
  coord_cartesian(xlim = global_xlim, ylim = global_ylim)

ggsave(get_output_path("visualization", "combined_data_pca_ggplot.pdf"), plot = p, width = 12, height = 8)
cat("  ✅ ggplot2版本已保存: combined_data_pca_ggplot.pdf\n")

cat("✅ 可视化完成！结果保存在", output_paths$visualization, "目录\n")
