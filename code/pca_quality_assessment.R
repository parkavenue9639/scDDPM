# PCA质量评估脚本 - 提供客观的数字化评估指标

library(SingleCellExperiment)
library(scater)

cat("🔍 PCA质量评估分析\n")
cat("==========================================\n")

# 加载目录结构创建函数
source("code/create_output_structure.R")

# 获取数据集路径信息
dataset_paths <- get_dataset_paths()
dataset_name <- dataset_paths$dataset_name

# 创建输出目录结构
output_paths <- create_output_structure()

# ========== 检查文件是否存在 ==========
real_data_file <- dataset_paths$real_data_file
gen_data_file <- dataset_paths$generated_data_file

if (!file.exists(real_data_file) || !file.exists(gen_data_file)) {
    cat("❌ 数据文件不存在，请先运行完整分析流程\n")
    stop("文件不存在")
}

# ========== 加载数据 ==========
org_counts <- read.csv(real_data_file, row.names = 1)
gen_counts <- read.csv(gen_data_file, row.names = 1, header = TRUE)

# 处理标签
org_label <- as.factor(org_counts$label)
gen_label <- as.factor(gen_counts$label)

# 移除标签列
org_expr <- org_counts[, !colnames(org_counts) %in% c("label", "Cell")]
gen_expr <- gen_counts[, !colnames(gen_counts) %in% c("label", "Cell")]

# 确保基因一致
common_genes <- intersect(colnames(org_expr), colnames(gen_expr))
org_expr <- org_expr[, common_genes]
gen_expr <- gen_expr[, common_genes]

# ========== 计算PCA ==========
cat("📊 计算PCA结果...\n")

# 真实数据PCA
counts_real <- t(org_expr[, 1:min(1000, ncol(org_expr))])
sce_real <- SingleCellExperiment(assays = list(counts = as.matrix(counts_real), 
                                              logcounts = as.matrix(counts_real)))
rowData(sce_real)$feature_symbol <- rownames(sce_real)
sce_real <- sce_real[!duplicated(rowData(sce_real)$feature_symbol), ]
pca_real <- runPCA(sce_real)
pca_data_real <- reducedDim(pca_real, "PCA")

# 生成数据PCA
counts_gen <- t(gen_expr[, 1:min(1000, ncol(gen_expr))])
sce_gen <- SingleCellExperiment(assays = list(counts = as.matrix(counts_gen),
                                             logcounts = as.matrix(counts_gen)))
rowData(sce_gen)$feature_symbol <- rownames(sce_gen)
sce_gen <- sce_gen[!duplicated(rowData(sce_gen)$feature_symbol), ]
pca_gen <- runPCA(sce_gen)
pca_data_gen <- reducedDim(pca_gen, "PCA")

# 合并数据PCA
all_expr <- rbind(org_expr, gen_expr)
counts_all <- t(all_expr[, 1:min(1000, ncol(all_expr))])
sce_all <- SingleCellExperiment(assays = list(counts = as.matrix(counts_all),
                                             logcounts = as.matrix(counts_all)))
rowData(sce_all)$feature_symbol <- rownames(sce_all)
sce_all <- sce_all[!duplicated(rowData(sce_all)$feature_symbol), ]
pca_all <- runPCA(sce_all)
pca_data_all <- reducedDim(pca_all, "PCA")

# ========== 评估指标计算 ==========

# 1. PCA解释方差比例
cat("\n📈 PCA解释方差分析:\n")
variance_real <- attr(pca_data_real, "percentVar")[1:2]
variance_gen <- attr(pca_data_gen, "percentVar")[1:2]
variance_all <- attr(pca_data_all, "percentVar")[1:2]

if (is.null(variance_real)) {
    # 如果没有percentVar属性，手动计算
    pca_real_full <- prcomp(t(counts_real), scale. = TRUE)
    variance_real <- summary(pca_real_full)$importance[2, 1:2] * 100
    
    pca_gen_full <- prcomp(t(counts_gen), scale. = TRUE)
    variance_gen <- summary(pca_gen_full)$importance[2, 1:2] * 100
    
    pca_all_full <- prcomp(t(counts_all), scale. = TRUE)
    variance_all <- summary(pca_all_full)$importance[2, 1:2] * 100
}

cat(sprintf("  真实数据 - PC1: %.1f%%, PC2: %.1f%%, 累计: %.1f%%\n", 
            variance_real[1], variance_real[2], sum(variance_real)))
cat(sprintf("  生成数据 - PC1: %.1f%%, PC2: %.1f%%, 累计: %.1f%%\n", 
            variance_gen[1], variance_gen[2], sum(variance_gen)))
cat(sprintf("  合并数据 - PC1: %.1f%%, PC2: %.1f%%, 累计: %.1f%%\n", 
            variance_all[1], variance_all[2], sum(variance_all)))

# 2. 聚类分离度评估（轮廓系数）
library(cluster)

calculate_silhouette <- function(pca_data, labels) {
    if (length(unique(labels)) > 1) {
        dist_matrix <- dist(pca_data[, 1:2])
        sil <- silhouette(as.numeric(as.factor(labels)), dist_matrix)
        return(mean(sil[, 3]))
    } else {
        return(NA)
    }
}

cat("\n🎯 聚类质量评估（轮廓系数）:\n")
cat("  范围: -1到1，越接近1越好，>0.5为良好\n")

sil_real <- calculate_silhouette(pca_data_real, org_label[1:ncol(counts_real)])
sil_gen <- calculate_silhouette(pca_data_gen, gen_label[1:ncol(counts_gen)])
all_labels <- c(org_label, gen_label)[1:ncol(counts_all)]
sil_all <- calculate_silhouette(pca_data_all, all_labels)

cat(sprintf("  真实数据轮廓系数: %.3f %s\n", sil_real, 
            ifelse(sil_real > 0.5, "✅ 优秀", ifelse(sil_real > 0.25, "⚠️  一般", "❌ 较差"))))
cat(sprintf("  生成数据轮廓系数: %.3f %s\n", sil_gen, 
            ifelse(sil_gen > 0.5, "✅ 优秀", ifelse(sil_gen > 0.25, "⚠️  一般", "❌ 较差"))))
cat(sprintf("  合并数据轮廓系数: %.3f %s\n", sil_all, 
            ifelse(sil_all > 0.5, "✅ 优秀", ifelse(sil_all > 0.25, "⚠️  一般", "❌ 较差"))))

# 3. 分布相似性评估
cat("\n📏 分布相似性评估:\n")

# 计算每个细胞类型的中心点
cell_types <- unique(org_label)
real_centers <- matrix(0, nrow = length(cell_types), ncol = 2)
gen_centers <- matrix(0, nrow = length(cell_types), ncol = 2)

for (i in seq_along(cell_types)) {
    real_idx <- which(org_label[1:ncol(counts_real)] == cell_types[i])
    gen_idx <- which(gen_label[1:ncol(counts_gen)] == cell_types[i])
    
    if (length(real_idx) > 0) {
        real_centers[i, ] <- colMeans(pca_data_real[real_idx, 1:2, drop = FALSE])
    }
    if (length(gen_idx) > 0) {
        gen_centers[i, ] <- colMeans(pca_data_gen[gen_idx, 1:2, drop = FALSE])
    }
}

# 计算中心点距离
center_distances <- sqrt(rowSums((real_centers - gen_centers)^2))
names(center_distances) <- cell_types

cat("  聚类中心距离（越小越好）:\n")
for (i in seq_along(center_distances)) {
    cat(sprintf("    %s: %.2f\n", names(center_distances)[i], center_distances[i]))
}

avg_center_distance <- mean(center_distances, na.rm = TRUE)
cat(sprintf("  平均中心距离: %.2f %s\n", avg_center_distance,
            ifelse(avg_center_distance < 5, "✅ 优秀", ifelse(avg_center_distance < 10, "⚠️  一般", "❌ 较差"))))

# 4. 数据重叠度评估
cat("\n🔄 真实与生成数据重叠度:\n")

# 计算重叠区域
real_range_pc1 <- range(pca_data_real[, 1])
real_range_pc2 <- range(pca_data_real[, 2])
gen_range_pc1 <- range(pca_data_gen[, 1])
gen_range_pc2 <- range(pca_data_gen[, 2])

overlap_pc1 <- max(0, min(real_range_pc1[2], gen_range_pc1[2]) - max(real_range_pc1[1], gen_range_pc1[1])) /
               max(real_range_pc1[2] - real_range_pc1[1], gen_range_pc1[2] - gen_range_pc1[1])

overlap_pc2 <- max(0, min(real_range_pc2[2], gen_range_pc2[2]) - max(real_range_pc2[1], gen_range_pc2[1])) /
               max(real_range_pc2[2] - real_range_pc2[1], gen_range_pc2[2] - gen_range_pc2[1])

avg_overlap <- (overlap_pc1 + overlap_pc2) / 2

cat(sprintf("  PC1维度重叠度: %.1f%%\n", overlap_pc1 * 100))
cat(sprintf("  PC2维度重叠度: %.1f%%\n", overlap_pc2 * 100))
cat(sprintf("  平均重叠度: %.1f%% %s\n", avg_overlap * 100,
            ifelse(avg_overlap > 0.8, "✅ 优秀", ifelse(avg_overlap > 0.6, "⚠️  一般", "❌ 较差"))))

# ========== 综合评分 ==========
cat("\n🏆 综合质量评分:\n")

score_silhouette <- pmax(0, pmin(1, (sil_real + sil_gen) / 2))
score_center_dist <- pmax(0, pmin(1, 1 - avg_center_distance / 20))
score_overlap <- avg_overlap
score_variance <- pmin(1, sum(variance_real) / 50)  # 假设50%为满分

total_score <- (score_silhouette + score_center_dist + score_overlap + score_variance) / 4

cat(sprintf("  聚类质量得分: %.2f/1.0\n", score_silhouette))
cat(sprintf("  中心距离得分: %.2f/1.0\n", score_center_dist))
cat(sprintf("  数据重叠得分: %.2f/1.0\n", score_overlap))
cat(sprintf("  方差解释得分: %.2f/1.0\n", score_variance))
cat("  ", paste(rep("-", 25), collapse = ""), "\n")
cat(sprintf("  📊 总体质量得分: %.2f/1.0 %s\n", total_score,
            ifelse(total_score > 0.8, "🌟 优秀", ifelse(total_score > 0.6, "👍 良好", ifelse(total_score > 0.4, "⚠️  一般", "❌ 需改进")))))

# ========== 保存评估结果 ==========
assessment_results <- data.frame(
    Metric = c("Real_Silhouette", "Generated_Silhouette", "Combined_Silhouette",
               "Average_Center_Distance", "PC1_Overlap", "PC2_Overlap", "Average_Overlap",
               "Real_PC1_Variance", "Real_PC2_Variance", "Total_Score"),
    Value = c(sil_real, sil_gen, sil_all, avg_center_distance, 
              overlap_pc1, overlap_pc2, avg_overlap,
              variance_real[1], variance_real[2], total_score)
)

write.csv(assessment_results, get_output_path("quality_assessment", "pca_quality_assessment.csv"), row.names = FALSE)

cat("\n💾 评估结果已保存到:", get_output_path("quality_assessment", "pca_quality_assessment.csv"), "\n")
cat("==========================================\n")
cat("✅ PCA质量评估完成！\n") 