# KEGG富集分析 - 详细比较分析
# 用于深入分析真实数据和生成数据的差异

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(enrichplot)
library(VennDiagram)
library(pheatmap)

# 加载目录结构创建函数
source("code/create_output_structure.R")

# 获取数据集路径信息
dataset_paths <- get_dataset_paths()
dataset_name <- dataset_paths$dataset_name

# 创建输出目录结构
output_paths <- create_output_structure()

# ========== 1. 智能加载KEGG结果 ==========
# 优先级: 文件 > 内存 > 重新计算

kegg_real_file <- get_output_path("kegg_analysis", "kegg_real_results.rds")
kegg_gen_file <- get_output_path("kegg_analysis", "kegg_generated_results.rds")

# 检查是否有保存的文件
if (file.exists(kegg_real_file) && file.exists(kegg_gen_file)) {
  cat("📂 从文件加载KEGG分析结果...\n")
  kegg_real <- readRDS(kegg_real_file)
  kegg_gen <- readRDS(kegg_gen_file)
  cat("✅ KEGG结果加载成功!\n")
  cat("  - 真实数据结果: ", kegg_real_file, "\n")
  cat("  - 生成数据结果: ", kegg_gen_file, "\n")
} else if (exists("kegg_real") && exists("kegg_gen")) {
  cat("💾 从内存中获取KEGG结果...\n")
  cat("✅ 使用内存中的KEGG结果\n")
} else {
  cat("⚠️  未找到KEGG结果文件，正在重新运行分析...\n")
  source("code/kegg.R")
  cat("✅ KEGG分析完成并已保存!\n")
}

# ========== 2. 详细比较分析 ==========

# 提取通路信息
real_pathways <- kegg_real@result
gen_pathways <- kegg_gen@result

cat("📊 详细统计信息:\n")
cat("真实数据富集通路数:", nrow(real_pathways), "\n")
cat("生成数据富集通路数:", nrow(gen_pathways), "\n")

# 获取通路名称
real_terms <- real_pathways$Description
gen_terms <- gen_pathways$Description

# 计算差异
common_terms <- intersect(real_terms, gen_terms)
only_real <- setdiff(real_terms, gen_terms)
only_gen <- setdiff(gen_terms, real_terms)

cat("共同通路数:", length(common_terms), "\n")
cat("仅在真实数据中的通路数:", length(only_real), "\n")
cat("仅在生成数据中的通路数:", length(only_gen), "\n")

# ========== 3. 保存详细比较结果 ==========

# 保存仅在真实数据中的通路
if (length(only_real) > 0) {
  write.csv(data.frame(Real_Only_Pathways = only_real), 
            get_output_path("kegg_analysis", "real_only_pathways.csv"), row.names = FALSE)
  cat("💾 仅在真实数据中的通路已保存到:", get_output_path("kegg_analysis", "real_only_pathways.csv"), "\n")
}

# 保存仅在生成数据中的通路  
if (length(only_gen) > 0) {
  write.csv(data.frame(Generated_Only_Pathways = only_gen), 
            get_output_path("kegg_analysis", "generated_only_pathways.csv"), row.names = FALSE)
  cat("💾 仅在生成数据中的通路已保存到:", get_output_path("kegg_analysis", "generated_only_pathways.csv"), "\n")
}

# ========== 4. 通路富集强度比较 ==========

# 合并通路信息进行比较
real_pathways$Data_Type <- "Real"
gen_pathways$Data_Type <- "Generated"

# 只比较共同通路
common_real <- real_pathways[real_pathways$Description %in% common_terms, ]
common_gen <- gen_pathways[gen_pathways$Description %in% common_terms, ]

# 创建比较数据框
comparison_df <- data.frame(
  Pathway = common_real$Description,
  Real_pvalue = -log10(common_real$pvalue),
  Generated_pvalue = -log10(common_gen$pvalue[match(common_real$Description, common_gen$Description)]),
  Real_Count = common_real$Count,
  Generated_Count = common_gen$Count[match(common_real$Description, common_gen$Description)]
)

# 保存比较结果
write.csv(comparison_df, get_output_path("kegg_analysis", "pathway_enrichment_comparison.csv"), row.names = FALSE)

# ========== 5. 创建比较图表 ==========

pdf(get_output_path("kegg_analysis", "detailed_kegg_comparison.pdf"), width = 16, height = 10)

# 图1: 富集强度散点图
p1 <- ggplot(comparison_df, aes(x = Real_pvalue, y = Generated_pvalue)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Real Data -log10(p-value)", 
       y = "Generated Data -log10(p-value)",
       title = "KEGG通路富集强度比较") +
  theme_minimal() +
  coord_equal()

print(p1)

# 图2: 基因数量比较
p2 <- ggplot(comparison_df, aes(x = Real_Count, y = Generated_Count)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Real Data Gene Count", 
       y = "Generated Data Gene Count",
       title = "KEGG通路中基因数量比较") +
  theme_minimal()

print(p2)

# 图3: 顶级通路比较柱状图
top_pathways <- head(comparison_df[order(-comparison_df$Real_pvalue), ], 15)
top_long <- reshape2::melt(top_pathways[, c("Pathway", "Real_pvalue", "Generated_pvalue")], 
                          id.vars = "Pathway")

p3 <- ggplot(top_long, aes(x = reorder(Pathway, value), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "KEGG Pathway", 
       y = "-log10(p-value)",
       title = "顶级KEGG通路富集强度比较",
       fill = "数据类型") +
  scale_fill_manual(values = c("Real_pvalue" = "skyblue", "Generated_pvalue" = "orange"),
                   labels = c("真实数据", "生成数据")) +
  theme_minimal()

print(p3)

# 图4: 神经系统相关通路特写
neuro_keywords <- c("neuro", "synap", "axon", "alzheimer", "parkinson", 
                   "huntington", "amyotrophic", "dopamin", "serotonin", "gaba")

neuro_pathways <- comparison_df[
  grepl(paste(neuro_keywords, collapse = "|"), 
        comparison_df$Pathway, ignore.case = TRUE), ]

if (nrow(neuro_pathways) > 0) {
  neuro_long <- reshape2::melt(neuro_pathways[, c("Pathway", "Real_pvalue", "Generated_pvalue")], 
                              id.vars = "Pathway")
  
  p4 <- ggplot(neuro_long, aes(x = reorder(Pathway, value), y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(x = "神经系统相关通路", 
         y = "-log10(p-value)",
         title = "神经系统相关KEGG通路比较",
         fill = "数据类型") +
    scale_fill_manual(values = c("Real_pvalue" = "lightgreen", "Generated_pvalue" = "salmon"),
                     labels = c("真实数据", "生成数据")) +
    theme_minimal()
  
  print(p4)
}

dev.off()

# ========== 6. 生成总结统计 ==========

summary_stats <- data.frame(
  Metric = c("总通路数_真实", "总通路数_生成", "共同通路数", 
             "仅真实数据", "仅生成数据", "相似度百分比"),
  Value = c(nrow(real_pathways), nrow(gen_pathways), length(common_terms),
            length(only_real), length(only_gen), 
            round(length(common_terms) / max(nrow(real_pathways), nrow(gen_pathways)) * 100, 2))
)

write.csv(summary_stats, get_output_path("kegg_analysis", "kegg_summary_statistics.csv"), row.names = FALSE)

# ========== 7. 输出重要发现 ==========
cat("\n========================================\n")
cat("🎯 关键发现总结:\n")  
cat("========================================\n")

cat("📊 数据相似度:", round(length(common_terms) / max(nrow(real_pathways), nrow(gen_pathways)) * 100, 2), "%\n")

if (length(only_real) > 0) {
  cat("\n🔍 仅在真实数据中发现的通路 (前5个):\n")
  cat(paste("  -", head(only_real, 5), collapse = "\n"), "\n")
}

if (length(only_gen) > 0) {
  cat("\n🔍 仅在生成数据中发现的通路 (前5个):\n")
  cat(paste("  -", head(only_gen, 5), collapse = "\n"), "\n")
}

# 找出富集差异最大的通路
comparison_df$pvalue_diff <- abs(comparison_df$Real_pvalue - comparison_df$Generated_pvalue)
top_diff <- head(comparison_df[order(-comparison_df$pvalue_diff), ], 5)

cat("\n📈 富集强度差异最大的通路:\n")
for (i in 1:nrow(top_diff)) {
  cat(sprintf("  - %s (差异: %.2f)\n", 
              top_diff$Pathway[i], top_diff$pvalue_diff[i]))
}

cat("\n✅ 详细分析完成! 所有结果已保存到", output_paths$kegg_analysis, "目录\n")
cat("📁 生成的文件:\n")
cat("  -", get_output_path("kegg_analysis", "detailed_kegg_comparison.pdf"), "(详细比较图表)\n")
cat("  -", get_output_path("kegg_analysis", "pathway_enrichment_comparison.csv"), "(富集强度比较)\n") 
cat("  -", get_output_path("kegg_analysis", "kegg_summary_statistics.csv"), "(总结统计)\n")
if (length(only_real) > 0) cat("  -", get_output_path("kegg_analysis", "real_only_pathways.csv"), "(仅真实数据通路)\n")
if (length(only_gen) > 0) cat("  -", get_output_path("kegg_analysis", "generated_only_pathways.csv"), "(仅生成数据通路)\n") 