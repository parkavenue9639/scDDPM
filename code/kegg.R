# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(enrichplot)
library(VennDiagram)

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

if (!file.exists(real_data_file)) {
    cat("❌ 真实数据文件不存在:", real_data_file, "\n")
    cat("请先运行 Preprocess.R 生成预处理数据\n")
    stop("文件不存在")
}

if (!file.exists(gen_data_file)) {
    cat("❌ 生成数据文件不存在:", gen_data_file, "\n")
    cat("请先运行 scDDPM.py 生成数据\n")
    stop("文件不存在")
}

# ========== Load Data ==========
# Load real and generated expression data
expr_real <- read.csv(real_data_file, check.names = FALSE)
rownames(expr_real) <- make.unique(as.character(expr_real[, 1]))
expr_real <- expr_real[, -1]

expr_gen <- read.csv(gen_data_file, check.names = FALSE)
rownames(expr_gen) <- make.unique(as.character(expr_gen[, 1]))
expr_gen <- expr_gen[, -1]

# ========== Extract Labels ==========
labels_real <- as.character(expr_real$label)
expr_real$label <- NULL
expr_real <- t(expr_real)

labels_gen <- as.character(expr_gen$label)
expr_gen$label <- NULL
expr_gen <- t(expr_gen)

# 清理可能的引号
labels_real <- gsub('^"|"$', '', labels_real)
labels_gen <- gsub('^"|"$', '', labels_gen)

# 显示实际的标签值用于调试
cat("🔍 真实数据中的唯一标签:", paste(unique(labels_real), collapse = ", "), "\n")
cat("🔍 生成数据中的唯一标签:", paste(unique(labels_gen), collapse = ", "), "\n")

# ========== 数据集细胞类型动态分析 ==========
# 动态获取和分析数据集的细胞类型组成
unique_labels_real <- unique(labels_real)
unique_labels_gen <- unique(labels_gen)

# 找到两个数据集共同的标签
common_labels <- intersect(unique_labels_real, unique_labels_gen)

# 分析数据集基本信息
cat("📊 ========== 数据集细胞类型组成分析 ==========\n")
cat("🔬 真实数据包含", length(unique_labels_real), "种细胞类型:", paste(unique_labels_real, collapse = ", "), "\n")
cat("🧬 生成数据包含", length(unique_labels_gen), "种细胞类型:", paste(unique_labels_gen, collapse = ", "), "\n")
cat("🔄 共同细胞类型", length(common_labels), "种:", paste(common_labels, collapse = ", "), "\n")

# 显示每种细胞类型的样本数分布
cat("\n📈 样本数分布:\n")
for (label in common_labels) {
    real_count <- sum(labels_real == label)
    gen_count <- sum(labels_gen == label)
    cat("   -", label, ": 真实数据", real_count, "个, 生成数据", gen_count, "个\n")
}

# 动态验证数据集是否适合进行分析
if (length(common_labels) < 2) {
    cat("❌ 数据集验证失败: 共同标签少于2个，无法进行差异表达分析\n")
    cat("建议检查数据标签是否一致或数据预处理是否正确\n")
    stop("需要至少2个共同的细胞类型标签来进行比较分析")
}

# 计算总样本数和数据集规模
total_real_samples <- length(labels_real)
total_gen_samples <- length(labels_gen)
cat("\n📊 数据集规模: 真实数据", total_real_samples, "个样本, 生成数据", total_gen_samples, "个样本\n")

# 根据数据集规模动态调整分析参数
min_samples_per_group <- max(3, ceiling(total_real_samples / (length(common_labels) * 50)))  # 动态最小样本数
min_deg_threshold <- ifelse(length(common_labels) >= 5, 0.3, 0.5)  # 动态logFC阈值
min_genes_for_kegg <- max(5, ceiling(length(common_labels)))  # 动态KEGG分析基因数阈值

cat("🔧 动态分析参数:\n")
cat("   - 最小组样本数:", min_samples_per_group, "\n")
cat("   - 差异基因logFC阈值:", min_deg_threshold, "\n")
cat("   - KEGG分析最小基因数:", min_genes_for_kegg, "\n")
cat("========================================================\n\n")

# ========== 动态多细胞类型差异表达分析函数 ==========
# 每种细胞类型 vs 其他所有细胞类型，使用动态参数
get_deg_one_vs_rest <- function(expr_mat, labels, target_celltype, min_samples = min_samples_per_group, logfc_threshold = min_deg_threshold) {
    target_idx <- which(labels == target_celltype)
    other_idx <- which(labels != target_celltype)
    
    # 动态样本数验证
    if (length(target_idx) < min_samples || length(other_idx) < min_samples) {
        cat("⚠️  细胞类型", target_celltype, "样本数不足 (", length(target_idx), "vs", length(other_idx), ", 需要≥", min_samples, ")\n")
        return(NULL)
    }
    
    cat("🔍 分析", target_celltype, ": ", length(target_idx), "个目标样本 vs", length(other_idx), "个对照样本...\n")
    
    # 差异表达分析
    pvals <- apply(expr_mat, 1, function(x) {
        tryCatch(wilcox.test(x[target_idx], x[other_idx])$p.value, error = function(e) NA)
    })
    logfc <- rowMeans(expr_mat[, target_idx, drop = FALSE]) - rowMeans(expr_mat[, other_idx, drop = FALSE])
    
    deg_df <- data.frame(
        gene = rownames(expr_mat), 
        logFC = logfc, 
        pval = pvals,
        celltype = target_celltype,
        stringsAsFactors = FALSE
    )
    
    # 使用动态阈值过滤显著差异基因
    deg_df <- deg_df[!is.na(deg_df$pval) & deg_df$pval < 0.05 & abs(deg_df$logFC) > logfc_threshold, ]
    
    cat("   ✅", target_celltype, ": 发现", nrow(deg_df), "个差异表达基因 (logFC >", logfc_threshold, ")\n")
    return(deg_df)
}

# 智能分析进度跟踪函数
analyze_celltype_progress <- function(current, total, celltype) {
    progress_percent <- round((current / total) * 100, 1)
    progress_bar <- paste0("[", paste(rep("=", floor(progress_percent / 5)), collapse = ""), 
                          paste(rep("-", 20 - floor(progress_percent / 5)), collapse = ""), "]")
    cat("🚀 进度", progress_bar, progress_percent, "% - 当前:", celltype, "\n")
}

# ========== 智能多细胞类型差异表达分析主流程 ==========
cat("🧬 开始", length(common_labels), "种细胞类型的差异表达分析...\n\n")

deg_real_list <- list()
deg_gen_list <- list()

# 使用进度跟踪的分析循环
for (i in seq_along(common_labels)) {
    celltype <- common_labels[i]
    
    # 显示分析进度
    analyze_celltype_progress(i, length(common_labels), celltype)
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    
    # 真实数据分析
    cat("📊 [真实数据] 分析", celltype, "差异表达基因...\n")
    deg_real_ct <- get_deg_one_vs_rest(expr_real, labels_real, celltype, min_samples_per_group, min_deg_threshold)
    if (!is.null(deg_real_ct) && nrow(deg_real_ct) > 0) {
        deg_real_list[[celltype]] <- deg_real_ct
        cat("   ✅ 真实数据:", celltype, "成功识别", nrow(deg_real_ct), "个差异基因\n")
    } else {
        cat("   ⚠️  真实数据:", celltype, "未发现足够的差异基因\n")
    }
    
    # 生成数据分析
    cat("🧬 [生成数据] 分析", celltype, "差异表达基因...\n")
    deg_gen_ct <- get_deg_one_vs_rest(expr_gen, labels_gen, celltype, min_samples_per_group, min_deg_threshold)
    if (!is.null(deg_gen_ct) && nrow(deg_gen_ct) > 0) {
        deg_gen_list[[celltype]] <- deg_gen_ct
        cat("   ✅ 生成数据:", celltype, "成功识别", nrow(deg_gen_ct), "个差异基因\n")
    } else {
        cat("   ⚠️  生成数据:", celltype, "未发现足够的差异基因\n")
    }
    
    cat("\n")
}

# 智能分析结果验证
successful_real <- length(deg_real_list)
successful_gen <- length(deg_gen_list)
successful_both <- length(intersect(names(deg_real_list), names(deg_gen_list)))

cat("📈 ========== 差异表达分析结果汇总 ==========\n")
cat("🔬 成功分析的细胞类型 (真实数据):", successful_real, "/", length(common_labels), "\n")
cat("🧬 成功分析的细胞类型 (生成数据):", successful_gen, "/", length(common_labels), "\n")
cat("🎯 两类数据都成功的细胞类型:", successful_both, "/", length(common_labels), "\n")

if (successful_both == 0) {
    cat("❌ 数据质量检查失败: 没有细胞类型在两类数据中都发现足够的差异基因\n")
    cat("建议:\n")
    cat("   1. 检查数据质量和标签一致性\n")
    cat("   2. 降低差异基因阈值 (当前logFC >", min_deg_threshold, ")\n")
    cat("   3. 增加样本数或检查数据预处理\n")
    stop("无法进行后续KEGG分析")
}

cat("✅ 差异表达分析完成，将对", successful_both, "种细胞类型进行KEGG富集分析\n")
cat("============================================\n\n")

# ========== 多细胞类型基因ID转换和KEGG富集分析 ==========
# 对每种细胞类型进行基因ID转换和KEGG分析
kegg_real_list <- list()
kegg_gen_list <- list()
gene_real_list <- list()
gene_gen_list <- list()

common_celltypes <- intersect(names(deg_real_list), names(deg_gen_list))
cat("🔬 开始对", length(common_celltypes), "种细胞类型进行KEGG富集分析\n")

for (celltype in common_celltypes) {
    cat("🧬 处理细胞类型:", celltype, "\n")
    
    # 真实数据的基因ID转换和KEGG分析
    if (celltype %in% names(deg_real_list)) {
        genes_real <- deg_real_list[[celltype]]$gene
        gene_real_ct <- bitr(genes_real, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        gene_real_list[[celltype]] <- gene_real_ct
        
        cat("   📊 真实数据: ", length(genes_real), "个DEG ->", nrow(gene_real_ct), "个ENTREZ ID\n")
        
        if (nrow(gene_real_ct) >= min_genes_for_kegg) {  # 动态基因数阈值
            kegg_real_ct <- enrichKEGG(gene = gene_real_ct$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
            if (!is.null(kegg_real_ct) && nrow(kegg_real_ct@result) > 0) {
                kegg_real_list[[celltype]] <- kegg_real_ct
                cat("   ✅ 真实数据KEGG通路:", nrow(kegg_real_ct@result), "条\n")
            } else {
                cat("   ⚠️  真实数据未发现显著KEGG通路\n")
            }
        } else {
            cat("   ⚠️  真实数据基因数太少，跳过KEGG分析\n")
        }
    }
    
    # 生成数据的基因ID转换和KEGG分析
    if (celltype %in% names(deg_gen_list)) {
        genes_gen <- deg_gen_list[[celltype]]$gene
        gene_gen_ct <- bitr(genes_gen, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        gene_gen_list[[celltype]] <- gene_gen_ct
        
        cat("   📊 生成数据: ", length(genes_gen), "个DEG ->", nrow(gene_gen_ct), "个ENTREZ ID\n")
        
        if (nrow(gene_gen_ct) >= min_genes_for_kegg) {  # 动态基因数阈值
            kegg_gen_ct <- enrichKEGG(gene = gene_gen_ct$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
            if (!is.null(kegg_gen_ct) && nrow(kegg_gen_ct@result) > 0) {
                kegg_gen_list[[celltype]] <- kegg_gen_ct
                cat("   ✅ 生成数据KEGG通路:", nrow(kegg_gen_ct@result), "条\n")
            } else {
                cat("   ⚠️  生成数据未发现显著KEGG通路\n")
            }
        } else {
            cat("   ⚠️  生成数据基因数太少，跳过KEGG分析\n")
        }
    }
}

# 检查KEGG分析结果
if (length(kegg_real_list) == 0 && length(kegg_gen_list) == 0) {
    stop("❌ 所有细胞类型都未能进行KEGG富集分析")
}

# ========== 保存多细胞类型KEGG结果到文件 ==========
# 保存差异基因列表
saveRDS(deg_real_list, get_output_path("kegg_analysis", "deg_real_by_celltype.rds"))
saveRDS(deg_gen_list, get_output_path("kegg_analysis", "deg_gen_by_celltype.rds"))

# 保存KEGG结果
saveRDS(kegg_real_list, get_output_path("kegg_analysis", "kegg_real_by_celltype.rds"))
saveRDS(kegg_gen_list, get_output_path("kegg_analysis", "kegg_gen_by_celltype.rds"))

cat("💾 多细胞类型KEGG分析结果已保存到:\n")
cat("  - 差异基因:", get_output_path("kegg_analysis", "deg_real_by_celltype.rds"), "\n")
cat("  - 差异基因:", get_output_path("kegg_analysis", "deg_gen_by_celltype.rds"), "\n")
cat("  - KEGG结果:", get_output_path("kegg_analysis", "kegg_real_by_celltype.rds"), "\n")
cat("  - KEGG结果:", get_output_path("kegg_analysis", "kegg_gen_by_celltype.rds"), "\n")

# ========== 多细胞类型通路比较分析 ==========
# 汇总所有细胞类型的通路比较结果
pathway_comparison_summary <- list()
all_common_pathways <- list()

# 创建汇总报告
cat("📋 生成多细胞类型KEGG通路比较报告...\n")

for (celltype in common_celltypes) {
    if (celltype %in% names(kegg_real_list) && celltype %in% names(kegg_gen_list)) {
        real_terms <- kegg_real_list[[celltype]]@result$Description
        gen_terms  <- kegg_gen_list[[celltype]]@result$Description
        
        common_terms <- intersect(real_terms, gen_terms)
        only_real    <- setdiff(real_terms, gen_terms)
        only_gen     <- setdiff(gen_terms, real_terms)
        
        pathway_comparison_summary[[celltype]] <- list(
            celltype = celltype,
            real_pathways = length(real_terms),
            gen_pathways = length(gen_terms),
            common_pathways = length(common_terms),
            only_real = length(only_real),
            only_gen = length(only_gen),
            consistency_rate = ifelse(length(real_terms) + length(gen_terms) > 0, 
                                    length(common_terms) / (length(real_terms) + length(gen_terms) - length(common_terms)) * 100, 0)
        )
        
        all_common_pathways[[celltype]] <- common_terms
        
        cat("🧬", celltype, ": 真实", length(real_terms), "条, 生成", length(gen_terms), "条, 共同", length(common_terms), "条\n")
    }
}

# 保存通路比较结果
write.csv(do.call(rbind, lapply(pathway_comparison_summary, function(x) data.frame(x))), 
          get_output_path("kegg_analysis", "pathway_comparison_summary.csv"), row.names = FALSE)

# 保存每个细胞类型的共同通路
for (celltype in names(all_common_pathways)) {
    if (length(all_common_pathways[[celltype]]) > 0) {
        write.csv(data.frame(Common_Pathways = all_common_pathways[[celltype]]), 
                  get_output_path("kegg_analysis", paste0("common_pathways_", gsub(" ", "_", celltype), ".csv")), 
                  row.names = FALSE)
    }
}

# ========== 智能多细胞类型可视化 ==========
cat("🎨 生成", length(common_celltypes), "种细胞类型的KEGG可视化图表...\n")

# 智能计算图表布局和尺寸
n_real <- length(kegg_real_list)
n_gen <- length(kegg_gen_list)
n_total <- n_real + n_gen

cat("📊 可视化统计: 真实数据", n_real, "个图表, 生成数据", n_gen, "个图表, 总计", n_total, "个\n")

if (n_total > 0) {
    # 智能布局算法：根据细胞类型数量优化布局
    if (n_total <= 2) {
        n_cols <- n_total; n_rows <- 1
    } else if (n_total <= 4) {
        n_cols <- 2; n_rows <- ceiling(n_total / 2)
    } else if (n_total <= 9) {
        n_cols <- 3; n_rows <- ceiling(n_total / 3)
    } else {
        n_cols <- 4; n_rows <- ceiling(n_total / 4)
    }
    
    # 动态调整图表尺寸
    chart_width <- max(8, n_cols * 5)
    chart_height <- max(6, n_rows * 4)
    
    cat("📐 图表布局:", n_rows, "行 x", n_cols, "列, 尺寸:", chart_width, "x", chart_height, "\n")
    
    pdf(get_output_path("kegg_analysis", "KEGG_multicelltype_analysis.pdf"), width = chart_width, height = chart_height)
    
    # 设置多面板布局
    par(mfrow = c(n_rows, n_cols))
    
    # 绘制每个细胞类型的真实数据KEGG图
    for (celltype in names(kegg_real_list)) {
        if (nrow(kegg_real_list[[celltype]]@result) > 0) {
            tryCatch({
                print(dotplot(kegg_real_list[[celltype]], showCategory = 10, 
                             title = paste("Real Data -", celltype)))
            }, error = function(e) {
                cat("⚠️  绘制", celltype, "真实数据图表时出错:", e$message, "\n")
            })
        }
    }
    
    # 绘制每个细胞类型的生成数据KEGG图
    for (celltype in names(kegg_gen_list)) {
        if (nrow(kegg_gen_list[[celltype]]@result) > 0) {
            tryCatch({
                print(dotplot(kegg_gen_list[[celltype]], showCategory = 10, 
                             title = paste("Generated Data -", celltype)))
            }, error = function(e) {
                cat("⚠️  绘制", celltype, "生成数据图表时出错:", e$message, "\n")
            })
        }
    }
    
    dev.off()
    
    # 生成通路一致性汇总图
    if (length(pathway_comparison_summary) > 0) {
        pdf(get_output_path("kegg_analysis", "KEGG_consistency_summary.pdf"), width = 12, height = 8)
        
        # 准备数据
        summary_df <- do.call(rbind, lapply(pathway_comparison_summary, function(x) data.frame(x)))
        
        # 一致性率柱状图
        if (nrow(summary_df) > 0) {
            barplot(summary_df$consistency_rate, 
                   names.arg = summary_df$celltype,
                   main = "KEGG Pathway Consistency Rate by Cell Type",
                   ylab = "Consistency Rate (%)",
                   xlab = "Cell Type",
                   col = rainbow(nrow(summary_df)),
                   las = 2)
            
            # 添加数值标签
            text(1:nrow(summary_df) - 0.5, summary_df$consistency_rate + 2, 
                 paste0(round(summary_df$consistency_rate, 1), "%"), 
                 pos = 3, cex = 0.8)
        }
        
        dev.off()
    }
    
    cat("\n🎉 ========== 动态多细胞类型KEGG分析完成 ==========\n")
    cat("📊 最终分析汇总:\n")
    cat("   🔬 数据集规模: 真实", total_real_samples, "样本, 生成", total_gen_samples, "样本\n")
    cat("   🧬 细胞类型总数:", length(common_celltypes), "种\n")
    cat("   ✅ 成功KEGG分析:", length(pathway_comparison_summary), "种细胞类型\n")
    cat("   📈 分析成功率:", round(length(pathway_comparison_summary) / length(common_celltypes) * 100, 1), "%\n")
    
    if (length(pathway_comparison_summary) > 0) {
        avg_consistency <- mean(sapply(pathway_comparison_summary, function(x) x$consistency_rate))
        max_consistency <- max(sapply(pathway_comparison_summary, function(x) x$consistency_rate))
        min_consistency <- min(sapply(pathway_comparison_summary, function(x) x$consistency_rate))
        
        cat("   🎯 通路一致性:\n")
        cat("      - 平均:", round(avg_consistency, 1), "%\n")
        cat("      - 最高:", round(max_consistency, 1), "%\n")
        cat("      - 最低:", round(min_consistency, 1), "%\n")
        
        # 找出表现最好和最差的细胞类型
        best_celltype <- names(which.max(sapply(pathway_comparison_summary, function(x) x$consistency_rate)))
        worst_celltype <- names(which.min(sapply(pathway_comparison_summary, function(x) x$consistency_rate)))
        
        cat("   🏆 最佳一致性:", best_celltype, "(", round(max_consistency, 1), "%)\n")
        cat("   ⚠️  最低一致性:", worst_celltype, "(", round(min_consistency, 1), "%)\n")
    }
    
    # 显示使用的动态参数
    cat("\n🔧 本次分析使用的动态参数:\n")
    cat("   - 最小组样本数:", min_samples_per_group, "\n")
    cat("   - 差异基因logFC阈值:", min_deg_threshold, "\n")
    cat("   - KEGG分析最小基因数:", min_genes_for_kegg, "\n")
    
    cat("===============================================\n")
    
} else {
    cat("⚠️  没有足够的KEGG结果进行可视化\n")
    cat("💡 建议降低分析阈值或检查数据质量\n")
}

cat("🎯 动态多细胞类型KEGG分析完成！\n")
cat("📁 所有结果已保存到:", output_paths$kegg_analysis, "目录\n")
cat("📊 此分析框架可自动适配任意数量的细胞类型 (2-50+ 种)\n")
