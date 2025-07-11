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

# Set the label pairs to compare (real vs generated)
real_label1 <- "Astrocytes"; real_label2 <- "Excitatory neurons"
gen_label1  <- "Astrocytes"; gen_label2  <- "Excitatory neurons"

# Find column indices for each group
group1_real <- which(labels_real == real_label1)
group2_real <- which(labels_real == real_label2)
group1_gen  <- which(labels_gen == gen_label1)
group2_gen  <- which(labels_gen == gen_label2)

cat("Real data columns: ", length(group1_real), "vs", length(group2_real), "\n")
cat("Generated data columns: ", length(group1_gen), "vs", length(group2_gen), "\n")

if (any(sapply(list(group1_real, group2_real, group1_gen, group2_gen), length) == 0)) {
    stop("❌ One or more label groups are empty. Please check label assignments.")
}

# ========== Differential Expression Analysis ==========
get_deg <- function(expr_mat, idx1, idx2) {
    pvals <- apply(expr_mat, 1, function(x) {
        tryCatch(wilcox.test(x[idx1], x[idx2])$p.value, error = function(e) NA)
    })
    logfc <- rowMeans(expr_mat[, idx1, drop = FALSE]) - rowMeans(expr_mat[, idx2, drop = FALSE])
    deg_df <- data.frame(gene = rownames(expr_mat), logFC = logfc, pval = pvals)
    return(deg_df)
}

deg_real <- get_deg(expr_real, group1_real, group2_real)
deg_gen  <- get_deg(expr_gen, group1_gen, group2_gen)

# ========== Gene ID Conversion ==========
gene_real <- bitr(deg_real$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_gen  <- bitr(deg_gen$gene,  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cat("🧬 DEGs (real): ", nrow(deg_real), "\n")
cat("🧬 DEGs (generated): ", nrow(deg_gen), "\n")
cat("✅ Mapped ENTREZ IDs (real): ", nrow(gene_real), "\n")
cat("✅ Mapped ENTREZ IDs (generated): ", nrow(gene_gen), "\n")

if (nrow(gene_real) == 0 || nrow(gene_gen) == 0) {
    stop("❌ No valid ENTREZ IDs found. Please verify gene names.")
}

# ========== KEGG Enrichment ==========
kegg_real <- enrichKEGG(gene = gene_real$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
kegg_gen  <- enrichKEGG(gene = gene_gen$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)

if (is.null(kegg_real) || is.null(kegg_gen)) {
    stop("❌ enrichKEGG returned NULL. Check if valid genes are recognized by KEGG.")
}

# ========== 保存KEGG结果到文件 ==========
saveRDS(kegg_real, get_output_path("kegg_analysis", "kegg_real_results.rds"))
saveRDS(kegg_gen, get_output_path("kegg_analysis", "kegg_generated_results.rds"))
cat("💾 KEGG分析结果已保存到:\n")
cat("  -", get_output_path("kegg_analysis", "kegg_real_results.rds"), "\n")
cat("  -", get_output_path("kegg_analysis", "kegg_generated_results.rds"), "\n")

# ========== Pathway Comparison ==========
real_terms <- kegg_real@result$Description
gen_terms  <- kegg_gen@result$Description

common_terms <- intersect(real_terms, gen_terms)
only_real    <- setdiff(real_terms, gen_terms)
only_gen     <- setdiff(gen_terms, real_terms)

write.csv(data.frame(Common = common_terms), get_output_path("kegg_analysis", "common_kegg_pathways.csv"), row.names = FALSE)

# ========== Visualization ==========
pdf(get_output_path("kegg_analysis", "KEGG_comparison_and_venn.pdf"), width = 14, height = 8)

# Side-by-side dotplots
par(mfrow = c(1, 2))
print(dotplot(kegg_real, showCategory = 15, title = "Real Data KEGG"))
print(dotplot(kegg_gen, showCategory = 15, title = "Generated Data KEGG"))

# Venn diagram
grid.newpage()
draw.pairwise.venn(
  area1       = length(real_terms),
  area2       = length(gen_terms),
  cross.area  = length(common_terms),
  category    = c("Real", "Generated"),
  fill        = c("skyblue", "pink"),
  alpha       = 0.6,
  cat.pos     = c(-20, 20),
  cat.dist    = 0.05,
  scaled      = FALSE
)

dev.off()
cat("✅ Analysis complete. Results saved to", output_paths$kegg_analysis, "directory\n")
