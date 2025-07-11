# Clear environment
rm(list = ls())

# 加载目录结构创建函数
source("code/create_output_structure.R")

# 获取数据集路径信息
dataset_paths <- get_dataset_paths()
dataset_name <- dataset_paths$dataset_name

# 创建输出目录结构
output_paths <- create_output_structure()

# 检查文件是否存在
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

# Load original and generated data
data <- read.csv(real_data_file, header = TRUE, row.names = 1)
gen_data <- read.csv(gen_data_file, header = TRUE, row.names = 1)

# Load necessary libraries
library(Seurat)
library(SeuratObject)
library(magrittr)
library(dplyr)

# Check sparsity
cat("真实数据稀疏度:", sum(data == 0), "/", nrow(data) * ncol(data), "\n")
cat("真实数据维度:", dim(data), "\n")

# Select the first 100 cells and transpose: Seurat requires features as rows and cells as columns
if (ncol(data) > 100) {
    data <- data[, 1:100]
    cat("选择前100个细胞进行分析\n")
}
data <- t(data)

# Create Seurat object for original data
scell <- CreateSeuratObject(counts = data, project = "Alzheimer0205")

# Check Seurat object's internal structure (optional)
cat("Seurat对象维度:", dim(scell@assays$RNA@layers$counts), "\n")

# Identify top 100 highly variable genes using VST method
data1000 <- FindVariableFeatures(scell, selection.method = "vst", nfeatures = 100)
top <- head(VariableFeatures(data1000), 100)

# Save top genes to CSV
write.csv(top, get_output_path("differential_expression", "real_top100_genes.csv"), row.names = FALSE)
cat("✅ 真实数据前100个高变异基因已保存\n")

# ------------------ Repeat for Generated Data ------------------
cat("处理生成数据...\n")
cat("生成数据维度:", dim(gen_data), "\n")

# Select and transpose generated data
if (ncol(gen_data) > 100) {
    gen_data <- gen_data[, 1:100]
    cat("选择前100个细胞进行分析\n")
}
gen_data <- t(gen_data)

# Create Seurat object for generated data
gen_scell <- CreateSeuratObject(counts = gen_data, project = "Alzheimer0205")

# Check structure
cat("生成数据Seurat对象维度:", dim(gen_scell@assays$RNA@layers$counts), "\n")

# Identify top 100 variable genes
gen_data100 <- FindVariableFeatures(gen_scell, selection.method = "vst", nfeatures = 100)
top <- head(VariableFeatures(gen_data100), 100)

# Save top genes
write.csv(top, get_output_path("differential_expression", "generated_top100_genes.csv"), row.names = FALSE)
cat("✅ 生成数据前100个高变异基因已保存\n")

# 比较两个基因列表
real_genes <- read.csv(get_output_path("differential_expression", "real_top100_genes.csv"))$x
gen_genes <- read.csv(get_output_path("differential_expression", "generated_top100_genes.csv"))$x

common_genes <- intersect(real_genes, gen_genes)
cat("共同的高变异基因数量:", length(common_genes), "\n")
cat("真实数据特有基因数量:", length(setdiff(real_genes, gen_genes)), "\n")
cat("生成数据特有基因数量:", length(setdiff(gen_genes, real_genes)), "\n")

# 保存比较结果
comparison_df <- data.frame(
    metric = c("共同基因", "真实数据特有", "生成数据特有"),
    count = c(length(common_genes), length(setdiff(real_genes, gen_genes)), length(setdiff(gen_genes, real_genes)))
)
write.csv(comparison_df, get_output_path("differential_expression", "gene_comparison.csv"), row.names = FALSE)
write.csv(data.frame(gene = common_genes), get_output_path("differential_expression", "common_genes.csv"), row.names = FALSE)

cat("✅ 差异表达分析完成！结果保存在", output_paths$differential_expression, "目录\n")
