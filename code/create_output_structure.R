# 自动创建输出目录结构
# 用于在分析脚本中调用，确保目录存在

create_output_structure <- function(dataset_name = NULL) {
    # 如果没有指定数据集名称，从环境变量中读取
    if (is.null(dataset_name)) {
        dataset_name <- Sys.getenv("DATASET_NAME", "AD01103")
    }
    # 主输出目录
    base_dir <- "output"
    dataset_dir <- file.path(base_dir, dataset_name)
    
    # 功能分类子目录
    subdirs <- c(
        "generated_data",       # 生成的数据文件
        "visualization",        # 可视化结果（PCA, t-SNE等）
        "clustering",          # 聚类分析结果
        "differential_expression", # 差异表达分析
        "kegg_analysis",       # KEGG富集分析
        "quality_assessment",  # 质量评估结果
        "reports"             # 报告和文档
    )
    
    # 创建目录
    if (!dir.exists(base_dir)) {
        dir.create(base_dir, recursive = TRUE)
    }
    
    if (!dir.exists(dataset_dir)) {
        dir.create(dataset_dir, recursive = TRUE)
    }
    
    for (subdir in subdirs) {
        full_path <- file.path(dataset_dir, subdir)
        if (!dir.exists(full_path)) {
            dir.create(full_path, recursive = TRUE)
            cat("📁 创建目录:", full_path, "\n")
        }
    }
    
    # 返回目录路径列表
    return(list(
        base = base_dir,
        dataset = dataset_dir,
        generated_data = file.path(dataset_dir, "generated_data"),
        visualization = file.path(dataset_dir, "visualization"),
        clustering = file.path(dataset_dir, "clustering"),
        differential_expression = file.path(dataset_dir, "differential_expression"),
        kegg_analysis = file.path(dataset_dir, "kegg_analysis"),
        quality_assessment = file.path(dataset_dir, "quality_assessment"),
        reports = file.path(dataset_dir, "reports")
    ))
}

# 获取标准化文件路径
get_output_path <- function(category, filename, dataset_name = NULL) {
    # 如果没有指定数据集名称，从环境变量中读取
    if (is.null(dataset_name)) {
        dataset_name <- Sys.getenv("DATASET_NAME", "AD01103")
    }
    paths <- create_output_structure(dataset_name)
    return(file.path(paths[[category]], filename))
}

# 获取数据集相关的输入文件路径
get_dataset_paths <- function(dataset_name = NULL) {
    # 如果没有指定数据集名称，从环境变量中读取
    if (is.null(dataset_name)) {
        dataset_name <- Sys.getenv("DATASET_NAME", "AD01103")
    }
    
    return(list(
        dataset_name = dataset_name,
        real_data_file = paste0("FD1000/", dataset_name, "PreProLabel1000.csv"),
        generated_data_file = get_output_path("generated_data", paste0(dataset_name, "_generated.csv"), dataset_name)
    ))
} 