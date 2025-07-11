# 智能数据集检测函数
detect_current_dataset <- function() {
    # 1. 首先尝试从环境变量获取
    env_dataset <- Sys.getenv("DATASET_NAME", "")
    if (env_dataset != "") {
        cat("📌 从环境变量检测到数据集:", env_dataset, "\n")
        return(env_dataset)
    }
    
    # 2. 检查FD1000目录中的预处理文件
    if (dir.exists("FD1000")) {
        preprocess_files <- list.files("FD1000", pattern = ".*PreProLabel1000\\.csv$", full.names = FALSE)
        if (length(preprocess_files) > 0) {
            # 提取数据集名称
            dataset_names <- gsub("PreProLabel1000\\.csv$", "", preprocess_files)
            if (length(dataset_names) == 1) {
                cat("📊 从预处理文件检测到数据集:", dataset_names[1], "\n")
                return(dataset_names[1])
            } else if (length(dataset_names) > 1) {
                # 如果有多个数据集，优先选择最新的
                file_times <- file.mtime(file.path("FD1000", preprocess_files))
                latest_idx <- which.max(file_times)
                selected_dataset <- dataset_names[latest_idx]
                cat("📊 检测到多个数据集，选择最新的:", selected_dataset, "\n")
                cat("   可用数据集:", paste(dataset_names, collapse = ", "), "\n")
                return(selected_dataset)
            }
        }
    }
    
    # 3. 检查output目录中的生成数据
    if (dir.exists("output")) {
        output_dirs <- list.dirs("output", full.names = FALSE, recursive = FALSE)
        output_dirs <- output_dirs[output_dirs != ""]  # 移除空字符串
        if (length(output_dirs) > 0) {
            # 检查哪些目录包含generated_data
            valid_datasets <- c()
            for (dir_name in output_dirs) {
                gen_data_dir <- file.path("output", dir_name, "generated_data")
                if (dir.exists(gen_data_dir)) {
                    gen_files <- list.files(gen_data_dir, pattern = ".*_generated\\.csv$")
                    if (length(gen_files) > 0) {
                        valid_datasets <- c(valid_datasets, dir_name)
                    }
                }
            }
            
            if (length(valid_datasets) == 1) {
                cat("📁 从输出目录检测到数据集:", valid_datasets[1], "\n")
                return(valid_datasets[1])
            } else if (length(valid_datasets) > 1) {
                # 选择最新修改的
                dir_times <- sapply(valid_datasets, function(x) {
                    max(file.mtime(list.files(file.path("output", x), recursive = TRUE, full.names = TRUE)))
                })
                latest_dataset <- names(which.max(dir_times))
                cat("📁 检测到多个输出数据集，选择最新的:", latest_dataset, "\n")
                cat("   可用数据集:", paste(valid_datasets, collapse = ", "), "\n")
                return(latest_dataset)
            }
        }
    }
    
    # 4. 如果都无法检测到，返回默认值但给出警告
    cat("⚠️  无法自动检测数据集，使用默认值: AD01103\n")
    cat("💡 建议:\n")
    cat("   1. 设置环境变量: export DATASET_NAME=your_dataset\n")
    cat("   2. 确保预处理文件存在于 FD1000/ 目录\n")
    cat("   3. 确保生成数据存在于 output/dataset_name/ 目录\n")
    return("AD01103")
}

# 自动创建输出目录结构
create_output_structure <- function(dataset_name = NULL) {
    # 如果没有指定数据集名称，智能检测当前数据集
    if (is.null(dataset_name)) {
        dataset_name <- detect_current_dataset()
    }
    
    cat("🎯 使用数据集:", dataset_name, "\n")
    
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
    # 如果没有指定数据集名称，智能检测当前数据集
    if (is.null(dataset_name)) {
        dataset_name <- detect_current_dataset()
    }
    paths <- create_output_structure(dataset_name)
    return(file.path(paths[[category]], filename))
}

# 获取数据集相关的输入文件路径
get_dataset_paths <- function(dataset_name = NULL) {
    # 如果没有指定数据集名称，智能检测当前数据集
    if (is.null(dataset_name)) {
        dataset_name <- detect_current_dataset()
    }
    
    return(list(
        dataset_name = dataset_name,
        real_data_file = paste0("FD1000/", dataset_name, "PreProLabel1000.csv"),
        generated_data_file = get_output_path("generated_data", paste0(dataset_name, "_generated.csv"), dataset_name)
    ))
}
