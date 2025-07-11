# KEGG结果加载和管理脚本
# 用于加载已保存的KEGG分析结果到R环境中

# 加载目录结构创建函数
source("code/create_output_structure.R")

# ========== 加载KEGG结果函数 ==========
load_kegg_results <- function(verbose = TRUE, dataset_name = NULL) {
  # 创建输出目录结构（如果不存在）
  output_paths <- create_output_structure(dataset_name)
  
  kegg_real_file <- get_output_path("kegg_analysis", "kegg_real_results.rds")
  kegg_gen_file <- get_output_path("kegg_analysis", "kegg_generated_results.rds")
  
  if (file.exists(kegg_real_file) && file.exists(kegg_gen_file)) {
    if (verbose) cat("📂 加载KEGG分析结果...\n")
    
    kegg_real <<- readRDS(kegg_real_file)
    kegg_gen <<- readRDS(kegg_gen_file)
    
    if (verbose) {
      cat("✅ KEGG结果加载成功!\n")
      cat("📊 结果概要:\n")
      cat("  - 真实数据富集通路数:", nrow(kegg_real@result), "\n")
      cat("  - 生成数据富集通路数:", nrow(kegg_gen@result), "\n")
      cat("  - 文件更新时间:\n")
      cat("    * 真实数据:", format(file.mtime(kegg_real_file), "%Y-%m-%d %H:%M:%S"), "\n")
      cat("    * 生成数据:", format(file.mtime(kegg_gen_file), "%Y-%m-%d %H:%M:%S"), "\n")
    }
    
    return(TRUE)
  } else {
    if (verbose) {
      cat("❌ KEGG结果文件不存在:\n")
      if (!file.exists(kegg_real_file)) cat("  缺失:", kegg_real_file, "\n")
      if (!file.exists(kegg_gen_file)) cat("  缺失:", kegg_gen_file, "\n")
      cat("💡 请先运行 'Rscript code/kegg.R' 生成结果\n")
    }
    return(FALSE)
  }
}

# ========== 检查KEGG结果状态函数 ==========
check_kegg_status <- function(dataset_name = NULL) {
  cat("🔍 KEGG分析结果状态检查:\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  
  # 创建输出目录结构（如果不存在）
  output_paths <- create_output_structure(dataset_name)
  
  # 检查文件
  kegg_real_file <- get_output_path("kegg_analysis", "kegg_real_results.rds")
  kegg_gen_file <- get_output_path("kegg_analysis", "kegg_generated_results.rds")
  
  cat("📁 文件状态:\n")
  if (file.exists(kegg_real_file)) {
    cat("  ✅ 真实数据结果:", kegg_real_file, "\n")
    cat("     大小:", format(file.size(kegg_real_file)/1024, digits=2), "KB\n")
    cat("     修改时间:", format(file.mtime(kegg_real_file), "%Y-%m-%d %H:%M:%S"), "\n")
  } else {
    cat("  ❌ 真实数据结果: 文件不存在\n")
  }
  
  if (file.exists(kegg_gen_file)) {
    cat("  ✅ 生成数据结果:", kegg_gen_file, "\n")
    cat("     大小:", format(file.size(kegg_gen_file)/1024, digits=2), "KB\n")
    cat("     修改时间:", format(file.mtime(kegg_gen_file), "%Y-%m-%d %H:%M:%S"), "\n")
  } else {
    cat("  ❌ 生成数据结果: 文件不存在\n")
  }
  
  # 检查内存
  cat("\n💾 内存状态:\n")
  if (exists("kegg_real")) {
    cat("  ✅ kegg_real: 已加载到内存\n")
  } else {
    cat("  ❌ kegg_real: 未在内存中\n")
  }
  
  if (exists("kegg_gen")) {
    cat("  ✅ kegg_gen: 已加载到内存\n")
  } else {
    cat("  ❌ kegg_gen: 未在内存中\n")
  }
  
  # 给出建议
  cat("\n💡 建议操作:\n")
  files_exist <- file.exists(kegg_real_file) && file.exists(kegg_gen_file)
  memory_exists <- exists("kegg_real") && exists("kegg_gen")
  
  if (!files_exist && !memory_exists) {
    cat("  🚀 运行: Rscript code/kegg.R (首次分析)\n")
  } else if (files_exist && !memory_exists) {
    cat("  📂 运行: load_kegg_results() (从文件加载)\n")
  } else if (!files_exist && memory_exists) {
    cat("  💾 结果在内存中，建议重新运行kegg.R保存到文件\n")
  } else {
    cat("  ✅ 结果已准备就绪，可以进行分析\n")
  }
}

# ========== 清理KEGG结果函数 ==========
clear_kegg_results <- function(remove_files = FALSE, remove_memory = TRUE, dataset_name = NULL) {
  if (remove_memory) {
    if (exists("kegg_real")) {
      rm(kegg_real, envir = .GlobalEnv)
      cat("🗑️  已从内存清除: kegg_real\n")
    }
    if (exists("kegg_gen")) {
      rm(kegg_gen, envir = .GlobalEnv)
      cat("🗑️  已从内存清除: kegg_gen\n")
    }
  }
  
  if (remove_files) {
    # 创建输出目录结构（如果不存在）
    output_paths <- create_output_structure(dataset_name)
    
    kegg_real_file <- get_output_path("kegg_analysis", "kegg_real_results.rds")
    kegg_gen_file <- get_output_path("kegg_analysis", "kegg_generated_results.rds")
    
    if (file.exists(kegg_real_file)) {
      file.remove(kegg_real_file)
      cat("🗑️  已删除文件:", kegg_real_file, "\n")
    }
    if (file.exists(kegg_gen_file)) {
      file.remove(kegg_gen_file)
      cat("🗑️  已删除文件:", kegg_gen_file, "\n")
    }
  }
}

# ========== 主函数调用 ==========
if (!interactive()) {
  # 如果作为脚本运行，显示状态
  check_kegg_status()
} else {
  # 如果在交互式环境中，提供帮助信息
  cat("📋 KEGG结果管理命令:\n")
  cat("  - load_kegg_results()     # 加载已保存的结果\n")
  cat("  - check_kegg_status()     # 检查结果状态\n")
  cat("  - clear_kegg_results()    # 清理结果\n")
} 