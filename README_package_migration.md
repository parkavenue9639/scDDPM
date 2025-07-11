# R包迁移指南

## 📦 导出当前设备的R包

### 1. 运行导出脚本
```bash
Rscript export_packages.R
```

这将生成以下文件：
- `packages_detailed.csv` - 详细包信息
- `packages_simple.csv` - 简化包信息（推荐使用）
- `session_info.rds` - R session信息
- `install_packages.R` - 基础安装脚本
- `install_packages_smart.R` - 智能安装脚本

### 2. 查看包信息
```r
# 查看session信息
session_info <- readRDS("session_info.rds")
print(session_info)

# 查看包列表
packages <- read.csv("packages_simple.csv")
head(packages)
```

## 🚀 在新设备上安装包

### 方法1：使用智能安装脚本（推荐）
```bash
Rscript install_packages_smart.R
```

### 方法2：使用基础安装脚本
```bash
Rscript install_packages.R
```

### 方法3：手动安装
```r
# 读取包列表
packages <- read.csv("packages_simple.csv")

# 安装所有包
install.packages(packages$Package)

# 对于Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 安装Bioconductor包
BiocManager::install(c("SC3", "SingleCellExperiment", "scater", "scran"))
```

## 📋 针对scDDPM项目的包

根据项目需求，以下包是必需的：

### CRAN包
```r
install.packages(c(
    "ggplot2", "dplyr", "tidyr", "readr", "stringr",
    "magrittr", "purrr", "tibble", "forcats", "rlang",
    "glue", "fs", "here", "usethis", "devtools",
    "roxygen2", "testthat", "knitr", "rmarkdown",
    "DT", "plotly", "shiny", "flexdashboard"
))
```

### Bioconductor包
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SC3", "SingleCellExperiment", "scater", "scran",
    "scRNAseq", "Biobase", "BiocGenerics", "S4Vectors",
    "IRanges", "GenomicRanges", "SummarizedExperiment",
    "DelayedArray", "MatrixGenerics", "XVector",
    "GenomeInfoDb", "AnnotationDbi", "org.Hs.eg.db",
    "clusterProfiler", "enrichplot", "DOSE", "pathview"
))
```

## 🔧 故障排除

### 常见问题

1. **包安装失败**
   ```r
   # 更新R和包
   update.packages()
   
   # 检查包依赖
   install.packages("package_name", dependencies = TRUE)
   ```

2. **Bioconductor包问题**
   ```r
   # 更新BiocManager
   BiocManager::install()
   
   # 检查版本兼容性
   BiocManager::valid()
   ```

3. **权限问题**
   ```r
   # 检查库路径
   .libPaths()
   
   # 设置用户库
   .libPaths(c(.libPaths(), "~/R/library"))
   ```

### 验证安装
```r
# 检查关键包是否安装
required_packages <- c("SC3", "SingleCellExperiment", "ggplot2", "dplyr")
for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("✅", pkg, "已安装\n")
    } else {
        cat("❌", pkg, "未安装\n")
    }
}
```

## 📁 文件说明

- `export_packages.R` - 导出包列表脚本
- `install_packages_smart.R` - 智能安装脚本
- `packages_simple.csv` - 包名和版本列表
- `packages_detailed.csv` - 详细包信息
- `session_info.rds` - R环境信息

## 💡 最佳实践

1. **定期更新包列表**：在安装新包后重新导出
2. **使用版本控制**：将包列表文件加入版本控制
3. **测试安装**：在新环境中验证所有包正常工作
4. **记录依赖**：为项目创建专门的包管理文件 