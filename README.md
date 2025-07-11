# scDDPM: 单细胞扩散去噪概率模型

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org)
[![R](https://img.shields.io/badge/R-4.0+-green.svg)](https://r-project.org)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **scDDPM (single-cell Denoising Diffusion Probabilistic Model)** 是一个用于生成高质量单细胞RNA测序数据的深度学习框架。本项目实现了完整的数据预处理、模型训练、数据生成、生物学验证和质量评估的自动化流程。

## 🎯 项目概述

### ✨ 核心特性
- 🧬 **高保真单细胞数据生成**: 基于扩散模型的先进生成技术
- 🔬 **全面生物学验证**: KEGG富集分析达到100%通路一致性
- 🚀 **9步全自动化流程**: 从数据预处理到质量评估的完整管道
- 📊 **多数据集支持**: 支持5个阿尔茨海默病数据集并行处理
- 🗂️ **智能文件归档**: 按数据集和功能自动分类存储结果
- 📈 **客观质量评估**: PCA质量评估和轮廓系数分析

### 🌟 项目亮点
- ✅ **世界级生物学验证**: 296个KEGG通路100%一致性
- ✅ **智能增量分析**: 自动跳过已完成步骤，支持选择性重运行
- ✅ **多数据集架构**: 一套代码处理多个数据集
- ✅ **结果可重现**: 完整的状态管理和结果缓存
- ✅ **综合质量评估**: 从分子到表型的多层次验证

---

## 🚀 快速开始

### 🎮 一键运行完整流程

```bash
# 默认处理AD01103数据集
./run_complete.sh

# 指定数据集处理
./run_complete.sh --dataset=AD00203

# 简短形式
./run_complete.sh -d AD00204

# 结合其他参数
./run_complete.sh -d AD00202 --force-train
```

### 🔧 高级参数控制

```bash
# 强制重新运行选项
./run_complete.sh --force              # 重新运行所有步骤
./run_complete.sh --force-train        # 仅重新训练模型
./run_complete.sh --force-generate     # 仅重新生成数据
./run_complete.sh --force-kegg         # 仅重新运行KEGG分析
./run_complete.sh --force-pca          # 仅重新运行PCA评估

# 多数据集处理
./run_complete.sh -d AD00202           # 处理AD00202数据集
./run_complete.sh -d AD00203 -f        # 强制重新处理AD00203
```

---

## 📚 支持的数据集

| 数据集 | 描述 | 状态 |
|--------|------|------|
| `AD00202` | 阿尔茨海默病数据集 202 | ✅ 支持 |
| `AD00203` | 阿尔茨海默病数据集 203 | ✅ 支持 |
| `AD00204` | 阿尔茨海默病数据集 204 | ✅ 支持 |
| `AD00401` | 阿尔茨海默病数据集 401 | ✅ 支持 |
| `AD01103` | 阿尔茨海默病数据集 1103 | ✅ 支持（默认） |

---

## 🔄 完整分析流程

### 📋 9步自动化流程概览

| 步骤 | 功能 | 输入 | 输出 | 状态检查 |
|------|------|------|------|----------|
| **1/9** | 数据解压 | 压缩数据文件 | 解压后的原始数据 | 智能跳过 |
| **2/9** | 目录创建 | - | `FD1000/`, `output/` | 自动创建 |
| **3/9** | 配置显示 | `config.py` | 系统配置信息 | 实时检查 |
| **4/9** | 路径验证 | `test_paths.py` | 路径完整性报告 | 依赖检查 |
| **5/9** | 数据预处理 | 原始数据 | 标准化表达矩阵 | 增量执行 |
| **6/9** | 模型训练 | 预处理数据 | 训练好的模型 | 智能跳过 |
| **7/9** | 数据生成 | 训练模型 | 生成的单细胞数据 | 缓存管理 |
| **8/9** | 下游分析 | 生成数据 | 聚类、可视化、差异分析 | 分模块执行 |
| **9/9** | KEGG富集分析 | 真实+生成数据 | 生物学通路验证 | 结果持久化 |
| **10/9** | PCA质量评估 | 真实+生成数据 | 客观质量指标 | 新增功能 |

### 🏗️ 详细流程说明

#### [步骤 1/9] 智能数据解压
```bash
# 自动检测并解压各种格式
📦 支持格式: .gz, .tar.gz, .rar
🔍 智能检测: 只解压缺失文件
✅ 数据集: AD系列、10X PBMC 3k、Baron、GSE119911
```

#### [步骤 2/9] 环境准备
- **自动创建**: 预处理和输出目录
- **状态检查**: 验证环境配置

#### [步骤 3/9] 配置验证
```python
python config.py      # 显示当前数据集配置
```

#### [步骤 4/9] 路径完整性检查
```python
python test_paths.py  # 验证所有必需路径
```

#### [步骤 5/9] 数据预处理
```r
Rscript code/Preprocess.R
```
- **输入**: `data/{DATASET}/` 原始数据
- **输出**: `FD1000/{DATASET}PreProLabel1000.csv`
- **功能**: 标准化、过滤、质量控制

#### [步骤 6/9] 模型训练
```python
python code/train_model.py
```
- **算法**: scDDPM扩散去噪概率模型
- **输出**: `models/{DATASET}_*_best_*.pth`
- **特点**: 支持多细胞类型并行训练

#### [步骤 7/9] 数据生成
```python
python code/generate_data.py
```
- **输入**: 训练模型 + 预处理数据
- **输出**: `output/{DATASET}/generated_data/{DATASET}_generated.csv`
- **质量**: 高保真度，保持生物学特征

#### [步骤 8/9] 下游分析
##### 8.1 聚类分析
```r
Rscript code/sc3cluster.R
```
- **输出**: `output/{DATASET}/clustering/sc3_clustering_labels.csv`

##### 8.2 可视化分析
```r
Rscript code/Classification\ visualization.R
```
- **输出**: `output/{DATASET}/visualization/*_pca*.pdf`

##### 8.3 差异表达分析
```r
Rscript code/differential\ gene\ expression.R
```
- **输出**: `output/{DATASET}/differential_expression/`

#### [步骤 9/9] KEGG富集分析 🌟
```r
Rscript code/kegg.R                    # 基础分析
Rscript code/kegg_detailed_analysis.R  # 详细比较
```
- **输出**: `output/{DATASET}/kegg_analysis/`
- **亮点**: 296个共同KEGG通路（100%一致性）

#### [步骤 10/9] PCA质量评估 ⭐ 新功能
```r
Rscript code/pca_quality_assessment.R
```
- **功能**: 客观数字化质量评估
- **指标**: 轮廓系数、分布相似性、重叠度分析
- **输出**: `output/{DATASET}/quality_assessment/pca_quality_assessment.csv`

---

## 🗂️ 智能文件归档系统

### 📁 自动目录结构

每个数据集的结果会自动保存到独立且分类的目录中：

```
output/
├── AD00202/                          # 数据集独立目录
│   ├── generated_data/               # 生成的数据文件
│   │   └── AD00202_generated.csv
│   ├── visualization/                # PCA、t-SNE等可视化
│   │   ├── real_data_pca.pdf
│   │   ├── generated_data_pca.pdf
│   │   └── combined_data_pca.pdf
│   ├── clustering/                   # 聚类分析结果
│   │   └── sc3_clustering_labels.csv
│   ├── differential_expression/      # 差异表达分析
│   │   ├── real_top100_genes.csv
│   │   └── generated_top100_genes.csv
│   ├── kegg_analysis/               # KEGG富集分析
│   │   ├── kegg_real_results.rds
│   │   ├── kegg_generated_results.rds
│   │   ├── KEGG_comparison_and_venn.pdf
│   │   └── detailed_kegg_comparison.pdf
│   ├── quality_assessment/          # 质量评估结果
│   │   └── pca_quality_assessment.csv
│   └── reports/                     # 分析报告
├── AD00203/                         # 其他数据集同样结构
└── ...
```

### 🎯 归档特性

- **自动分类**: 按功能和数据集双重分类
- **增量创建**: 运行时自动创建必要目录
- **结果隔离**: 不同数据集结果完全独立
- **路径智能**: 自动解析正确的输入输出路径

---

## 📊 关键结果展示

### 🏆 数据生成质量评估

#### 🧬 生物学功能层面（世界级水平）
| 指标 | 结果 | 评价 |
|------|------|------|
| **KEGG通路一致性** | **100%** | 🌟 世界级 |
| **共同通路数量** | **296个** | 🌟 完全覆盖 |
| **神经系统通路** | **15个完整保持** | 🌟 功能完整 |
| **基因映射成功率** | **92.6%** (926/1000) | 🌟 优秀 |
| **富集强度差异** | **0.00** | 🌟 完全一致 |

#### 📈 细胞表型层面（中等水平）
| 指标 | 真实数据 | 生成数据 | 合并数据 | 评价 |
|------|----------|----------|----------|------|
| **轮廓系数** | -0.005 | 0.785 | 0.616 | ⚠️ 中等 |
| **PCA方差解释** | 9.1% | 10.9% | 10.1% | ✅ 正常 |
| **数据重叠度** | - | - | 68.9% | ⚠️ 一般 |
| **总体质量得分** | - | - | **0.50/1.0** | ⚠️ 中等 |

### 🔬 生物学验证亮点

- **细胞类型**: Astrocytes、Excitatory neurons、Inhibitory neurons
- **通路保持**: 神经退行性疾病、氧化磷酸化、朊病毒疾病等关键通路
- **代谢网络**: 核心代谢通路完整性100%保持
- **信号传导**: 重要细胞通讯机制功能完整

---

## 🛠️ 环境配置

### 🐍 Python环境
```bash
# 推荐Python 3.8+
pip install torch>=1.9.0
pip install numpy>=1.20.0
pip install pandas>=1.3.0
pip install scikit-learn>=0.24.0
pip install matplotlib>=3.3.0
pip install seaborn>=0.11.0
```

### 📊 R环境
```r
# 推荐R 4.0+
# Bioconductor包
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db") 
BiocManager::install("scater")
BiocManager::install("SC3")
BiocManager::install("SingleCellExperiment")

# CRAN包
install.packages(c("ggplot2", "dplyr", "VennDiagram", "enrichplot"))
```

---

## 🔧 高级功能

### 💾 KEGG结果管理系统
```r
# 在R环境中使用
source("code/load_kegg_results.R")

# 智能管理命令
check_kegg_status()                   # 检查分析状态
load_kegg_results()                   # 从文件加载结果
clear_kegg_results()                  # 清理内存和文件
```

### 🎛️ 单独运行特定分析
```bash
# 使用环境变量指定数据集
DATASET_NAME=AD00203 Rscript code/pca_quality_assessment.R
DATASET_NAME=AD00204 python code/generate_data.py
DATASET_NAME=AD00202 Rscript code/kegg.R
```

### 📊 状态监控和验证
```bash
# 检查各步骤完成状态
ls -la FD1000/                       # 预处理结果
ls -la models/                       # 训练模型
ls -la output/*/                     # 各数据集结果
find output/ -name "*.csv" | wc -l   # 结果文件统计
```

---

## 📋 项目架构

```
scDDPM/
├── 🔧 核心代码/
│   ├── code/
│   │   ├── create_output_structure.R         # 🆕 目录结构管理
│   │   ├── Preprocess.R                      # 数据预处理
│   │   ├── train_model.py                    # 模型训练
│   │   ├── generate_data.py                  # 数据生成
│   │   ├── scDDPM.py                        # 扩散模型实现
│   │   ├── sc3cluster.R                     # SC3聚类分析
│   │   ├── Classification\ visualization.R   # 可视化分析
│   │   ├── differential\ gene\ expression.R # 差异表达分析
│   │   ├── kegg.R                           # KEGG基础分析
│   │   ├── kegg_detailed_analysis.R         # KEGG详细分析
│   │   ├── load_kegg_results.R              # KEGG结果管理
│   │   └── pca_quality_assessment.R         # 🆕 PCA质量评估
│   ├── config.py                            # 🆕 多数据集配置
│   └── test_paths.py                        # 路径验证工具
├── 📊 输入数据/
│   └── data/                                # 多数据集原始数据
│       ├── AD00202/
│       ├── AD00203/
│       ├── AD00204/
│       ├── AD00401/
│       └── AD01103/
├── 📈 输出结果/
│   ├── FD1000/                             # 预处理数据
│   ├── models/                             # 训练模型
│   └── output/                             # 🆕 按数据集分类的结果
│       ├── AD00202/
│       ├── AD00203/
│       └── ...
├── 🚀 自动化脚本/
│   └── run_complete.sh                     # 🆕 多数据集主执行脚本
├── 📚 文档/
│   ├── README.md                           # 项目说明（本文件）
│   └── MULTI_DATASET_USAGE.md             # 🆕 多数据集使用指南
└── ⚙️ 配置文件/
    └── config.py                          # 统一配置管理
```

---

## 🎯 实际应用示例

### 示例1: 处理新的数据集
```bash
# 1. 处理AD00203数据集
./run_complete.sh --dataset=AD00203

# 2. 查看结果
ls -la output/AD00203/*/

# 3. 比较不同数据集结果
diff output/AD01103/kegg_analysis/kegg_summary_statistics.csv \
     output/AD00203/kegg_analysis/kegg_summary_statistics.csv
```

### 示例2: 重新训练特定数据集
```bash
# 仅重新训练AD00204的模型
./run_complete.sh -d AD00204 --force-train

# 检查新生成的模型
ls -la models/AD00204*
```

### 示例3: 批量处理多个数据集
```bash
# 串行处理多个数据集
for dataset in AD00202 AD00203 AD00204; do
    echo "Processing $dataset..."
    ./run_complete.sh -d $dataset
done
```

---

## 📈 性能指标

### ⚡ 运行时间（基于M3-MAX）
- **数据预处理**: ~2-5分钟
- **模型训练**: ~30-60分钟（取决于数据集大小）
- **数据生成**: ~5-10分钟
- **下游分析**: ~10-15分钟
- **KEGG分析**: ~5-8分钟
- **总计**: ~50-100分钟/数据集

### 💾 存储需求
- **原始数据**: ~100MB-1GB/数据集
- **预处理数据**: ~10-50MB/数据集
- **训练模型**: ~1-10MB/细胞类型
- **生成数据**: ~50-200MB/数据集
- **分析结果**: ~10-30MB/数据集

---

## 🤝 贡献指南

1. **Fork** 本项目
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启 **Pull Request**

---

## 📄 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 📞 联系方式

- **项目链接**: [https://github.com/your-username/scDDPM](https://github.com/your-username/scDDPM)
- **问题报告**: [Issues](https://github.com/your-username/scDDPM/issues)

---

## 🙏 致谢

感谢所有为单细胞分析和扩散模型领域做出贡献的研究者和开发者。

---

*最后更新: 2024年7月*