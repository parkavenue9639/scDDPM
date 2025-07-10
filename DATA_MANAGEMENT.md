# 数据管理说明

## 📁 数据文件策略

### 🎯 推荐方案：只上传压缩文件

为了在不同设备间高效协作，我们采用以下策略：

1. **GitHub仓库中只包含压缩文件**（.gz, .tar.gz, .rar）
2. **解压后的文件不提交到Git**（通过.gitignore排除）
3. **使用自动化脚本解压数据**

### 📊 文件大小对比

| 文件类型 | 大小 | 是否上传 |
|---------|------|----------|
| AD01103_expr.txt.gz | ~50MB | ✅ 上传 |
| AD01103_expr.txt | 364MB | ❌ 不上传 |
| AD00203_expr.txt.gz | ~30MB | ✅ 上传 |
| AD00203_expr.txt | 178MB | ❌ 不上传 |

## 🚀 在不同设备上使用

### 首次克隆仓库
```bash
git clone https://github.com/你的用户名/scDDPM.git
cd scDDPM
```

### 获取数据文件
有两种方式：

#### 方式1：使用下载脚本（推荐）
```bash
./download_data.sh
# 然后手动下载数据文件到对应目录
```

#### 方式2：从其他设备复制
```bash
# 从其他设备复制压缩文件
scp -r user@other-device:/path/to/scDDPM/data/*.gz ./data/
```

### 运行完整流程
```bash
./run_complete.sh
```

## 📋 数据文件清单

### 必需文件
- `data/AD01103/AD01103_expr.txt.gz`
- `data/AD01103/AD01103_cell_label.txt.gz`
- `data/AD00202/AD00202_expr.txt.gz`
- `data/AD00202/AD00202_cell_label.txt.gz`
- `data/AD00203/AD00203_expr.txt.gz`
- `data/AD00203/AD00203_cell_label.txt.gz`
- `data/AD00204/AD00204_expr.txt.gz`
- `data/AD00204/AD00204_cell_label.txt.gz`
- `data/AD00401/AD00401_expr.txt.gz`
- `data/AD00401/AD00401_cell_label.txt.gz`
- `data/GSE119911/GSE119911_series_matrix.txt.gz`
- `data/10X PBMC 3k/pbmc3k_filtered_gene_bc_matrices.tar.gz`
- `data/10X PBMC 3k/pbmc3k_analysis.tar.gz`

### 可选文件
- Baron Human/Mouse 数据集（.rar格式）

## 🔧 技术细节

### Git LFS配置
项目已配置Git LFS来处理大文件：
- `.gitattributes` 文件定义了LFS跟踪规则
- 压缩文件使用LFS存储
- 解压文件被.gitignore排除

### 自动化解压
`run_complete.sh` 脚本会自动：
1. 检查压缩文件是否存在
2. 解压必要的文件
3. 运行完整的分析流程

## 💡 最佳实践

1. **定期备份**：将压缩文件备份到云存储
2. **版本控制**：只提交代码和配置文件，不提交大文件
3. **文档记录**：记录数据来源和处理步骤
4. **环境一致性**：在不同设备上使用相同的环境配置

## 🆘 常见问题

### Q: 为什么不上传解压后的文件？
A: 解压后的文件太大（100MB+），超过GitHub限制，且Git版本控制效率低。

### Q: 如何在新设备上快速开始？
A: 克隆仓库 → 下载数据文件 → 运行 `./run_complete.sh`

### Q: 数据文件丢失怎么办？
A: 使用 `./download_data.sh` 脚本重新下载，或从备份恢复。 