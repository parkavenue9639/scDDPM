#!/bin/bash

echo "================= scDDPM 数据下载脚本 ================="
echo "此脚本用于从GitHub仓库获取scDDPM项目所需的原始数据"
echo "=========================================================="

# 检查是否在Git仓库中
if [ ! -d ".git" ]; then
    echo "❌ 错误：请在Git仓库根目录下运行此脚本"
    exit 1
fi

# 检查远程仓库配置
if ! git remote get-url origin > /dev/null 2>&1; then
    echo "❌ 错误：未找到远程仓库配置"
    exit 1
fi

echo "🔍 检查GitHub仓库中的数据文件..."

# 获取远程仓库URL
REPO_URL=$(git remote get-url origin)
echo "📁 仓库地址: $REPO_URL"

# 从GitHub获取数据文件
echo "📥 从GitHub仓库下载数据文件..."

# 检查并下载AD00202数据
echo "📂 检查 AD00202 数据..."
if [ ! -f "data/AD00202/AD00202_expr.txt.gz" ]; then
    echo "  - 缺失: AD00202_expr.txt.gz"
    echo "  - 请从GitHub仓库下载或联系项目维护者"
else
    echo "  ✅ 已存在: AD00202_expr.txt.gz"
fi

if [ ! -f "data/AD00202/AD00202_cell_label.txt.gz" ]; then
    echo "  - 缺失: AD00202_cell_label.txt.gz"
    echo "  - 请从GitHub仓库下载或联系项目维护者"
else
    echo "  ✅ 已存在: AD00202_cell_label.txt.gz"
fi

# 检查并下载AD00203数据
echo "📂 检查 AD00203 数据..."
if [ ! -f "data/AD00203/AD00203_expr.txt.gz" ]; then
    echo "  - 缺失: AD00203_expr.txt.gz"
else
    echo "  ✅ 已存在: AD00203_expr.txt.gz"
fi

if [ ! -f "data/AD00203/AD00203_cell_label.txt.gz" ]; then
    echo "  - 缺失: AD00203_cell_label.txt.gz"
else
    echo "  ✅ 已存在: AD00203_cell_label.txt.gz"
fi

# 检查并下载AD00204数据
echo "📂 检查 AD00204 数据..."
if [ ! -f "data/AD00204/AD00204_expr.txt.gz" ]; then
    echo "  - 缺失: AD00204_expr.txt.gz"
else
    echo "  ✅ 已存在: AD00204_expr.txt.gz"
fi

if [ ! -f "data/AD00204/AD00204_cell_label.txt.gz" ]; then
    echo "  - 缺失: AD00204_cell_label.txt.gz"
else
    echo "  ✅ 已存在: AD00204_cell_label.txt.gz"
fi

# 检查并下载AD00401数据
echo "📂 检查 AD00401 数据..."
if [ ! -f "data/AD00401/AD00401_expr.txt.gz" ]; then
    echo "  - 缺失: AD00401_expr.txt.gz"
else
    echo "  ✅ 已存在: AD00401_expr.txt.gz"
fi

if [ ! -f "data/AD00401/AD00401_cell_label.txt.gz" ]; then
    echo "  - 缺失: AD00401_cell_label.txt.gz"
else
    echo "  ✅ 已存在: AD00401_cell_label.txt.gz"
fi

# 检查并下载AD01103数据
echo "📂 检查 AD01103 数据..."
if [ ! -f "data/AD01103/AD01103_expr.txt.gz" ]; then
    echo "  - 缺失: AD01103_expr.txt.gz"
else
    echo "  ✅ 已存在: AD01103_expr.txt.gz"
fi

if [ ! -f "data/AD01103/AD01103_cell_label.txt.gz" ]; then
    echo "  - 缺失: AD01103_cell_label.txt.gz"
else
    echo "  ✅ 已存在: AD01103_cell_label.txt.gz"
fi

# 检查并下载GSE119911数据
echo "📂 检查 GSE119911 数据..."
if [ ! -f "data/GSE119911/GSE119911_series_matrix.txt.gz" ]; then
    echo "  - 缺失: GSE119911_series_matrix.txt.gz"
else
    echo "  ✅ 已存在: GSE119911_series_matrix.txt.gz"
fi

# 检查并下载10X PBMC 3k数据
echo "📂 检查 10X PBMC 3k 数据..."
if [ ! -f "data/10X PBMC 3k/pbmc3k_filtered_gene_bc_matrices.tar.gz" ]; then
    echo "  - 缺失: pbmc3k_filtered_gene_bc_matrices.tar.gz"
else
    echo "  ✅ 已存在: pbmc3k_filtered_gene_bc_matrices.tar.gz"
fi

if [ ! -f "data/10X PBMC 3k/pbmc3k_analysis.tar.gz" ]; then
    echo "  - 缺失: pbmc3k_analysis.tar.gz"
else
    echo "  ✅ 已存在: pbmc3k_analysis.tar.gz"
fi

echo ""
echo "📊 数据状态统计："

# 统计缺失和存在的文件
missing_count=0
existing_count=0

for file in \
    "data/AD00202/AD00202_expr.txt.gz" \
    "data/AD00202/AD00202_cell_label.txt.gz" \
    "data/AD00203/AD00203_expr.txt.gz" \
    "data/AD00203/AD00203_cell_label.txt.gz" \
    "data/AD00204/AD00204_expr.txt.gz" \
    "data/AD00204/AD00204_cell_label.txt.gz" \
    "data/AD00401/AD00401_expr.txt.gz" \
    "data/AD00401/AD00401_cell_label.txt.gz" \
    "data/AD01103/AD01103_expr.txt.gz" \
    "data/AD01103/AD01103_cell_label.txt.gz" \
    "data/GSE119911/GSE119911_series_matrix.txt.gz" \
    "data/10X PBMC 3k/pbmc3k_filtered_gene_bc_matrices.tar.gz" \
    "data/10X PBMC 3k/pbmc3k_analysis.tar.gz"; do
    if [ -f "$file" ]; then
        ((existing_count++))
    else
        ((missing_count++))
    fi
done

echo "  - 已存在文件: $existing_count 个"
echo "  - 缺失文件: $missing_count 个"

if [ $missing_count -eq 0 ]; then
    echo ""
    echo "🎉 所有数据文件都已存在！"
    echo "✅ 可以直接运行 ./run_complete.sh 开始分析流程"
else
    echo ""
    echo "⚠️  部分数据文件缺失，请："
    echo "1. 从GitHub仓库下载缺失的文件"
    echo "2. 或联系项目维护者获取完整数据"
    echo "3. 或使用 git lfs pull 拉取LFS文件"
fi

echo ""
echo "🔧 获取缺失文件的方法："
echo "1. 使用 git lfs pull 拉取Git LFS文件"
echo "2. 从GitHub网页直接下载"
echo "3. 使用 git clone --recurse-submodules 重新克隆"

# 询问是否自动拉取LFS文件
if [ $missing_count -gt 0 ]; then
    echo ""
    read -p "🤔 是否自动拉取Git LFS文件？(y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "📥 正在拉取Git LFS文件..."
        git lfs pull
        echo "✅ Git LFS文件拉取完成"
        
        # 重新检查文件状态
        echo ""
        echo "🔄 重新检查文件状态..."
        missing_count=0
        existing_count=0
        
        for file in \
            "data/AD00202/AD00202_expr.txt.gz" \
            "data/AD00202/AD00202_cell_label.txt.gz" \
            "data/AD00203/AD00203_expr.txt.gz" \
            "data/AD00203/AD00203_cell_label.txt.gz" \
            "data/AD00204/AD00204_expr.txt.gz" \
            "data/AD00204/AD00204_cell_label.txt.gz" \
            "data/AD00401/AD00401_expr.txt.gz" \
            "data/AD00401/AD00401_cell_label.txt.gz" \
            "data/AD01103/AD01103_expr.txt.gz" \
            "data/AD01103/AD01103_cell_label.txt.gz" \
            "data/GSE119911/GSE119911_series_matrix.txt.gz" \
            "data/10X PBMC 3k/pbmc3k_filtered_gene_bc_matrices.tar.gz" \
            "data/10X PBMC 3k/pbmc3k_analysis.tar.gz"; do
            if [ -f "$file" ]; then
                ((existing_count++))
            else
                ((missing_count++))
            fi
        done
        
        echo "  - 已存在文件: $existing_count 个"
        echo "  - 缺失文件: $missing_count 个"
        
        if [ $missing_count -eq 0 ]; then
            echo ""
            echo "🎉 所有数据文件都已获取！"
        fi
    fi
fi
echo ""
echo "📋 下一步："
echo "运行 ./run_complete.sh 开始自动化分析流程" 