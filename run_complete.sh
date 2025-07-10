#!/bin/bash

set -e

# 检查是否有强制重新运行的参数
FORCE_RERUN=false
if [ "$1" = "--force" ] || [ "$1" = "-f" ]; then
    FORCE_RERUN=true
    echo "🔄 强制重新运行模式已启用"
fi

echo "================= scDDPM 完整自动化脚本 ================="
echo "当前时间: $(date)"
echo "工作目录: $(pwd)"
echo "使用方法: $0 [--force|-f]  # 添加 --force 或 -f 参数强制重新运行"
echo "----------------------------------------------------------"

# ================== 智能数据解压模块 ==================
echo "🔍 [1/7] 检查并解压数据..."

# 函数：检查文件是否存在
check_and_extract() {
    local compressed_file="$1"
    local extracted_file="$2"
    local extract_cmd="$3"
    
    if [ -f "$compressed_file" ]; then
        if [ ! -f "$extracted_file" ]; then
            echo "📦 解压文件: $compressed_file"
            eval "$extract_cmd"
            if [ $? -eq 0 ]; then
                echo "✅ 解压成功: $extracted_file"
            else
                echo "❌ 解压失败: $compressed_file"
                return 1
            fi
        else
            echo "✅ 文件已存在: $extracted_file"
        fi
    else
        echo "⚠️  文件不存在: $compressed_file"
    fi
}

# 检查并解压AD系列数据集
echo "📂 处理AD系列数据集..."
for dataset in AD00202 AD00203 AD00204 AD00401 AD01103; do
    if [ -d "data/$dataset" ]; then
        cd "data/$dataset"
        
        # 解压表达矩阵
        check_and_extract \
            "${dataset}_expr.txt.gz" \
            "${dataset}_expr.txt" \
            "gunzip -f ${dataset}_expr.txt.gz"
        
        # 解压细胞标签
        check_and_extract \
            "${dataset}_cell_label.txt.gz" \
            "${dataset}_cell_label.txt" \
            "gunzip -f ${dataset}_cell_label.txt.gz"
        
        cd ../..
    fi
done

# 检查并解压10X PBMC 3k数据集
echo "📂 处理10X PBMC 3k数据集..."
if [ -d "data/10X PBMC 3k" ]; then
    cd "data/10X PBMC 3k"
    
    # 检查tar.gz文件
    for tar_file in *.tar.gz; do
        if [ -f "$tar_file" ]; then
            dir_name="${tar_file%.tar.gz}"
            if [ ! -d "$dir_name" ]; then
                echo "📦 解压文件: $tar_file"
                tar -xzf "$tar_file"
                if [ $? -eq 0 ]; then
                    echo "✅ 解压成功: $dir_name"
                else
                    echo "❌ 解压失败: $tar_file"
                fi
            else
                echo "✅ 目录已存在: $dir_name"
            fi
        fi
    done
    
    cd ../..
fi

# 检查并解压Baron数据集
echo "📂 处理Baron数据集..."
for dataset in "Baron Human" "Baron Mouse"; do
    if [ -d "data/$dataset" ]; then
        cd "data/$dataset"
        
        # 检查rar文件
        for rar_file in *.rar; do
            if [ -f "$rar_file" ]; then
                dir_name="${rar_file%.rar}"
                if [ ! -d "$dir_name" ] && [ ! -f "$dir_name" ]; then
                    echo "📦 解压文件: $rar_file"
                    if command -v unrar &> /dev/null; then
                        unrar x "$rar_file"
                    elif command -v rar &> /dev/null; then
                        rar x "$rar_file"
                    else
                        echo "⚠️  未找到unrar或rar命令，请手动解压: $rar_file"
                        echo "   安装命令: brew install unrar (macOS) 或 apt-get install unrar (Ubuntu)"
                    fi
                else
                    echo "✅ 文件已存在: $dir_name"
                fi
            fi
        done
        
        cd ../..
    fi
done

# 检查GSE119911数据集
echo "📂 处理GSE119911数据集..."
if [ -d "data/GSE119911" ]; then
    cd "data/GSE119911"
    
    for gz_file in *.gz; do
        if [ -f "$gz_file" ]; then
            txt_file="${gz_file%.gz}"
            check_and_extract "$gz_file" "$txt_file" "gunzip -f $gz_file"
        fi
    done
    
    cd ../..
fi

echo "✅ 数据解压检查完成！"

# ================== 创建输出目录 ==================
echo "📁 [2/7] 创建输出目录..."
mkdir -p FD1000
mkdir -p output

# ================== 显示配置信息 ==================
echo "🔧 [3/7] 显示配置信息..."
python config.py

# ================== 路径检测 ==================
echo "🔍 [4/7] 检查路径和文件完整性..."
python test_paths.py

# ================== 数据预处理 ==================
echo "🧹 [5/7] 检查数据预处理..."

# 检查预处理后的文件是否已存在
preprocessed_file="FD1000/AD01103PreProLabel1000.csv"
if [ -f "$preprocessed_file" ] && [ "$FORCE_RERUN" = false ]; then
    echo "✅ 预处理数据已存在: $preprocessed_file"
    echo "   - 文件大小: $(du -h "$preprocessed_file" | cut -f1)"
    echo "   - 跳过预处理步骤"
else
    if [ "$FORCE_RERUN" = true ]; then
        echo "🔄 强制重新运行：开始预处理..."
    else
        echo "📊 预处理数据不存在，开始运行预处理..."
    fi
    Rscript code/Preprocess.R
    
    # 检查预处理结果
    if [ ! -f "$preprocessed_file" ]; then
        echo "❌ 预处理输出文件不存在，流程中止！"
        exit 1
    fi
    echo "✅ 预处理完成，文件已生成: $preprocessed_file"
fi

# ================== 数据生成 ==================
echo "🧠 [6/7] 检查数据生成..."

# 检查生成后的文件是否已存在
generated_file="output/AD01103_generated.csv"
if [ -f "$generated_file" ] && [ "$FORCE_RERUN" = false ]; then
    echo "✅ 生成数据已存在: $generated_file"
    echo "   - 文件大小: $(du -h "$generated_file" | cut -f1)"
    echo "   - 跳过生成步骤"
else
    if [ "$FORCE_RERUN" = true ]; then
        echo "🔄 强制重新运行：开始数据生成..."
    else
        echo "🧠 生成数据不存在，开始运行scDDPM生成模型..."
    fi
    python code/scDDPM.py
    
    # 检查生成结果
    if [ ! -f "$generated_file" ]; then
        echo "❌ 生成数据文件不存在，流程中止！"
        exit 1
    fi
    echo "✅ 数据生成完成，文件已生成: $generated_file"
fi

# ================== 下游分析 ==================
echo "📊 [7/7] 检查下游分析..."

# 检查下游分析结果文件
analysis_files=(
    "output/clustering_results.csv"
    "output/visualization_plots.pdf"
    "output/differential_expression.csv"
    "output/kegg_pathways.csv"
)

# 检查是否所有分析文件都存在
all_analysis_exist=true
for file in "${analysis_files[@]}"; do
    if [ ! -f "$file" ]; then
        all_analysis_exist=false
        break
    fi
done

if [ "$all_analysis_exist" = true ] && [ "$FORCE_RERUN" = false ]; then
    echo "✅ 下游分析结果已存在:"
    for file in "${analysis_files[@]}"; do
        if [ -f "$file" ]; then
            echo "   - $(basename "$file"): $(du -h "$file" | cut -f1)"
        fi
    done
    echo "   - 跳过下游分析步骤"
else
    if [ "$FORCE_RERUN" = true ]; then
        echo "🔄 强制重新运行：开始下游分析..."
    else
        echo "📊 开始运行下游分析脚本..."
    fi
    
    echo "  - 聚类分析"
    Rscript code/sc3cluster.R
    
    echo "  - 可视化分析"
    Rscript code/Classification\ visualization.R
    
    echo "  - 差异表达分析"
    Rscript code/differential\ gene\ expression.R
    
    echo "  - 通路富集分析"
    Rscript code/kegg.R
fi

echo "----------------------------------------------------------"
echo "🎉 全部流程完成！"
echo ""
echo "📊 执行统计："
echo "  - 数据解压: ✅ 完成"
echo "  - 路径检测: ✅ 完成"
echo "  - 数据预处理: $(if [ -f "$preprocessed_file" ] && [ "$FORCE_RERUN" = false ]; then echo "⏭️  跳过（已存在）"; else echo "✅ 完成"; fi)"
echo "  - 数据生成: $(if [ -f "$generated_file" ] && [ "$FORCE_RERUN" = false ]; then echo "⏭️  跳过（已存在）"; else echo "✅ 完成"; fi)"
echo "  - 下游分析: $(if [ "$all_analysis_exist" = true ] && [ "$FORCE_RERUN" = false ]; then echo "⏭️  跳过（已存在）"; else echo "✅ 完成"; fi)"
echo ""
echo "📁 结果目录："
echo "  - 预处理数据: FD1000/"
echo "  - 生成数据: output/"
echo "  - 分析结果: output/"
echo ""
echo "⏰ 当前时间: $(date)"
echo "==========================================================" 