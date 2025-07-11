#!/bin/bash

set -e

# ================== 参数解析 ==================
FORCE_RERUN_ALL=false
FORCE_RERUN_TRAIN=false
FORCE_RERUN_GENERATE=false
FORCE_RERUN_KEGG=false
FORCE_RERUN_PCA=false
DATASET_NAME=""  # 默认为空，表示处理所有数据集
PROCESS_ALL_DATASETS=true  # 默认处理所有数据集

# 针对特定数据集的强制重执行选项
FORCE_TRAIN_FOR=""
FORCE_GENERATE_FOR=""
FORCE_KEGG_FOR=""
FORCE_PCA_FOR=""
FORCE_ALL_FOR=""
SKIP_DATASETS=""

for arg in "$@"; do
    case $arg in
        --force|-f)
            FORCE_RERUN_ALL=true
            ;;
        --force-train|-ft)
            FORCE_RERUN_TRAIN=true
            ;;
        --force-generate|-fg)
            FORCE_RERUN_GENERATE=true
            ;;
        --force-kegg|-fk)
            FORCE_RERUN_KEGG=true
            ;;
        --force-pca|-fp)
            FORCE_RERUN_PCA=true
            ;;
        --dataset=*)
            DATASET_NAME="${arg#*=}"
            PROCESS_ALL_DATASETS=false
            ;;
        -d|--dataset)
            shift
            DATASET_NAME="$1"
            PROCESS_ALL_DATASETS=false
            ;;
        --force-train-for=*)
            FORCE_TRAIN_FOR="${arg#*=}"
            PROCESS_ALL_DATASETS=false
            ;;
        --force-generate-for=*)
            FORCE_GENERATE_FOR="${arg#*=}"
            PROCESS_ALL_DATASETS=false
            ;;
        --force-kegg-for=*)
            FORCE_KEGG_FOR="${arg#*=}"
            PROCESS_ALL_DATASETS=false
            ;;
        --force-pca-for=*)
            FORCE_PCA_FOR="${arg#*=}"
            PROCESS_ALL_DATASETS=false
            ;;
        --force-all-for=*)
            FORCE_ALL_FOR="${arg#*=}"
            PROCESS_ALL_DATASETS=false
            ;;
        --skip-dataset=*)
            SKIP_DATASETS="${SKIP_DATASETS} ${arg#*=}"
            ;;
    esac
done

# 所有支持的数据集
ALL_DATASETS="AD00202 AD00203 AD00204 AD00401 AD01103"

# 验证数据集名称的函数
validate_dataset() {
    local dataset=$1
    case $dataset in
        AD00202|AD00203|AD00204|AD00401|AD01103)
            return 0
            ;;
        *)
            echo "❌ 无效的数据集名称: $dataset"
            echo "   支持的数据集: AD00202, AD00203, AD00204, AD00401, AD01103"
            exit 1
            ;;
    esac
}

# 检查是否在跳过列表中
is_dataset_skipped() {
    local dataset=$1
    echo "$SKIP_DATASETS" | grep -q "$dataset"
}

# 处理参数逻辑
if [ "$PROCESS_ALL_DATASETS" = true ]; then
    if [ -n "$DATASET_NAME" ]; then
        validate_dataset "$DATASET_NAME"
        DATASETS_TO_PROCESS="$DATASET_NAME"
        echo "✅ 使用指定数据集: $DATASET_NAME"
    else
        DATASETS_TO_PROCESS="$ALL_DATASETS"
        echo "🔄 批量处理所有数据集: $ALL_DATASETS"
        if [ -n "$SKIP_DATASETS" ]; then
            echo "⏭️  跳过数据集:$SKIP_DATASETS"
        fi
    fi
else
    # 处理针对特定数据集的操作
    TARGET_DATASETS=""
    if [ -n "$FORCE_TRAIN_FOR" ]; then
        validate_dataset "$FORCE_TRAIN_FOR"
        TARGET_DATASETS="$TARGET_DATASETS $FORCE_TRAIN_FOR"
    fi
    if [ -n "$FORCE_GENERATE_FOR" ]; then
        validate_dataset "$FORCE_GENERATE_FOR"
        TARGET_DATASETS="$TARGET_DATASETS $FORCE_GENERATE_FOR"
    fi
    if [ -n "$FORCE_KEGG_FOR" ]; then
        validate_dataset "$FORCE_KEGG_FOR"
        TARGET_DATASETS="$TARGET_DATASETS $FORCE_KEGG_FOR"
    fi
    if [ -n "$FORCE_PCA_FOR" ]; then
        validate_dataset "$FORCE_PCA_FOR"
        TARGET_DATASETS="$TARGET_DATASETS $FORCE_PCA_FOR"
    fi
    if [ -n "$FORCE_ALL_FOR" ]; then
        validate_dataset "$FORCE_ALL_FOR"
        TARGET_DATASETS="$TARGET_DATASETS $FORCE_ALL_FOR"
    fi
    if [ -n "$DATASET_NAME" ]; then
        validate_dataset "$DATASET_NAME"
        TARGET_DATASETS="$TARGET_DATASETS $DATASET_NAME"
    fi
    
    # 去重并设置处理列表
    DATASETS_TO_PROCESS=$(echo "$TARGET_DATASETS" | tr ' ' '\n' | sort -u | tr '\n' ' ')
    echo "🎯 针对性处理数据集: $DATASETS_TO_PROCESS"
fi

echo "================= scDDPM 完整自动化脚本 ================="
echo "当前时间: $(date)"
echo "工作目录: $(pwd)"
echo "使用方法: $0 [选项]"
echo ""
echo "📊 数据集选项:"
echo "    无参数                   批量处理所有数据集 (默认)"
echo "    -d, --dataset DATASET   指定单个数据集"
echo "    --skip-dataset=DATASET  跳过指定数据集 (可多次使用)"
echo ""
echo "🔄 全局强制重运行选项:"
echo "    -f, --force            强制重新运行所有步骤"
echo "    --force-train, -ft     强制重新训练模型"
echo "    --force-generate, -fg  强制重新生成数据"
echo "    --force-kegg, -fk      强制重新运行KEGG分析"
echo "    --force-pca, -fp       强制重新运行PCA评估"
echo ""
echo "🎯 针对特定数据集的强制重运行选项:"
echo "    --force-train-for=DATASET     强制重新训练指定数据集"
echo "    --force-generate-for=DATASET  强制重新生成指定数据集"
echo "    --force-kegg-for=DATASET      强制重新运行指定数据集KEGG分析"
echo "    --force-pca-for=DATASET       强制重新运行指定数据集PCA评估"
echo "    --force-all-for=DATASET       强制重新运行指定数据集所有步骤"
echo ""
echo "💡 使用示例:"
echo "    $0                             # 处理所有数据集"
echo "    $0 -d AD00203                  # 只处理AD00203"
echo "    $0 --skip-dataset=AD00202      # 处理除AD00202外的所有数据集"
echo "    $0 --force-train-for=AD00203   # 只对AD00203重新训练"
echo "    $0 --force-kegg-for=AD00204    # 只对AD00204重新运行KEGG"
echo ""
echo "📚 支持的数据集: AD00202, AD00203, AD00204, AD00401, AD01103"
echo "----------------------------------------------------------"

# ================== 主处理函数 ==================
process_single_dataset() {
    local current_dataset=$1
    local dataset_force_train=$2
    local dataset_force_generate=$3
    local dataset_force_kegg=$4
    local dataset_force_pca=$5
    local dataset_force_all=$6
    
    echo ""
    echo "🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥"
    echo "              开始处理数据集: $current_dataset"
    echo "🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥"
    
    # 设置当前数据集的环境变量
    export DATASET_NAME=$current_dataset
    
    # 确定强制重运行选项
    local force_train_current=$FORCE_RERUN_TRAIN
    local force_generate_current=$FORCE_RERUN_GENERATE
    local force_kegg_current=$FORCE_RERUN_KEGG
    local force_pca_current=$FORCE_RERUN_PCA
    local force_all_current=$FORCE_RERUN_ALL
    
    # 应用针对特定数据集的强制选项
    if [ "$dataset_force_train" = true ] || [ "$dataset_force_all" = true ]; then
        force_train_current=true
    fi
    if [ "$dataset_force_generate" = true ] || [ "$dataset_force_all" = true ]; then
        force_generate_current=true
    fi
    if [ "$dataset_force_kegg" = true ] || [ "$dataset_force_all" = true ]; then
        force_kegg_current=true
    fi
    if [ "$dataset_force_pca" = true ] || [ "$dataset_force_all" = true ]; then
        force_pca_current=true
    fi
    if [ "$dataset_force_all" = true ]; then
        force_all_current=true
    fi

# ================== 智能数据解压模块 ==================
echo "🔍 [1/8] 检查并解压数据..."

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
echo "📁 [2/8] 创建输出目录..."
mkdir -p FD1000
mkdir -p output

# ================== 显示配置信息 ==================
echo "🔧 [3/8] 显示配置信息..."
python config.py

# ================== 路径检测 ==================
echo "🔍 [4/8] 检查路径和文件完整性..."
python test_paths.py

# ================== 数据预处理 ==================
echo "🧹 [5/8] 检查数据预处理..."

preprocessed_file="FD1000/${current_dataset}PreProLabel1000.csv"
if [ -f "$preprocessed_file" ] && [ "$force_all_current" = false ]; then
    echo "✅ 预处理数据已存在: $preprocessed_file"
    echo "   - 文件大小: $(du -h "$preprocessed_file" | cut -f1)"
    echo "   - 跳过预处理步骤"
else
    if [ "$force_all_current" = true ]; then
        echo "🔄 强制重新运行：开始预处理..."
    else
        echo "📊 预处理数据不存在，开始运行预处理..."
    fi
    Rscript code/Preprocess.R
    if [ ! -f "$preprocessed_file" ]; then
        echo "❌ 预处理输出文件不存在，流程中止！"
        exit 1
    fi
    echo "✅ 预处理完成，文件已生成: $preprocessed_file"
fi

# ================== 数据训练 ==================
echo "🧠 [6/8] 检查模型训练..."

model_dir="models"
model_pattern="${model_dir}/${current_dataset}*_best_*.pth"
model_count=$(ls $model_pattern 2>/dev/null | wc -l)

if { [ $model_count -gt 0 ] && [ "$force_all_current" = false ] && [ "$force_train_current" = false ]; }; then
    echo "✅ 最佳模型文件已存在: $model_count 个"
    echo "   - 跳过训练步骤"
else
    if [ "$force_all_current" = true ] || [ "$force_train_current" = true ]; then
        echo "🔄 强制重新运行：开始模型训练..."
        # 可选：删除旧模型
        rm -f $model_pattern
    else
        echo "🧠 最佳模型文件不存在，开始训练..."
    fi
    python code/train_model.py
    model_count=$(ls $model_pattern 2>/dev/null | wc -l)
    if [ $model_count -eq 0 ]; then
        echo "❌ 未生成最佳模型文件，流程中止！"
        exit 1
    fi
    echo "✅ 模型训练完成，生成 $model_count 个最佳模型文件"
fi

# ================== 数据生成 ==================
echo "🎨 [7/8] 检查数据生成..."

generated_file="output/${current_dataset}/generated_data/${current_dataset}_generated.csv"
if [ -f "$generated_file" ] && [ "$force_all_current" = false ] && [ "$force_generate_current" = false ]; then
    echo "✅ 生成数据已存在: $generated_file"
    echo "   - 文件大小: $(du -h "$generated_file" | cut -f1)"
    echo "   - 跳过生成步骤"
else
    if [ "$force_all_current" = true ] || [ "$force_generate_current" = true ]; then
        echo "🔄 强制重新运行：开始数据生成..."
        if [ -f "$generated_file" ]; then
            echo "🗑️  删除旧的生成数据文件: $generated_file"
            rm "$generated_file"
        fi
    else
        echo "🎨 生成数据不存在，开始运行generate_data.py..."
    fi
    python code/generate_data.py
    if [ ! -f "$generated_file" ]; then
        echo "❌ 生成数据文件不存在，流程中止！"
        exit 1
    fi
    echo "✅ 数据生成完成，文件已生成: $generated_file"
fi

# ================== 下游分析 ==================
echo "📊 [7/9] 检查下游分析..."

# 检查下游分析结果文件（不包括KEGG分析）
analysis_files=(
    "output/clustering_results.csv"
    "output/visualization_plots.pdf"
    "output/differential_expression.csv"
)

# 检查是否所有分析文件都存在
all_analysis_exist=true
for file in "${analysis_files[@]}"; do
    if [ ! -f "$file" ]; then
        all_analysis_exist=false
        break
    fi
done

if [ "$all_analysis_exist" = true ] && [ "$force_all_current" = false ]; then
    echo "✅ 下游分析结果已存在:"
    for file in "${analysis_files[@]}"; do
        if [ -f "$file" ]; then
            echo "   - $(basename "$file"): $(du -h "$file" | cut -f1)"
        fi
    done
    echo "   - 跳过下游分析步骤"
else
    if [ "$force_all_current" = true ]; then
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
fi

# ================== KEGG富集分析 ==================
echo "🧬 [8/9] 智能KEGG富集分析..."

# 检查必需的输入文件
if [ ! -f "$preprocessed_file" ]; then
    echo "❌ 预处理数据文件不存在: $preprocessed_file"
    echo "   请先运行数据预处理步骤"
    exit 1
fi

if [ ! -f "$generated_file" ]; then
    echo "❌ 生成数据文件不存在: $generated_file"
    echo "   请先运行数据生成步骤"
    exit 1
fi

# 检查KEGG结果文件
kegg_real_file="output/${current_dataset}/kegg_analysis/kegg_real_results.rds"
kegg_gen_file="output/${current_dataset}/kegg_analysis/kegg_generated_results.rds"
kegg_basic_results=(
    "output/${current_dataset}/kegg_analysis/KEGG_comparison_and_venn.pdf"
    "output/${current_dataset}/kegg_analysis/common_kegg_pathways.csv"
)
kegg_detailed_results=(
    "output/${current_dataset}/kegg_analysis/detailed_kegg_comparison.pdf"
    "output/${current_dataset}/kegg_analysis/pathway_enrichment_comparison.csv"
    "output/${current_dataset}/kegg_analysis/kegg_summary_statistics.csv"
)

# 检查基础KEGG分析是否完成
basic_kegg_complete=true
for file in "${kegg_basic_results[@]}"; do
    if [ ! -f "$file" ]; then
        basic_kegg_complete=false
        break
    fi
done

# 检查详细KEGG分析是否完成
detailed_kegg_complete=true
for file in "${kegg_detailed_results[@]}"; do
    if [ ! -f "$file" ]; then
        detailed_kegg_complete=false
        break
    fi
done

# 智能KEGG分析执行逻辑
if [ "$force_all_current" = true ] || [ "$force_kegg_current" = true ]; then
    echo "🔄 强制重新运行KEGG分析..."
    
    # 清理旧的KEGG结果
    if [ "$force_kegg_current" = true ]; then
        echo "🗑️  清理旧的KEGG分析结果..."
        rm -f "$kegg_real_file" "$kegg_gen_file"
        for file in "${kegg_basic_results[@]}" "${kegg_detailed_results[@]}"; do
            [ -f "$file" ] && rm -f "$file"
        done
    fi
    
    echo "  - 运行基础KEGG分析..."
    Rscript code/kegg.R
    
    echo "  - 运行详细KEGG比较分析..."
    Rscript code/kegg_detailed_analysis.R
    
    echo "✅ KEGG富集分析完成（强制重新运行）"
    
elif [ -f "$kegg_real_file" ] && [ -f "$kegg_gen_file" ] && [ "$detailed_kegg_complete" = true ]; then
    echo "✅ 完整的KEGG分析结果已存在:"
    echo "   - KEGG结果文件: $(du -h "$kegg_real_file" | cut -f1), $(du -h "$kegg_gen_file" | cut -f1)"
    for file in "${kegg_basic_results[@]}" "${kegg_detailed_results[@]}"; do
        if [ -f "$file" ]; then
            echo "   - $(basename "$file"): $(du -h "$file" | cut -f1)"
        fi
    done
    echo "   - 跳过KEGG分析步骤"
    
elif [ -f "$kegg_real_file" ] && [ -f "$kegg_gen_file" ] && [ "$basic_kegg_complete" = true ]; then
    echo "✅ 基础KEGG分析已完成，运行详细分析..."
    echo "   - 从已保存的结果运行详细比较分析..."
    Rscript code/kegg_detailed_analysis.R
    echo "✅ 详细KEGG分析完成"
    
else
    echo "🧬 开始完整的KEGG富集分析流程..."
    
    # 检查KEGG结果状态
    echo "  - 检查KEGG分析状态..."
    Rscript code/load_kegg_results.R
    
    echo "  - 运行基础KEGG分析..."
    Rscript code/kegg.R
    
    echo "  - 运行详细KEGG比较分析..."
    Rscript code/kegg_detailed_analysis.R
    
    # 验证分析结果
    if [ ! -f "$kegg_real_file" ] || [ ! -f "$kegg_gen_file" ]; then
        echo "❌ KEGG分析未能生成结果文件，流程中止！"
        exit 1
    fi
    
    echo "✅ KEGG富集分析流程完成"
fi

# ================== PCA质量评估 ==================
echo "📊 [9/9] PCA质量评估..."

# 检查必需的输入文件
if [ ! -f "$preprocessed_file" ] || [ ! -f "$generated_file" ]; then
    echo "❌ 输入数据文件不存在，无法进行PCA质量评估"
    echo "   需要文件: $preprocessed_file, $generated_file"
    exit 1
fi

# 检查PCA评估结果文件
pca_assessment_file="output/${current_dataset}/quality_assessment/pca_quality_assessment.csv"
pca_visualization_files=(
    "output/${current_dataset}/visualization/real_data_pca.pdf"
    "output/${current_dataset}/visualization/generated_data_pca.pdf"
    "output/${current_dataset}/visualization/combined_data_pca.pdf"
    "output/${current_dataset}/visualization/real_data_pca_ggplot.pdf"
    "output/${current_dataset}/visualization/generated_data_pca_ggplot.pdf"
    "output/${current_dataset}/visualization/combined_data_pca_ggplot.pdf"
)

# 检查PCA评估是否完成
pca_assessment_complete=true
if [ ! -f "$pca_assessment_file" ]; then
    pca_assessment_complete=false
fi

# 检查基础可视化是否完成
basic_viz_complete=true
for file in "${pca_visualization_files[@]}"; do
    if [ ! -f "$file" ]; then
        basic_viz_complete=false
        break
    fi
done

# 智能PCA评估执行逻辑
if [ "$force_all_current" = true ] || [ "$force_pca_current" = true ]; then
    echo "🔄 强制重新运行PCA质量评估..."
    
    # 清理旧的PCA评估结果
    if [ "$force_pca_current" = true ]; then
        echo "🗑️  清理旧的PCA评估结果..."
        rm -f "$pca_assessment_file"
        for file in "${pca_visualization_files[@]}"; do
            [ -f "$file" ] && rm -f "$file"
        done
    fi
    
    echo "  - 运行PCA可视化分析..."
    Rscript code/Classification\ visualization.R
    
    echo "  - 运行PCA质量评估..."
    Rscript code/pca_quality_assessment.R
    
    echo "✅ PCA质量评估完成（强制重新运行）"
    
elif [ -f "$pca_assessment_file" ] && [ "$basic_viz_complete" = true ]; then
    echo "✅ PCA质量评估结果已存在:"
    echo "   - 评估文件: $(du -h "$pca_assessment_file" | cut -f1)"
    for file in "${pca_visualization_files[@]}"; do
        if [ -f "$file" ]; then
            echo "   - $(basename "$file"): $(du -h "$file" | cut -f1)"
        fi
    done
    echo "   - 跳过PCA评估步骤"
    
elif [ "$basic_viz_complete" = true ] && [ ! -f "$pca_assessment_file" ]; then
    echo "✅ PCA可视化已完成，运行质量评估..."
    echo "  - 运行PCA质量评估..."
    Rscript code/pca_quality_assessment.R
    echo "✅ PCA质量评估完成"
    
else
    echo "📊 开始完整的PCA分析和质量评估流程..."
    
    echo "  - 运行PCA可视化分析..."
    Rscript code/Classification\ visualization.R
    
    echo "  - 运行PCA质量评估..."
    Rscript code/pca_quality_assessment.R
    
    # 验证评估结果
    if [ ! -f "$pca_assessment_file" ]; then
        echo "❌ PCA质量评估未能生成结果文件，流程中止！"
        exit 1
    fi
    
    echo "✅ PCA分析和质量评估流程完成"
fi

    echo ""
    echo "✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅"
    echo "           数据集 $current_dataset 处理完成！"
    echo "✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅"
    echo ""
}

# ================== 主执行逻辑 ==================

# 主循环处理数据集
TOTAL_DATASETS=0
SUCCESSFUL_DATASETS=0
FAILED_DATASETS=""

for current_dataset in $DATASETS_TO_PROCESS; do
    # 检查是否跳过当前数据集
    if is_dataset_skipped "$current_dataset"; then
        echo "⏭️  跳过数据集: $current_dataset"
        continue
    fi
    
    TOTAL_DATASETS=$((TOTAL_DATASETS + 1))
    
    # 确定针对当前数据集的强制选项
    dataset_force_train=false
    dataset_force_generate=false
    dataset_force_kegg=false
    dataset_force_pca=false
    dataset_force_all=false
    
    # 检查针对特定数据集的强制选项
    if [ -n "$FORCE_TRAIN_FOR" ] && [ "$FORCE_TRAIN_FOR" = "$current_dataset" ]; then
        dataset_force_train=true
    fi
    if [ -n "$FORCE_GENERATE_FOR" ] && [ "$FORCE_GENERATE_FOR" = "$current_dataset" ]; then
        dataset_force_generate=true
    fi
    if [ -n "$FORCE_KEGG_FOR" ] && [ "$FORCE_KEGG_FOR" = "$current_dataset" ]; then
        dataset_force_kegg=true
    fi
    if [ -n "$FORCE_PCA_FOR" ] && [ "$FORCE_PCA_FOR" = "$current_dataset" ]; then
        dataset_force_pca=true
    fi
    if [ -n "$FORCE_ALL_FOR" ] && [ "$FORCE_ALL_FOR" = "$current_dataset" ]; then
        dataset_force_all=true
    fi
    
    # 处理当前数据集
    if process_single_dataset "$current_dataset" "$dataset_force_train" "$dataset_force_generate" "$dataset_force_kegg" "$dataset_force_pca" "$dataset_force_all"; then
        SUCCESSFUL_DATASETS=$((SUCCESSFUL_DATASETS + 1))
    else
        FAILED_DATASETS="$FAILED_DATASETS $current_dataset"
    fi
done

echo "----------------------------------------------------------"
echo "🎉 批量处理流程完成！"
echo ""
echo "📊 总体执行统计："
echo "  - 总处理数据集数: $TOTAL_DATASETS"
echo "  - 成功处理数据集数: $SUCCESSFUL_DATASETS"
if [ -n "$FAILED_DATASETS" ]; then
    echo "  - 失败数据集:$FAILED_DATASETS"
fi
echo "  - 数据解压: ✅ 完成"
echo "  - 路径检测: ✅ 完成"
echo "  - 各数据集详细处理状态请查看上方日志"
echo ""
echo "📁 结果目录结构："
echo "  - 预处理数据: FD1000/"
echo "  - 训练模型: models/"
echo "  - 输出根目录: output/"
echo "  - 按数据集组织: output/[DATASET_NAME]/"
echo "    ├── generated_data/           # 生成的单细胞数据"
echo "    ├── visualization/            # 可视化图表"
echo "    ├── clustering/               # 聚类分析结果"
echo "    ├── differential_expression/  # 差异表达分析"
echo "    ├── kegg_analysis/            # KEGG富集分析"
echo "    ├── quality_assessment/       # 数据质量评估"
echo "    └── reports/                  # 分析报告"
echo ""
echo "⏰ 当前时间: $(date)"
echo "==========================================================" 