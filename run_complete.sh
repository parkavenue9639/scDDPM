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

# 日志相关选项
ENABLE_LOGGING=true
LOG_DIR="logs"
LOG_TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
MAIN_LOG_FILE=""
DISABLE_LOGGING=false

# 时间统计相关
TIMING_DATA_FILE=""
CURRENT_DATASET_TIMING=""

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
        --no-log)
            DISABLE_LOGGING=true
            ;;
        --log-dir=*)
            LOG_DIR="${arg#*=}"
            ;;
        --help|-h)
            # 这里只设置标志，实际的帮助信息在后面显示
            ;;
    esac
done

# ================== 日志系统初始化 ==================
# 初始化日志系统
init_logging() {
    if [ "$DISABLE_LOGGING" = true ]; then
        echo "📝 日志记录已禁用"
        return 0
    fi
    
    # 创建日志目录
    mkdir -p "$LOG_DIR"
    
    # 设置主日志文件
    MAIN_LOG_FILE="$LOG_DIR/scDDPM_run_${LOG_TIMESTAMP}.log"
    
    # 清理旧日志文件（保留最近10个）
    find "$LOG_DIR" -name "scDDPM_run_*.log" -type f | sort -r | tail -n +11 | xargs rm -f 2>/dev/null || true
    
    echo "📝 启用日志记录: $MAIN_LOG_FILE"
    echo "=== scDDPM 处理日志 - $(date) ===" > "$MAIN_LOG_FILE"
}

# 日志记录函数
log_output() {
    if [ "$DISABLE_LOGGING" = false ] && [ -n "$MAIN_LOG_FILE" ]; then
        tee -a "$MAIN_LOG_FILE"
    else
        cat
    fi
}

# 数据集专用日志记录函数
log_dataset_output() {
    local dataset_name=$1
    if [ "$DISABLE_LOGGING" = false ] && [ -n "$LOG_DIR" ]; then
        local dataset_log="$LOG_DIR/dataset_${dataset_name}_${LOG_TIMESTAMP}.log"
        tee -a "$dataset_log" | log_output
    else
        log_output
    fi
}

# ================== 时间统计系统 ==================
# 初始化时间统计
init_timing_system() {
    if [ "$DISABLE_LOGGING" = false ] && [ -n "$LOG_DIR" ]; then
        TIMING_DATA_FILE="$LOG_DIR/timing_data_${LOG_TIMESTAMP}.tmp"
    else
        TIMING_DATA_FILE="/tmp/scDDPM_timing_$$"
    fi
    > "$TIMING_DATA_FILE"  # 清空文件
}

# 记录模块开始时间
start_module_timer() {
    local dataset=$1
    local module=$2
    local start_time=$(date +%s)
    
    echo "START|${dataset}|${module}|${start_time}" >> "$TIMING_DATA_FILE"
    
    {
        echo "⏱️  [$(date +"%H:%M:%S")] 开始模块: $module"
    } | log_dataset_output "$dataset"
}

# 记录模块结束时间并计算耗时
end_module_timer() {
    local dataset=$1
    local module=$2
    local end_time=$(date +%s)
    
    # 查找对应的开始时间
    local start_time=$(grep "START|${dataset}|${module}|" "$TIMING_DATA_FILE" | tail -1 | cut -d'|' -f4)
    
    if [ -n "$start_time" ]; then
        local duration=$((end_time - start_time))
        local hours=$((duration / 3600))
        local minutes=$(((duration % 3600) / 60))
        local seconds=$((duration % 60))
        
        local time_str=""
        if [ $hours -gt 0 ]; then
            time_str="${hours}小时${minutes}分${seconds}秒"
        elif [ $minutes -gt 0 ]; then
            time_str="${minutes}分${seconds}秒"
        else
            time_str="${seconds}秒"
        fi
        
        # 记录完成时间
        echo "END|${dataset}|${module}|${end_time}|${time_str}" >> "$TIMING_DATA_FILE"
        
        {
            echo "✅ [$(date +"%H:%M:%S")] 完成模块: $module (耗时: $time_str)"
        } | log_dataset_output "$dataset"
    fi
}

# 输出数据集时间统计报告
output_timing_report() {
    local dataset=$1
    
    {
        echo ""
        echo "📊 数据集 $dataset 各模块耗时统计："
        echo "----------------------------------------"
        
        # 查找该数据集的所有完成记录
        local has_data=false
        while IFS='|' read -r type ds module end_time time_str; do
            if [ "$type" = "END" ] && [ "$ds" = "$dataset" ] && [ "$module" != "数据集处理总时间" ]; then
                printf "  %-20s %s\n" "$module:" "$time_str"
                has_data=true
            fi
        done < "$TIMING_DATA_FILE"
        
        if [ "$has_data" = false ]; then
            echo "  无时间统计数据"
        fi
        echo "----------------------------------------"
        echo ""
    } | log_dataset_output "$dataset"
}

# 输出总体时间统计报告
output_global_timing_report() {
    {
        echo ""
        echo "📊📊📊 总体各数据集处理时间统计 📊📊📊"
        echo "=============================================="
        
        for dataset in $DATASETS_TO_PROCESS; do
            # 跳过被跳过的数据集
            if is_dataset_skipped "$dataset"; then
                continue
            fi
            
            echo ""
            echo "📋 数据集: $dataset"
            echo "--------------------"
            
            local has_data=false
            local total_seconds=0
            
            while IFS='|' read -r type ds module end_time time_str; do
                if [ "$type" = "END" ] && [ "$ds" = "$dataset" ] && [ "$module" != "数据集处理总时间" ]; then
                    printf "  %-20s %s\n" "$module:" "$time_str"
                    has_data=true
                    
                    # 计算总秒数
                    local seconds=$(echo "$time_str" | grep -o '[0-9]*秒' | sed 's/秒//')
                    local minutes=$(echo "$time_str" | grep -o '[0-9]*分' | sed 's/分//')
                    local hours=$(echo "$time_str" | grep -o '[0-9]*小时' | sed 's/小时//')
                    
                    [ -n "$seconds" ] && total_seconds=$((total_seconds + seconds))
                    [ -n "$minutes" ] && total_seconds=$((total_seconds + minutes * 60))
                    [ -n "$hours" ] && total_seconds=$((total_seconds + hours * 3600))
                fi
            done < "$TIMING_DATA_FILE"
            
            if [ "$has_data" = false ]; then
                echo "  无时间统计数据"
            else
                echo "  --------------------"
                local total_hours=$((total_seconds / 3600))
                local total_minutes=$(((total_seconds % 3600) / 60))
                local remaining_seconds=$((total_seconds % 60))
                
                if [ $total_hours -gt 0 ]; then
                    echo "  数据集总耗时: ${total_hours}小时${total_minutes}分${remaining_seconds}秒"
                elif [ $total_minutes -gt 0 ]; then
                    echo "  数据集总耗时: ${total_minutes}分${remaining_seconds}秒"
                else
                    echo "  数据集总耗时: ${remaining_seconds}秒"
                fi
            fi
        done
        
        echo ""
        echo "=============================================="
        echo ""
    } | log_output
}

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

# 初始化日志系统
init_logging

# 初始化时间统计系统
init_timing_system

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
echo "📝 日志记录选项:"
echo "    --no-log                      禁用日志记录"
echo "    --log-dir=DIR                 指定日志目录 (默认: logs)"
echo ""
echo "💡 使用示例:"
echo "    $0                             # 处理所有数据集"
echo "    $0 -d AD00203                  # 只处理AD00203"
echo "    $0 --skip-dataset=AD00202      # 处理除AD00202外的所有数据集"
echo "    $0 --force-train-for=AD00203   # 只对AD00203重新训练"
echo "    $0 --force-kegg-for=AD00204    # 只对AD00204重新运行KEGG"
echo "    $0 --no-log                    # 禁用日志记录"
echo "    $0 --log-dir=my_logs           # 使用自定义日志目录"
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
    
    # 记录数据集处理开始时间
    start_module_timer "$current_dataset" "数据集处理总时间"
    
    {
        echo ""
        echo "🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥"
        echo "              开始处理数据集: $current_dataset"
        echo "🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥🔥"
        echo "⏰ 开始时间: $(date)"
    
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
start_module_timer "$current_dataset" "数据解压"
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
end_module_timer "$current_dataset" "数据解压"

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
start_module_timer "$current_dataset" "数据预处理"
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
end_module_timer "$current_dataset" "数据预处理"

# ================== 数据训练 ==================
start_module_timer "$current_dataset" "模型训练"
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
end_module_timer "$current_dataset" "模型训练"

# ================== 数据生成 ==================
start_module_timer "$current_dataset" "数据生成"
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
end_module_timer "$current_dataset" "数据生成"

# ================== 下游分析 ==================
start_module_timer "$current_dataset" "下游分析"
echo "📊 [7/9] 检查下游分析..."

# 检查下游分析结果文件（不包括KEGG分析）
analysis_files=(
    "output/${current_dataset}/clustering/clustering_results.csv"
    "output/${current_dataset}/visualization/visualization_plots.pdf"
    "output/${current_dataset}/differential_expression/differential_expression.csv"
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
end_module_timer "$current_dataset" "下游分析"

# ================== KEGG富集分析 ==================
start_module_timer "$current_dataset" "KEGG富集分析"
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
end_module_timer "$current_dataset" "KEGG富集分析"

# ================== PCA质量评估 ==================
start_module_timer "$current_dataset" "PCA质量评估"
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
end_module_timer "$current_dataset" "PCA质量评估"

# 输出当前数据集的时间统计报告
output_timing_report "$current_dataset"

# 记录数据集处理结束时间
end_module_timer "$current_dataset" "数据集处理总时间"

        echo ""
        echo "✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅"
        echo "           数据集 $current_dataset 处理完成！"
        echo "✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅✅"
        echo "⏰ 完成时间: $(date)"
        echo ""
    } | log_dataset_output "$current_dataset"
}

# ================== 主执行逻辑 ==================

# 主循环处理数据集
TOTAL_DATASETS=0
SUCCESSFUL_DATASETS=0
FAILED_DATASETS=""

{
    echo "🚀 开始批量处理流程 - $(date)"
    echo "📊 待处理数据集: $DATASETS_TO_PROCESS"
    if [ -n "$SKIP_DATASETS" ]; then
        echo "⏭️  跳过数据集:$SKIP_DATASETS"
    fi
    echo ""
} | log_output

for current_dataset in $DATASETS_TO_PROCESS; do
    # 检查是否跳过当前数据集
    if is_dataset_skipped "$current_dataset"; then
        {
            echo "⏭️  跳过数据集: $current_dataset"
        } | log_output
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
        {
            echo "✅ 数据集 $current_dataset 处理成功"
        } | log_output
    else
        FAILED_DATASETS="$FAILED_DATASETS $current_dataset"
        {
            echo "❌ 数据集 $current_dataset 处理失败"
        } | log_output
    fi
done

{
    echo "----------------------------------------------------------"
    echo "🎉 批量处理流程完成！"
    echo "⏰ 完成时间: $(date)"
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
    if [ "$DISABLE_LOGGING" = false ] && [ -n "$MAIN_LOG_FILE" ]; then
        echo "📝 详细日志已保存到: $MAIN_LOG_FILE"
        echo "📝 各数据集日志: $LOG_DIR/dataset_*_${LOG_TIMESTAMP}.log"
        echo ""
    fi
    echo "⏰ 总执行时间: 从 $LOG_TIMESTAMP 到 $(date +"%Y%m%d_%H%M%S")"
    echo "=========================================================="
} | log_output

# 输出总体时间统计报告
output_global_timing_report 