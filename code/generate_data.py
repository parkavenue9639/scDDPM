import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import os
import time
import glob
from train_model import UNet1D, DiffusionModule, CSVDataBalancer
from tqdm import trange

# --------------------- Post-processing of Generated Data ---------------------
def post_process(generated, original_df):
    """
    Post-process generated data:
    1. Zero out values smaller than the absolute minimum of each column from original data.
    2. Align the mean and standard deviation between generated and real data.
    3. Ensure no negative values (biological plausibility).

    :param generated: Generated data (Tensor)
    :param original_df: Original data (DataFrame)
    :return: Adjusted NumPy array
    """
    gen_np = generated.cpu().numpy()

    # Remove label column, retain only numeric values
    orig_values = original_df.iloc[:, 1:-1].values

    # Ensure dimension match
    assert gen_np.shape[1] == orig_values.shape[1], \
        f"Column mismatch: generated {gen_np.shape[1]}, original {orig_values.shape[1]}"

    # Zero out values below the column-wise minimum absolute value
    threshold = np.abs(orig_values).min(axis=0)
    for col in range(gen_np.shape[1]):
        gen_np[:, col][gen_np[:, col] < threshold[col]] = 0

    # Mean and std correction
    M1, S1 = orig_values.mean(axis=0), orig_values.std(axis=0)
    M2, S2 = gen_np.mean(axis=0), gen_np.std(axis=0)
    S2[S2 < 1e-8] = 1e-8  # avoid division by zero

    adjusted_np = np.clip((gen_np - M2) * (S1 / S2) + M1, 0, None)
    return adjusted_np


# --------------------- Load Trained Model ---------------------
def load_model(model_path, device):
    """
    Load a trained model from checkpoint
    
    :param model_path: Path to the model checkpoint
    :param device: Device to load the model on
    :return: model, diffusion, label, num_classes
    """
    checkpoint = torch.load(model_path, map_location=device)
    
    # Extract model parameters
    num_classes = checkpoint['num_classes']
    label = checkpoint['label']
    
    # Initialize model
    model = UNet1D(num_classes=num_classes).to(device)
    diffusion = DiffusionModule(model).to(device)
    
    # Load state dicts
    model.load_state_dict(checkpoint['model_state_dict'])
    diffusion.load_state_dict(checkpoint['diffusion_state_dict'])
    
    # Set to evaluation mode
    model.eval()
    diffusion.eval()
    
    return model, diffusion, label, num_classes


# --------------------- Generate Data for Single Class ---------------------
def generate_for_class(model_path, csv_path, output_path, device, samples_per_batch=700, num_batches=6):
    """
    Generate data for a single cell type
    
    :param model_path: Path to the trained model
    :param csv_path: Path to the original CSV data
    :param output_path: Path to save generated data
    :param device: Device to run generation on
    :param samples_per_batch: Number of samples to generate per batch
    :param num_batches: Number of batches to generate
    """
    print(f"🎨 开始生成数据...")
    print(f"   - 模型文件: {model_path}")
    print(f"   - 原始数据: {csv_path}")
    print(f"   - 输出文件: {output_path}")
    
    # Load model
    model, diffusion, label, num_classes = load_model(model_path, device)
    print(f"   - 细胞类型: {label}")
    print(f"   - 类别数量: {num_classes}")
    
    # Load original data for post-processing
    original_df = pd.read_csv(csv_path)
    original_columns = original_df.columns.tolist()
    
    # Get label mapping
    balancer = CSVDataBalancer(csv_path)
    unique_labels = np.unique(balancer.df[balancer.label_col])
    label_to_num = {label: i for i, label in enumerate(unique_labels)}
    label_groups = {label: balancer.df[balancer.label_col] == label for label in unique_labels}
    
    # Find the label index
    unique_labels = list(label_groups.keys())
    label_to_num = {label: i for i, label in enumerate(unique_labels)}
    label_idx = label_to_num[label]
    
    # Get the data for the specific cell type
    cell_type_data = label_groups[label]
    
    print(f"   - 标签索引: {label_idx}")
    
    with torch.no_grad():
        generated_list = []
        total_samples = samples_per_batch * num_batches
        
        for batch_idx in trange(num_batches, desc="生成批次", unit="batch"):
            print(f"   - 生成批次 {batch_idx + 1}/{num_batches}...")
            labels = torch.full((samples_per_batch,), label_idx, device=device, dtype=torch.long)
            batch_generated = diffusion.sample(labels, samples_per_batch, device)
            generated_list.append(batch_generated)
            print(f"     ✅ 批次 {batch_idx + 1} 完成")
        
        print(f"🔗 合并生成数据...")
        generated = torch.cat(generated_list, dim=0)
        print(f"   - 生成数据形状: {generated.shape}")

        print(f"🔧 后处理数据...")
        processed = post_process(generated, cell_type_data)
        print(f"   - 处理后形状: {processed.shape}")

        print(f"💾 保存数据...")
        gen_df = pd.DataFrame(processed, columns=original_columns[1:-1])
        gen_df.insert(0, "Cell", [f"cell_{label.replace(' ', '_')}_{i}" for i in range(total_samples)])
        gen_df["label"] = label

        gen_df.to_csv(output_path, mode="a", header=not os.path.exists(output_path), index=False)
        print(f"   ✅ 数据已保存到: {output_path}")


# --------------------- Load Unified Model ---------------------
def load_unified_model(model_path, device):
    """
    Load a unified trained model from checkpoint
    
    :param model_path: Path to the unified model checkpoint
    :param device: Device to load the model on
    :return: model, diffusion, all_labels, num_classes, label_counts
    """
    checkpoint = torch.load(model_path, map_location=device)
    
    # Extract model parameters
    num_classes = checkpoint['num_classes']
    all_labels = checkpoint['all_labels']
    label_counts = checkpoint['label_counts']
    
    # Initialize model
    model = UNet1D(num_classes=num_classes).to(device)
    diffusion = DiffusionModule(model).to(device)
    
    # Load state dicts
    model.load_state_dict(checkpoint['model_state_dict'])
    diffusion.load_state_dict(checkpoint['diffusion_state_dict'])
    
    # Set to evaluation mode
    model.eval()
    diffusion.eval()
    
    return model, diffusion, all_labels, num_classes, label_counts


# --------------------- Generate Data for All Classes ---------------------
def generate_for_all_classes(model_path, csv_path, output_path, device, samples_per_batch=700, num_batches=6):
    """
    Generate data for all cell types using unified model
    
    :param model_path: Path to the unified trained model
    :param csv_path: Path to the original CSV data
    :param output_path: Path to save generated data
    :param device: Device to run generation on
    :param samples_per_batch: Number of samples to generate per batch
    :param num_batches: Number of batches to generate
    """
    print(f"🎨 开始使用统一模型生成数据...")
    print(f"   - 模型文件: {model_path}")
    print(f"   - 原始数据: {csv_path}")
    print(f"   - 输出文件: {output_path}")
    
    # Load unified model
    model, diffusion, all_labels, num_classes, label_counts = load_unified_model(model_path, device)
    print(f"   - 类别数量: {num_classes}")
    print(f"   - 所有类别: {all_labels}")
    
    # Load original data for post-processing
    original_df = pd.read_csv(csv_path)
    original_columns = original_df.columns.tolist()
    
    # Get label mapping
    balancer = CSVDataBalancer(csv_path)
    label_groups = {label: balancer.df[balancer.df[balancer.label_col] == label].copy() for label in unique_labels}
    
    # Create label to index mapping
    unique_labels = list(label_groups.keys())
    label_to_num = {label: i for i, label in enumerate(unique_labels)}
    
    print(f"   - 标签映射: {label_to_num}")
    
    # Generate data for each class
    total_classes = len(all_labels)
    
    for class_idx, label in enumerate(all_labels):
        print(f"\n{'='*50}")
        print(f"🧬 生成细胞类型 {class_idx + 1}/{total_classes}: {label}")
        print(f"   - 原始样本数: {label_counts.get(label, 0)}")
        print(f"   - 标签索引: {label_to_num[label]}")
        
        # Calculate generation count
        gen_count = label_counts.get(label, 0)
        if gen_count > 7000:
            gen_count = gen_count // 3  # Limit samples for large classes
        
        if gen_count > 0:
            with torch.no_grad():
                generated_list = []
                total_samples = samples_per_batch * num_batches
                
                for batch_idx in range(num_batches):
                    print(f"   - 生成批次 {batch_idx + 1}/{num_batches}...")
                    labels = torch.full((samples_per_batch,), label_to_num[label], device=device, dtype=torch.long)
                    batch_generated = diffusion.sample(labels, samples_per_batch, device)
                    generated_list.append(batch_generated)
                    print(f"     ✅ 批次 {batch_idx + 1} 完成")
                
                print(f"🔗 合并生成数据...")
                generated = torch.cat(generated_list, dim=0)
                print(f"   - 生成数据形状: {generated.shape}")

                print(f"🔧 后处理数据...")
                processed = post_process(generated, label_groups[label])
                print(f"   - 处理后形状: {processed.shape}")

                print(f"💾 保存数据...")
                gen_df = pd.DataFrame(processed, columns=original_columns[1:-1])
                gen_df.insert(0, "Cell", [f"cell_{label.replace(' ', '_')}_{i}" for i in range(total_samples)])
                gen_df["label"] = label

                gen_df.to_csv(output_path, mode="a", header=not os.path.exists(output_path), index=False)
                print(f"   ✅ 数据已保存到: {output_path}")
        else:
            print(f"   ⚠️  跳过生成（样本数为0）")


# --------------------- Main Generation Function ---------------------
def main():
    print("🎨 启动 scDDPM 数据生成...")
    print("=" * 60)
    
    # 设备检测
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print(f"🎯 使用 GPU: {torch.cuda.get_device_name()}")
    elif torch.backends.mps.is_available():
        device = torch.device("mps")
        print(f"🎯 使用 Apple Silicon GPU (MPS)")
    else:
        device = torch.device("cpu")
        print(f"🎯 使用 CPU")
    
    # ================== Configuration ==================
    dataset_name = "AD01103"  # 可以修改为: AD00202, AD00203, AD00204, AD00401, AD01103
    csv_path = f"FD1000/{dataset_name}PreProLabel1000.csv"
    model_dir = "models"
    output_path = f"output/{dataset_name}_generated.csv"
    
    print(f"📁 数据集: {dataset_name}")
    print(f"📂 输入文件: {csv_path}")
    print(f"📂 模型目录: {model_dir}")
    print(f"📂 输出文件: {output_path}")
    print("-" * 60)
    
    # 检查输入文件是否存在
    if not os.path.exists(csv_path):
        print(f"❌ 输入文件不存在: {csv_path}")
        print("请先运行 Preprocess.R 生成预处理数据")
        return
    
    # 检查模型目录是否存在
    if not os.path.exists(model_dir):
        print(f"❌ 模型目录不存在: {model_dir}")
        print("请先运行 train_model.py 训练模型")
        return
    
    # 创建输出目录
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # 移除之前的输出文件
    if os.path.exists(output_path):
        os.remove(output_path)
        print(f"🗑️  已删除之前的输出文件: {output_path}")
    
    # 获取所有最佳模型文件
    best_model_pattern = os.path.join(model_dir, f"{dataset_name}_*_best_*.pth")
    best_model_files = glob.glob(best_model_pattern)
    
    if not best_model_files:
        print(f"❌ 未找到最佳模型文件: {best_model_pattern}")
        print("请先运行 train_model.py 训练模型")
        return
    
    print(f"📊 找到 {len(best_model_files)} 个最佳模型文件")
    
    # 按细胞类型分组模型文件
    cell_type_models = {}
    for model_file in best_model_files:
        # 从文件名提取细胞类型
        filename = os.path.basename(model_file)
        # 格式: AD01103_{cell_type}_best_{number}.pth
        parts = filename.replace('.pth', '').split('_')
        if len(parts) >= 4 and parts[-2] == 'best':
            cell_type = '_'.join(parts[1:-2])  # 提取细胞类型部分
            if cell_type not in cell_type_models:
                cell_type_models[cell_type] = []
            cell_type_models[cell_type].append(model_file)
    
    print(f"📊 找到 {len(cell_type_models)} 个细胞类型的最佳模型")
    
    # 为每个细胞类型选择最佳模型并生成数据
    for cell_type, model_files in cell_type_models.items():
        print(f"\n{'='*60}")
        print(f"🧬 处理细胞类型: {cell_type}")
        print(f"   - 找到 {len(model_files)} 个模型文件")
        print(f"{'='*60}")
        
        # 选择损失最小的模型
        best_model = None
        best_loss = float('inf')
        
        for model_file in model_files:
            try:
                checkpoint = torch.load(model_file, map_location='cpu')
                loss = checkpoint['loss']
                if loss < best_loss:
                    best_loss = loss
                    best_model = model_file
            except Exception as e:
                print(f"⚠️  无法读取模型文件 {model_file}: {e}")
                continue
        
        if best_model:
            print(f"🎯 选择最佳模型: {os.path.basename(best_model)} (Loss: {best_loss:.4f})")
            
            try:
                generate_for_class(
                    model_path=best_model,
                    csv_path=csv_path,
                    output_path=output_path,
                    device=device,
                    samples_per_batch=700,
                    num_batches=6
                )
            except Exception as e:
                print(f"❌ 生成失败: {e}")
                continue
        else:
            print(f"❌ 无法找到有效的模型文件")
            continue
    
    print("\n" + "=" * 60)
    print("🎉 数据生成完成！")
    print(f"📁 生成数据位置: {output_path}")
    print("=" * 60)


# --------------------- Entry Point ---------------------
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

    print(f"⏰ 生成总耗时: {end_time - start_time:.2f} 秒") 