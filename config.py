# ================== scDDPM 配置文件 ==================
# 统一管理数据集配置，避免在多个文件中重复修改

# 当前要处理的数据集
# 可选值: "AD00202", "AD00203", "AD00204", "AD00401", "AD01103"
DATASET_NAME = "AD01103"

# 数据目录配置
DATA_DIR = "data"
OUTPUT_DIR = "output"
FD1000_DIR = "FD1000"

# 文件路径配置
def get_input_path(dataset_name=None):
    """获取输入文件路径"""
    if dataset_name is None:
        dataset_name = DATASET_NAME
    return f"{DATA_DIR}/{dataset_name}/{dataset_name}_expr.txt"

def get_label_path(dataset_name=None):
    """获取标签文件路径"""
    if dataset_name is None:
        dataset_name = DATASET_NAME
    return f"{DATA_DIR}/{dataset_name}/{dataset_name}_cell_label.txt"

def get_preprocessed_path(dataset_name=None):
    """获取预处理后的文件路径"""
    if dataset_name is None:
        dataset_name = DATASET_NAME
    return f"{FD1000_DIR}/{dataset_name}PreProLabel1000.csv"

def get_generated_path(dataset_name=None):
    """获取生成数据的输出路径"""
    if dataset_name is None:
        dataset_name = DATASET_NAME
    return f"{OUTPUT_DIR}/{dataset_name}_generated.csv"

# 模型配置
MODEL_CONFIG = {
    "timesteps": 400,
    "beta_start": 1e-4,
    "beta_end": 2e-2,
    "base_channels": 128,
    "time_emb_dim": 32,
    "learning_rate": 1e-4,
    "batch_size": 128,
    "epochs": 100
}

# 数据预处理配置
PREPROCESS_CONFIG = {
    "min_cells_percent": 0.05,  # 最小细胞百分比
    "n_features": 1000,         # 选择的基因数量
    "normalization_method": "LogNormalize",
    "scale_factor": 10000
}

# 打印当前配置
if __name__ == "__main__":
    print("🔧 scDDPM 配置信息:")
    print(f"数据集: {DATASET_NAME}")
    print(f"输入文件: {get_input_path()}")
    print(f"标签文件: {get_label_path()}")
    print(f"预处理文件: {get_preprocessed_path()}")
    print(f"生成文件: {get_generated_path()}")
    print(f"模型配置: {MODEL_CONFIG}")
    print(f"预处理配置: {PREPROCESS_CONFIG}") 