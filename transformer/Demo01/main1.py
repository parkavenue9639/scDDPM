import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import os

from b import UNet1D, DiffusionModule
from data import CSVDataBalancer


class GeneExpressionDataset(Dataset):
    def __init__(self, dataframe):
        """
        dataframe: 包含特征和标签的DataFrame
        """
        # 提取特征（假设前1000列是基因表达数据）
        self.features = dataframe.iloc[:, 1:-1].values.astype(np.float32)  # 跳过第一列名称
        self.labels = dataframe.iloc[:, -1].values

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        return {
            "expression": torch.tensor(self.features[idx]),
            "label": torch.tensor(self.labels[idx], dtype=torch.long)
        }


def train_and_generate(model, diffusion, dataloader, device, num_epochs=100):
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)

    model.train()
    for epoch in range(num_epochs):
        total_loss = 0
        for batch in dataloader:
            optimizer.zero_grad()

            x = batch["expression"].to(device)
            labels = batch["label"].to(device)

            loss = diffusion(x, labels)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        print(f"Epoch {epoch + 1}/{num_epochs} Loss: {total_loss / len(dataloader):.4f}")


def post_process(generated, original_df):
    # 转换为numpy并处理负值
    generated_np = generated.cpu().numpy()
    threshold = np.abs(original_df.iloc[:, 1:-1].values.min())
    generated_np[generated_np < -threshold] = 0

    # 数据缩放
    original_mean = original_df.iloc[:, 1:-1].mean().values
    original_std = original_df.iloc[:, 1:-1].std().values
    generated_mean = generated_np.mean(axis=0)
    generated_std = generated_np.std(axis=0)

    # 应用缩放
    scaled_data = (generated_np - generated_mean) * (original_std / (generated_std + 1e-8)) + original_mean
    return scaled_data


def main():
    # 配置参数
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    num_epochs = 100
    batch_size = 128
    num_samples_per_class = 1000  # 每个类别生成样本数

    # 初始化数据平衡器
    csv_path = "../数据/ALZHEIMER/AD01103/预处理数据/FD1000/AD01103PreProLabel1000.csv"
    data_balancer = CSVDataBalancer(csv_path)
    data_balancer.balance_data()
    balanced_dict = data_balancer.get_balanced_dict()

    # 准备输出文件
    output_path = "generate.csv"
    header_written = False

    # 遍历每个类别
    for label, df in balanced_dict.items():
        print(f"\nProcessing class {label}...")

        # 创建数据集和数据加载器
        dataset = GeneExpressionDataset(df)
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        # 初始化模型
        num_classes = len(balanced_dict)
        model = UNet1D(num_classes=num_classes).to(device)
        diffusion = DiffusionModule(model).to(device)

        # 训练模型
        train_and_generate(model, diffusion, dataloader, device, num_epochs)

        # 生成数据
        with torch.no_grad():
            labels_tensor = torch.full((num_samples_per_class,), label, device=device)
            generated = diffusion.sample(
                cell_types=labels_tensor,
                num_samples=num_samples_per_class,
                device=device
            )

        processed_data = post_process(generated, df)

        gen_df = pd.DataFrame(processed_data, columns=df.columns[1:-1])

        start_idx = len(pd.read_csv(output_path)) if os.path.exists(output_path) else 0
        gen_df.insert(0, "Cell", [f"cell{i + start_idx + 1}" for i in range(num_samples_per_class)])
        gen_df["label"] = label

        gen_df.to_csv(output_path, mode="a", header=not header_written, index=False)
        header_written = True

        print(f"Generated {num_samples_per_class} samples for class {label}")


