import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import os


class CSVDataBalancer:
    def __init__(self, csv_path):
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"文件 {csv_path} 不存在")
        self.df = pd.read_csv(csv_path)
        self.label_col = self.df.columns[-1]
        self.label_groups = {label: self.df[self.df[self.label_col] == label].copy()
                             for label in self.df[self.label_col].unique()}
        self.max_count = max(len(group) for group in self.label_groups.values())
        self.balanced_data = {}

    def balance_data(self):
        for label, group in self.label_groups.items():
            n_samples = len(group)
            if n_samples < self.max_count:
                extra_samples = group.sample(self.max_count - n_samples, replace=True, random_state=42)
                self.balanced_data[label] = pd.concat([group, extra_samples], ignore_index=True)
            else:
                self.balanced_data[label] = group.copy()

    def get_balanced_dict(self):
        return self.balanced_data


class GeneExpressionDataset(Dataset):
    def __init__(self, dataframe):
        self.features = dataframe.iloc[:, 1:-1].values.astype(np.float32)  # 确保float32
        self.labels = dataframe.iloc[:, -1].values.astype(np.int64)

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        return {
            "expression": torch.tensor(self.features[idx], dtype=torch.float32),
            "label": torch.tensor(self.labels[idx], dtype=torch.long)
        }


class TimeEmbedding(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.time_mlp = nn.Sequential(
            nn.Linear(1, dim),
            nn.GELU(),
            nn.Linear(dim, dim)
        )

    def forward(self, t):
        return self.time_mlp(t.float().unsqueeze(-1))  # 处理整数输入


class ConvBlock(nn.Module):
    def __init__(self, ch, time_emb_dim):
        super().__init__()
        self.time_mlp = TimeEmbedding(time_emb_dim)
        self.conv = nn.Sequential(
            nn.Conv1d(ch, ch, 3, padding=1),
            nn.BatchNorm1d(ch),
            nn.GELU(),
            nn.Conv1d(ch, ch, 3, padding=1),
            nn.BatchNorm1d(ch),
            nn.GELU()
        )
        self.condition = nn.Linear(time_emb_dim, ch)

    def forward(self, x, t):
        t_emb = self.time_mlp(t)
        scale = self.condition(t_emb).unsqueeze(-1)
        return self.conv(x) * (scale + 1)


class UNet1D(nn.Module):
    def __init__(self, in_channels=1, base_channels=128, time_emb_dim=32, num_classes=10):
        super().__init__()
        self.class_emb = nn.Embedding(num_classes, base_channels)

        # 初始投影
        self.initial = nn.Conv1d(in_channels, base_channels, kernel_size=1)

        # 编码器
        self.enc1 = ConvBlock(base_channels, time_emb_dim)
        self.pool1 = nn.MaxPool1d(2)
        self.enc2 = ConvBlock(base_channels, time_emb_dim)
        self.pool2 = nn.MaxPool1d(2)
        self.enc3 = ConvBlock(base_channels, time_emb_dim)
        self.pool3 = nn.MaxPool1d(2)

        # 瓶颈层
        self.bottleneck_norm = nn.GroupNorm(8, base_channels)
        self.bottleneck_relu = nn.ReLU()
        self.bottleneck_conv = ConvBlock(base_channels, time_emb_dim)

        # 解码器（修正通道数）
        self.up1 = nn.ConvTranspose1d(base_channels, base_channels, kernel_size=2, stride=2)
        self.dec1 = ConvBlock(base_channels * 2, time_emb_dim)  # 输入256

        self.up2 = nn.ConvTranspose1d(base_channels, base_channels, kernel_size=2, stride=2)
        self.dec2 = ConvBlock(base_channels * 2, time_emb_dim)  # 输入256

        self.up3 = nn.ConvTranspose1d(base_channels, base_channels, kernel_size=2, stride=2)
        self.dec3 = ConvBlock(base_channels * 2, time_emb_dim)  # 输入256

        self.final = nn.Conv1d(base_channels, 1, kernel_size=1)

    def forward(self, x, t, cell_type):
        x = x.unsqueeze(1)  # [B,1,1000]
        cell_emb = self.class_emb(cell_type).unsqueeze(-1)  # [B,128,1]

        # 编码器
        e1 = self.enc1(self.initial(x), t)  # [B,128,1000]
        p1 = self.pool1(e1)  # [B,128,500]

        e2 = self.enc2(p1, t)  # [B,128,500]
        p2 = self.pool2(e2)  # [B,128,250]

        e3 = self.enc3(p2, t)  # [B,128,250]
        p3 = self.pool3(e3)  # [B,128,125]

        # 瓶颈层
        b = self.bottleneck_norm(p3)
        b = self.bottleneck_relu(b)
        b = self.bottleneck_conv(b, t)  # [B,128,125]
        b = b * (cell_emb + 1)

        # 解码器
        d1 = self.up1(b)  # [B,128,250]
        d1 = torch.cat([d1, e3], dim=1)  # [B,256,250]
        d1 = self.dec1(d1, t)  # [B,128,250]

        d2 = self.up2(d1)  # [B,128,500]
        d2 = torch.cat([d2, e2], dim=1)  # [B,256,500]
        d2 = self.dec2(d2, t)  # [B,128,500]

        d3 = self.up3(d2)  # [B,128,1000]
        d3 = torch.cat([d3, e1], dim=1)  # [B,256,1000]
        d3 = self.dec3(d3, t)  # [B,128,1000]

        return self.final(d3).squeeze(1)  # [B,1000]


class DiffusionModule(nn.Module):
    def __init__(self, model, timesteps=400, beta_start=1e-4, beta_end=2e-2):
        super().__init__()
        self.model = model
        self.timesteps = timesteps

        # 注册缓冲区
        betas = torch.linspace(beta_start, beta_end, timesteps, dtype=torch.float32)
        alphas = 1. - betas
        alphas_cumprod = torch.cumprod(alphas, 0)

        self.register_buffer('betas', betas)
        self.register_buffer('alphas', alphas)
        self.register_buffer('alphas_cumprod', alphas_cumprod)
        self.register_buffer('sqrt_alphas_cumprod', torch.sqrt(alphas_cumprod))
        self.register_buffer('sqrt_one_minus_alphas_cumprod', torch.sqrt(1. - alphas_cumprod))

    def forward(self, x, cell_types):
        batch_size = x.size(0)
        t = torch.randint(0, self.timesteps, (batch_size,), device=x.device).long()

        noise = torch.randn_like(x, dtype=torch.float32)
        sqrt_alpha = self.sqrt_alphas_cumprod[t].view(-1, 1).to(x.dtype)
        sqrt_one_minus = self.sqrt_one_minus_alphas_cumprod[t].view(-1, 1).to(x.dtype)

        x_noisy = sqrt_alpha * x + sqrt_one_minus * noise
        pred_noise = self.model(x_noisy, t, cell_types)

        return F.mse_loss(noise, pred_noise)

    @torch.no_grad()
    def sample(self, cell_types, num_samples, device):
        x = torch.randn((num_samples, 1000), device=device, dtype=torch.float32)

        for t in reversed(range(self.timesteps)):
            ts = torch.full((num_samples,), t, device=device, dtype=torch.long)
            pred_noise = self.model(x, ts, cell_types)

            alpha = self.alphas[t]
            alpha_cumprod = self.alphas_cumprod[t]
            beta = self.betas[t]

            if t > 0:
                noise = torch.randn_like(x)
            else:
                noise = 0

            x = (x - beta * pred_noise / torch.sqrt(1 - alpha_cumprod)) / torch.sqrt(alpha)
            x += torch.sqrt(beta) * noise

        return x


def train_model(model, diffusion, dataloader, device, epochs=100):
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    model.train()

    for epoch in range(epochs):
        total_loss = 0
        for batch in dataloader:
            optimizer.zero_grad()
            x = batch["expression"].to(device)
            labels = batch["label"].to(device)

            loss = diffusion(x, labels)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        print(f"Epoch {epoch + 1}/{epochs} Loss: {total_loss / len(dataloader):.4f}")


def post_process(generated, original_df):
    gen_np = generated.cpu().numpy()
    threshold = original_df.iloc[:, 1:-1].abs().min().min()
    gen_np[gen_np < -threshold] = 0

    # 数据缩放
    orig_mean = original_df.iloc[:, 1:-1].mean().values
    orig_std = original_df.iloc[:, 1:-1].std().values
    gen_mean = gen_np.mean(0)
    gen_std = gen_np.std(0)

    scaled = (gen_np - gen_mean) * (orig_std / (gen_std + 1e-8)) + orig_mean
    return np.clip(scaled, 0, None)


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    csv_path = "AD01103.csv"
    output_path = "generated.csv"

    # 数据平衡
    balancer = CSVDataBalancer(csv_path)
    balancer.balance_data()
    balanced_dict = balancer.get_balanced_dict()

    # 准备输出文件
    if os.path.exists(output_path):
        os.remove(output_path)
    original_columns = pd.read_csv(csv_path).columns.tolist()

    # 逐类别处理
    for label, df in balanced_dict.items():
        print(f"\nProcessing class {label}...")

        # 数据加载
        dataset = GeneExpressionDataset(df)
        dataloader = DataLoader(dataset, batch_size=128, shuffle=True)

        # 初始化模型
        model = UNet1D(num_classes=len(balanced_dict)).to(device)
        diffusion = DiffusionModule(model).to(device)

        # 训练
        train_model(model, diffusion, dataloader, device, epochs=100)

        # 生成
        with torch.no_grad():
            labels = torch.full((1000,), label, device=device, dtype=torch.long)
            generated = diffusion.sample(labels, 1000, device)

        # 后处理
        processed = post_process(generated, df)

        # 构建DataFrame
        gen_df = pd.DataFrame(processed, columns=original_columns[1:-1])
        gen_df.insert(0, "Cell", [f"cell_{label}_{i}" for i in range(1000)])
        gen_df["label"] = label

        # 保存
        gen_df.to_csv(output_path, mode="a", header=not os.path.exists(output_path), index=False)

