import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import QuantileTransformer
import os
import matplotlib.pyplot as plt
import time

# --------------------- Data Loading Module ---------------------
class CSVDataBalancer:
    def __init__(self, csv_path):
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"File {csv_path} does not exist")
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


# --------------------- Dataset Definition ---------------------
class GeneExpressionDataset(Dataset):
    def __init__(self, dataframe):
        self.features = dataframe.iloc[:, 1:-1].values.astype(np.float32)  # ensure float32 type
        
        # Handle string labels by converting to numeric
        label_col = dataframe.iloc[:, -1].values
        if label_col.dtype == 'object':  # string labels
            unique_labels = np.unique(label_col)
            label_to_num = {label: i for i, label in enumerate(unique_labels)}
            self.labels = np.array([label_to_num[label] for label in label_col], dtype=np.int64)
            print(f"✅ 标签映射: {label_to_num}")
        else:
            self.labels = label_col.astype(np.int64) - 1  # label starts from 0

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        return {
            "expression": torch.tensor(self.features[idx], dtype=torch.float32),
            "label": torch.tensor(self.labels[idx], dtype=torch.long)
        }


# --------------------- Model Definition ---------------------
class TimeEmbedding(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.time_mlp = nn.Sequential(
            nn.Linear(1, dim),
            nn.GELU(),
            nn.Linear(dim, dim)
        )

    def forward(self, t):
        return self.time_mlp(t.float().unsqueeze(-1))  # convert integer to float


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

        # Initial projection
        self.initial = nn.Conv1d(in_channels, base_channels, kernel_size=1)

        # Encoder
        self.enc1 = ConvBlock(base_channels, time_emb_dim)
        self.pool1 = nn.MaxPool1d(2)
        self.enc2 = ConvBlock(base_channels, time_emb_dim)
        self.pool2 = nn.MaxPool1d(2)
        self.enc3 = ConvBlock(base_channels, time_emb_dim)
        self.pool3 = nn.MaxPool1d(2)

        # Bottleneck
        self.bottleneck_norm = nn.GroupNorm(8, base_channels)
        self.bottleneck_relu = nn.ReLU()
        self.bottleneck_conv = ConvBlock(base_channels, time_emb_dim)

        # Decoder
        self.up1 = nn.ConvTranspose1d(base_channels, base_channels, kernel_size=2, stride=2)
        self.dec1 = ConvBlock(base_channels, time_emb_dim)

        self.up2 = nn.ConvTranspose1d(base_channels, base_channels, kernel_size=2, stride=2)
        self.dec2 = ConvBlock(base_channels, time_emb_dim)

        self.up3 = nn.ConvTranspose1d(base_channels, base_channels, kernel_size=2, stride=2)
        self.dec3 = ConvBlock(base_channels, time_emb_dim)

        self.final = nn.Conv1d(base_channels, 1, kernel_size=1)

    def forward(self, x, t, cell_type):
        x = x.unsqueeze(1)  # [B,1,1000]
        cell_emb = self.class_emb(cell_type).unsqueeze(-1)  # [B,128,1]

        # Encoder
        e1 = self.enc1(self.initial(x), t)     # [B,128,1000]
        p1 = self.pool1(e1)                    # [B,128,500]

        e2 = self.enc2(p1, t)                  # [B,128,500]
        p2 = self.pool2(e2)                    # [B,128,250]

        e3 = self.enc3(p2, t)                  # [B,128,250]
        p3 = self.pool3(e3)                    # [B,128,125]

        # Bottleneck
        b = self.bottleneck_norm(p3)
        b = self.bottleneck_relu(b)
        b = self.bottleneck_conv(b, t)         # [B,128,125]
        b = b * (cell_emb + 1)

        # Decoder
        d1 = self.up1(b)                       # [B,128,250]
        d1 = d1 + e3
        d1 = self.dec1(d1, t)

        d2 = self.up2(d1)                      # [B,128,500]
        d2 = d2 + e2
        d2 = self.dec2(d2, t)

        d3 = self.up3(d2)                      # [B,128,1000]
        d3 = d3 + e1
        d3 = self.dec3(d3, t)

        return self.final(d3).squeeze(1)       # [B,1000]


class DiffusionModule(nn.Module):
    def __init__(self, model, timesteps=400, beta_start=1e-4, beta_end=2e-2):
        super().__init__()
        self.model = model
        self.timesteps = timesteps

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


# --------------------- Training Function ---------------------
def train_model(model, diffusion, dataloader, device, epochs=100):
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    model.train()
    loss_history = []
    
    print(f"🚀 开始训练模型...")
    print(f"   - 设备: {device}")
    print(f"   - 总轮次: {epochs}")
    print(f"   - 批次数量: {len(dataloader)}")
    print(f"   - 学习率: 1e-4")
    print("=" * 50)

    for epoch in range(epochs):
        total_loss = 0
        batch_count = 0
        
        # 每10个epoch显示一次详细进度
        if epoch % 10 == 0:
            print(f"📊 Epoch {epoch + 1}/{epochs} - 开始训练...")
        
        for batch_idx, batch in enumerate(dataloader):
            optimizer.zero_grad()
            x = batch["expression"].to(device)
            labels = batch["label"].to(device)

            loss = diffusion(x, labels)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            batch_count += 1
            
            # 每50个batch显示一次进度
            if batch_idx % 50 == 0 and epoch % 10 == 0:
                current_loss = loss.item()
                progress = (batch_idx + 1) / len(dataloader) * 100
                print(f"   Batch {batch_idx + 1}/{len(dataloader)} ({progress:.1f}%) - Loss: {current_loss:.4f}")

        avg_loss = total_loss / len(dataloader)
        loss_history.append(avg_loss)
        
        # 每10个epoch显示一次平均损失
        if epoch % 10 == 0:
            print(f"✅ Epoch {epoch + 1}/{epochs} 完成 - 平均损失: {avg_loss:.4f}")
            print("-" * 30)
        else:
            print(f"Epoch {epoch + 1}/{epochs} Loss: {avg_loss:.4f}")
    
    print("🎉 训练完成！")
    print(f"最终损失: {loss_history[-1]:.4f}")
    print("=" * 50)
    
    return loss_history


# --------------------- Main Training Function ---------------------
def main():
    print("🚀 启动 scDDPM 模型训练...")
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
    dataset_name = os.environ.get('DATASET_NAME', 'AD01103')  # 从环境变量读取数据集名称，默认AD01103
    csv_path = f"FD1000/{dataset_name}PreProLabel1000.csv"
    model_save_dir = "models"
    
    print(f"📁 数据集: {dataset_name}")
    print(f"📂 输入文件: {csv_path}")
    print(f"📂 模型保存目录: {model_save_dir}")
    print("-" * 60)
    
    # 检查输入文件是否存在
    if not os.path.exists(csv_path):
        print(f"❌ 输入文件不存在: {csv_path}")
        print("请先运行 Preprocess.R 生成预处理数据")
        return

    # 创建模型保存目录
    os.makedirs(model_save_dir, exist_ok=True)

    # Balance the dataset
    print(f"⚖️  数据平衡处理...")
    balancer = CSVDataBalancer(csv_path)
    print(f"   - 原始数据形状: {balancer.df.shape}")
    print(f"   - 标签列: {balancer.label_col}")
    
    balancer.balance_data()
    balanced_dict = balancer.get_balanced_dict()
    print(f"   - 平衡后类别数: {len(balanced_dict)}")

    # Get total number of classes
    all_labels = set(balancer.df[balancer.label_col].unique())
    num_classes = len(all_labels)
    print(f"   - 总类别数: {num_classes}")

    # Count real sample number per class
    real_counts = balancer.df[balancer.label_col].value_counts().to_dict()
    print(f"   - 各类别样本数:")
    for label, count in real_counts.items():
        print(f"     * {label}: {count}")
    print("-" * 60)

    # Train per class
    total_classes = len(balanced_dict)
    current_class = 0
    
    for label, df in balanced_dict.items():
        current_class += 1
        gen_count = real_counts.get(label, 0)
        if gen_count > 7000:
            gen_count = gen_count // 3  # Limit samples for large classes
        
        print(f"\n{'='*60}")
        print(f"🧬 处理细胞类型 {current_class}/{total_classes}: {label}")
        print(f"   - 目标生成数量: {gen_count} 样本")
        print(f"   - 原始数据形状: {df.shape}")
        print(f"   - 设备: {device}")
        print(f"{'='*60}")

        if gen_count > 0:
            print(f"📊 创建数据集...")
            dataset = GeneExpressionDataset(df)
            dataloader = DataLoader(dataset, batch_size=128, shuffle=True)
            print(f"   - 数据集大小: {len(dataset)}")
            print(f"   - 批次数量: {len(dataloader)}")

            print(f"🏗️  初始化模型...")
            model = UNet1D(num_classes=num_classes).to(device)
            diffusion = DiffusionModule(model).to(device)
            print(f"   - 模型参数数量: {sum(p.numel() for p in model.parameters()):,}")
            print(f"   - 可训练参数: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}")

            # 训练模型并保存最佳检查点
            print(f"🎯 开始训练模型...")
            print("=" * 50)
            
            # 保存最佳模型的变量
            best_losses = []  # 存储最佳损失值
            best_models = []  # 存储最佳模型状态
            
            optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
            model.train()
            epochs = 100
            
            print(f"🚀 开始训练模型...")
            print(f"   - 设备: {device}")
            print(f"   - 总轮次: {epochs}")
            print(f"   - 批次数量: {len(dataloader)}")
            print(f"   - 学习率: 1e-4")
            print(f"   - 保存最佳前3个模型")
            print("=" * 50)

            for epoch in range(epochs):
                total_loss = 0
                batch_count = 0
                
                # 每10个epoch显示一次详细进度
                if epoch % 10 == 0:
                    print(f"📊 Epoch {epoch + 1}/{epochs} - 开始训练...")
                
                for batch_idx, batch in enumerate(dataloader):
                    optimizer.zero_grad()
                    x = batch["expression"].to(device)
                    labels = batch["label"].to(device)

                    loss = diffusion(x, labels)
                    loss.backward()
                    optimizer.step()

                    total_loss += loss.item()
                    batch_count += 1
                    
                    # 每50个batch显示一次进度
                    if batch_idx % 50 == 0 and epoch % 10 == 0:
                        current_loss = loss.item()
                        progress = (batch_idx + 1) / len(dataloader) * 100
                        print(f"   Batch {batch_idx + 1}/{len(dataloader)} ({progress:.1f}%) - Loss: {current_loss:.4f}")

                avg_loss = total_loss / len(dataloader)
                
                # 每10个epoch显示一次平均损失
                if epoch % 10 == 0:
                    print(f"✅ Epoch {epoch + 1}/{epochs} 完成 - 平均损失: {avg_loss:.4f}")
                    print("-" * 30)
                else:
                    print(f"Epoch {epoch + 1}/{epochs} Loss: {avg_loss:.4f}")
                
                # 保存最佳模型逻辑
                if len(best_losses) < 3:
                    # 前3个epoch直接保存
                    best_losses.append(avg_loss)
                    best_models.append({
                        'epoch': epoch + 1,
                        'loss': avg_loss,
                        'model_state_dict': model.state_dict().copy(),
                        'diffusion_state_dict': diffusion.state_dict().copy()
                    })
                    print(f"💾 保存第{len(best_models)}个最佳模型 (Epoch {epoch + 1}, Loss: {avg_loss:.4f})")
                else:
                    # 检查是否比最差的模型更好
                    worst_loss = max(best_losses)
                    if avg_loss < worst_loss:
                        # 找到最差模型的索引
                        worst_idx = best_losses.index(worst_loss)
                        # 替换最差模型
                        best_losses[worst_idx] = avg_loss
                        best_models[worst_idx] = {
                            'epoch': epoch + 1,
                            'loss': avg_loss,
                            'model_state_dict': model.state_dict().copy(),
                            'diffusion_state_dict': diffusion.state_dict().copy()
                        }
                        print(f"💾 更新最佳模型 (Epoch {epoch + 1}, Loss: {avg_loss:.4f})")
            
            print("🎉 当前类别训练完成！")
            print(f"最终损失: {avg_loss:.4f}")
            print("=" * 50)
            
            # 保存当前类别的最佳模型
            print(f"💾 保存 {label} 的最佳模型...")
            for i, (best_model, best_loss) in enumerate(zip(best_models, best_losses)):
                model_save_path = os.path.join(model_save_dir, f"{dataset_name}_{label.replace(' ', '_')}_best_{i+1}.pth")
                torch.save({
                    'model_state_dict': best_model['model_state_dict'],
                    'diffusion_state_dict': best_model['diffusion_state_dict'],
                    'epoch': best_model['epoch'],
                    'loss': best_model['loss'],
                    'label': label,
                    'num_classes': num_classes,
                    'epochs': 100
                }, model_save_path)
                print(f"   - 最佳模型 {i+1}: Epoch {best_model['epoch']}, Loss: {best_model['loss']:.4f}")
                print(f"     📁 保存到: {model_save_path}")

    print("\n" + "=" * 60)
    print("🎉 所有模型训练完成！")
    print("📊 每个类别都已训练并保存最佳模型")
    print(f"📁 模型文件位置: {model_save_dir}/")
    print("=" * 60)


# --------------------- Entry Point ---------------------
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()

    print(f"⏰ 训练总耗时: {end_time - start_time:.2f} 秒") 