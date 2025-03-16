import pandas as pd
import numpy as np
import os

class CSVDataBalancer:
    def __init__(self, csv_path):
        """
        初始化 CSV 读取路径，并读取数据。
        :param csv_path: str, CSV 文件路径
        """
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"文件 {csv_path} 不存在，请检查路径！")

        self.csv_path = csv_path
        self.df = pd.read_csv(csv_path)

        # 提取列名
        self.column_names = self.df.columns.tolist()
        self.label_col = self.column_names[-1]  # 假设最后一列是 label
        self.name_col = self.column_names[0]  # 假设第一列是名字

        # 按 label 分类存入字典
        self.label_groups = {label: self.df[self.df[self.label_col] == label].copy()
                             for label in self.df[self.label_col].unique()}

        # 计算最大类别样本数
        self.max_count = max(len(group) for group in self.label_groups.values())

        # 存储平衡后的数据
        self.balanced_data = {}

    def balance_data(self):
        """
        对数据进行平衡，每个类别扩充到相同数量。
        """
        for label, group in self.label_groups.items():
            n_samples = len(group)
            if n_samples < self.max_count:
                extra_samples = group.sample(self.max_count - n_samples, replace=True, random_state=42)
                balanced_group = pd.concat([group, extra_samples], ignore_index=True)
            else:
                balanced_group = group.copy()

            self.balanced_data[label] = balanced_group  # 以 label 作为 key 存储 DataFrame

    def get_balanced_data(self, label=None):
        """
        获取平衡后的数据。
        :param label: str 或 None, 指定类别时返回该类别的 DataFrame，否则返回完整 DataFrame
        :return: DataFrame
        """
        if label:
            return self.balanced_data.get(label, None)
        else:
            return pd.concat(self.balanced_data.values(), ignore_index=True)

    def get_balanced_dict(self):
        """
        返回一个字典，其中 key 是类别，value 是该类别的 DataFrame。
        :return: dict
        """
        return self.balanced_data

    def save_balanced_data(self, save_path="balanced_data.csv"):
        """
        保存平衡后的数据到 CSV 文件。
        :param save_path: str, 保存的文件路径
        """
        balanced_df = self.get_balanced_data()
        balanced_df.to_csv(save_path, index=False)
        print(f"平衡后的数据已保存至 {save_path}")

    def show_balanced_info(self):
        """
        打印每个类别的样本数量。
        """
        for label, df in self.balanced_data.items():
            print(f"类别 {label} 数据大小: {df.shape}")
            print(df.head())  # 仅显示前 5 行数据

