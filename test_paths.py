#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
路径测试脚本
验证所有文件路径是否正确配置
"""

import os
import sys
from config import *

def test_paths():
    """测试所有路径配置"""
    print("🔍 测试文件路径配置...")
    
    # 测试数据目录
    print(f"\n📂 数据目录检查:")
    print(f"  data目录: {'✅' if os.path.exists(DATA_DIR) else '❌'} {DATA_DIR}")
    print(f"  output目录: {'✅' if os.path.exists(OUTPUT_DIR) else '❌'} {OUTPUT_DIR}")
    print(f"  FD1000目录: {'✅' if os.path.exists(FD1000_DIR) else '❌'} {FD1000_DIR}")
    
    # 测试输入文件
    print(f"\n📄 输入文件检查:")
    input_path = get_input_path()
    label_path = get_label_path()
    
    print(f"  表达矩阵: {'✅' if os.path.exists(input_path) else '❌'} {input_path}")
    print(f"  细胞标签: {'✅' if os.path.exists(label_path) else '❌'} {label_path}")
    
    # 检查压缩文件
    compressed_expr = input_path.replace('.txt', '.txt.gz')
    compressed_label = label_path.replace('.txt', '.txt.gz')
    
    print(f"\n📦 压缩文件检查:")
    print(f"  压缩表达矩阵: {'✅' if os.path.exists(compressed_expr) else '❌'} {compressed_expr}")
    print(f"  压缩细胞标签: {'✅' if os.path.exists(compressed_label) else '❌'} {compressed_label}")
    
    # 测试输出文件路径
    print(f"\n📤 输出文件路径检查:")
    preprocessed_path = get_preprocessed_path()
    generated_path = get_generated_path()
    
    print(f"  预处理文件: {preprocessed_path}")
    print(f"  生成文件: {generated_path}")
    
    # 检查是否需要解压
    print(f"\n🔧 操作建议:")
    if not os.path.exists(input_path) and os.path.exists(compressed_expr):
        print("  ⚠️  需要解压数据文件，请运行: ./run.sh")
    elif not os.path.exists(input_path) and not os.path.exists(compressed_expr):
        print("  ❌ 数据文件不存在，请检查数据目录")
    else:
        print("  ✅ 数据文件已准备就绪")
    
    if not os.path.exists(preprocessed_path):
        print("  ⚠️  需要运行预处理，请运行: Rscript code/Preprocess.R")
    else:
        print("  ✅ 预处理文件已存在")
    
    return True

if __name__ == "__main__":
    test_paths() 