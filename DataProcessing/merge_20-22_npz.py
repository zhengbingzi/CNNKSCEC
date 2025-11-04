import numpy as np
import os

# 定义输入输出目录
input_dir = "/path/of/Dnase-Hic/"
output_dir = "/path/of/Predict/"

# 确保输出目录存在
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 获取所有生成的 .npz 文件列表
npz_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".npz")])

# 初始化合并的列表
merged_position_data = []
merged_dnase_data = []
merged_hic_data = []

# 遍历每个 .npz 文件并加载数据
for npz_file in npz_files:
    file_path = os.path.join(input_dir, npz_file)
    data = np.load(file_path)

    # 加载 position_data, dnase_data, hic_data
    position_data = data['position_data']
    dnase_data = data['dnase_data']
    hic_data = data['hic_data']

    # 将数据添加到合并列表
    merged_position_data.append(position_data)
    merged_dnase_data.append(dnase_data)
    merged_hic_data.append(hic_data)

    print(f"Loaded: {npz_file}")

# 合并所有的数据
merged_position_data = np.concatenate(merged_position_data, axis=0)
merged_dnase_data = np.concatenate(merged_dnase_data, axis=0)
merged_hic_data = np.concatenate(merged_hic_data, axis=0)

# 保存合并后的数据为一个新的 .npz 文件
output_file = os.path.join(output_dir, "Dnase_Hic_all_merged.npz")
np.savez(output_file, position_data=merged_position_data, dnase_data=merged_dnase_data, hic_data=merged_hic_data)

print(f"Saved: {output_file}")
