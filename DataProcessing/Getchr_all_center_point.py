import numpy as np
import os
import re

def extract_center_points_from_npy(folder, output_folder):
    # 确保输出目录存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 遍历文件夹中的 .npy 文件
    for filename in os.listdir(folder):
        if filename.endswith(".npy") and not filename.endswith("_tmp.npy"):
            filepath = os.path.join(folder, filename)

            # 用正则提取染色体号
            match = re.search(r"(chr\d+)", filename)
            if not match:
                print(f"文件 {filename} 未匹配到染色体号，跳过")
                continue
            chrom = match.group(1)

            # 读取 numpy 文件并替换 NaN
            data = np.load(filepath)
            data = np.nan_to_num(data, nan=0.0)

            # 提取前两列
            positions = data[:, :2]

            # 输出 txt 文件路径
            output_file = os.path.join(output_folder, filename.replace(".npy", "_centerpoints.txt"))

            # 保存
            with open(output_file, "w") as f:
                for pos in positions:
                    f.write(f"{chrom}\t{int(pos[0])}\t{int(pos[1])}\n")

            # 校验行数
            print(f"{filename}: 原始行数 = {data.shape[0]}, 输出行数 = {positions.shape[0]}")
            print(f"已保存: {output_file}\n")


# ===== 使用示例 =====
input_folder = "/path/of/chr_all_sample"   # 存放 .npy 文件的文件夹
output_folder = "/path/of/chr_all_centerpoint_txt" # 输出目录

extract_center_points_from_npy(input_folder, output_folder)
