#合并Dnase+Hic  ->  ([21×21],[21×21],lable).npz
#positive
import os
import numpy as np

dnase_dir = "/path/of/DNase-ALL/"
hic_dir = "/path/of/chr_all_sample/"
output_dir = "/path/of/Predict/"
# 检查输出目录是否存在，不存在则创建
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 获取所有的 DNase 和 Hi-C 文件列表
dnase_files = sorted([f for f in os.listdir(dnase_dir) if f.endswith(".npy")])
hic_files = sorted([f for f in os.listdir(hic_dir) if f.endswith(".npy")])

# 匹配文件并处理
for dnase_file in dnase_files:
    # 提取染色体号
    chrom = dnase_file.split("_")[1].replace(".npy","")
    
    # 找到对应染色体号的 Hi-C 文件
    hic_file = f"{chrom}_matrixsize21.npy"
    
    if hic_file in hic_files:
        # 读取 DNase 和 Hi-C 数据
        dnase_data = np.load(os.path.join(dnase_dir, dnase_file))  # 形状为 (n_samples, 21, 21)
        hic_data = np.load(os.path.join(hic_dir, hic_file))  # 形状为 (n_samples, 442)
        
        # 将 Hi-C 数据的前441列 reshape 为 (n_samples, 21, 21) 矩阵，最后一列为标签
        # hic_matrices = hic_data[:, :441].reshape(-1, 21, 21)
        # labels = hic_data[:, 441]
        hic_matrices = hic_data[:, 2:].reshape(-1, 21, 21)
        # labels = hic_data[:, -1]
        # 检查行数是否一致
        if dnase_data.shape[0] == hic_matrices.shape[0]:
            # 将 DNase 和 Hi-C 数据组合为字典
            data_dict = {
                'dnase_data': dnase_data,
                'hic_data': hic_matrices,
                # 'labels': labels
            }
            
            # 保存为 .npz 文件
            output_file = os.path.join(output_dir, f"Dnase_Hic_{chrom}.npz")
            np.savez(output_file, **data_dict)
            
            # 打印数据形状和第一行数据
            # print(f"Chromosome: {chrom}")
            # print(f"DNase data shape: {dnase_data.shape}")
            # print(f"Hi-C data shape: {hic_matrices.shape}")
            # print(f"Labels shape: {labels.shape}")
            # print(f"First row of DNase data:\n{dnase_data[0]}")
            # print(f"First row of Hi-C data:\n{hic_matrices[0]}")
            # print(f"First label: {labels[0]}")
            
            print(f"Saved: {output_file}")
        else:
            print(f"Row mismatch for chromosome {chrom}: DNase rows {dnase_data.shape[0]}, Hi-C rows {hic_matrices.shape[0]}")
    else:
        print(f"No matching Hi-C file found for {chrom}")