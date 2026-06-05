# 生成与“合并正负样本”格式完全一致的预测集（无标签版）
import numpy as np
import os
import re

def extract_chromosome_number(filename):

    match = re.search(r'chr(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

def process_chr_matrices(chrom, dnase_npy_path, hic_npy_path, center_txt_path):

    data_list = []
    # 1. 读中心点
    coords = np.loadtxt(center_txt_path, dtype=str)  # (M,2)

    # 2. 读矩阵
    dnase = np.load(dnase_npy_path)   # (M,21,21)
    hic   = np.load(hic_npy_path)[:, 2:444].reshape(-1, 21, 21)      # (M,21,21)

    M = len(coords)
    assert dnase.shape[0] == hic.shape[0] == M, '样本数不一致'

    for i in range(M):
        # 与合并脚本保持同一三元组结构
        data_list.append((coords[i], dnase[i], hic[i]))
    return data_list

def save_combined_data_noLabel(output_path, chromosomes,
                               dnase_dir, hic_dir, center_dir):

    all_data = []
    for chrom in chromosomes:
        chr_name = f'chr{chrom}'

        center_file = os.path.join(center_dir, f'{chr_name}_matrixsize21_centerpoints.txt')
        hic_file    = os.path.join(hic_dir,   f'{chr_name}_matrixsize21.npy')
        dnase_file  = os.path.join(dnase_dir, f'dnase_{chr_name}.npy')
        

        if not (os.path.exists(center_file) and os.path.exists(dnase_file) and os.path.exists(hic_file)):
            print(f'[Warn] 跳过 {chr_name} ：文件不全')
            continue

        chr_data = process_chr_matrices(chrom, dnase_file, hic_file, center_file)
        all_data.extend(chr_data)

    # 与合并脚本完全相同的保存方式
    np.savez(output_path,
             position_data=np.array([x[0] for x in all_data]),
             hic_data=np.array([x[2] for x in all_data]),
             dnase_data=np.array([x[1] for x in all_data]))
             
    print(f'预测集已生成 {output_path}  总样本 {len(all_data)}')
    # ===== 打印每个 key 的形状 =====
    with np.load(output_path) as f:
        for k in f.files:
            print(f"{k} {f[k].shape} {f[k].dtype}")

# ---------------- main ----------------
if __name__ == '__main__':
    chromosomes = [20, 21, 22]

    center_dir = '/path/of/chr_all_centerpoint_txt'
    dnase_dir  = '/path/of/DNase-All'
    hic_dir    = '/path/of/chr_all_sample'

    output_path = '/path/of/chr20-22_predict.npz'

    save_combined_data_noLabel(output_path, chromosomes,
                               dnase_dir, hic_dir, center_dir)