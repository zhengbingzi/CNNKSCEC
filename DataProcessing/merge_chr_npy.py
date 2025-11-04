import os
import argparse
import numpy as np
from tqdm import tqdm
import re

def get_prefix(fname):
    # 从文件名中提取不带编号的前缀，例如：KR_5kb_matrix_chr1_negative_0.npy → KR_5kb_matrix_chr1_negative
    return "_".join(fname.split("_")[:-1])

def group_files_by_prefix(folder):
    groups = {}
    for fname in os.listdir(folder):
        if fname.endswith(".npy") and "chr" in fname:
            prefix = get_prefix(fname)
            groups.setdefault(prefix, []).append(os.path.join(folder, fname))
    return groups

def merge_npy_files(file_list):
    file_list.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]))
    return np.concatenate([np.load(f) for f in tqdm(file_list, desc="合并中")], axis=0)

def main(folder, output):
    os.makedirs(output, exist_ok=True)
    file_groups = group_files_by_prefix(folder)

    for prefix, files in file_groups.items():
        merged = merge_npy_files(files)
        save_path = os.path.join(output, f"{prefix}.npy")
        np.save(save_path, merged)
        print(f"保存: {save_path}  shape: {merged.shape}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", required=True, help="输入 .npy 文件所在目录")
    parser.add_argument("-o", "--output", default="merged_output", help="合并后输出目录")
    args = parser.parse_args()
    main(args.folder, args.output)
